/*
 * equation.h
 *
 *  Created on: Mar 8, 2013
 *      Author: Marco Agnese
 */

#ifndef EQUATION_H_
#define EQUATION_H_

#include "modesmapping.h"
#include "plasmastate.h"
#include "utilities.h"

using namespace std;

template <class M, class B, class G>
class equation {

public:
    equation(M* mesh_ptr)
        : mesh(mesh_ptr)
        , basis()
        , integration()
        , map(0)
        , temp_order(0)
        , phi_values(0)
        , der_phi_values(0)
        , J_values(4)
        , value(0)
        , dirichlet_coeff(mesh->num_lines())
        , tollerance(0.0)
        , num_iteration(0)
        , projection_status(false)
        , null_value(0)
    {
        compute_basis();
        compute_dirichlet_coeff();
    }

    inline B* get_basis(void) { return &basis; }
    inline G* get_integration(void) { return &integration; }
    inline eq_type& get_equation_type(void) { return equation_type; }

    inline void set_ptr_mode_mapping(modes_mapping<M, B, plasma_state<M, B>>* ptr) { map = ptr; }

    void disp_basis(const int& rank = 0);
    void disp_der_basis(const int& rank = 0);
    void disp_dir_coeff(const int& rank = 0);

    inline double& dirichlet(const int& i, const int& l, const int& ele_pos) { return dirichlet_coeff[-mesh->n(l, ele_pos) - 1][i]; };

    inline int& get_iter(void) { return num_iteration; }
    inline void inc_iter(void) { ++num_iteration; }
    inline double& get_toll(void) { return tollerance; }

    inline bool& is_initial_guess_projected(void) { return projection_status; }
    inline void set_initial_guess_projected(void) { projection_status = true; }

    virtual ~equation() {};

protected:
    M* mesh;
    B basis;
    G integration;
    modes_mapping<M, B, plasma_state<M, B>>* map;

    int temp_order;

    vector<vector<vector<double>>> phi_values;
    vector<vector<vector<double>>> der_phi_values;
    vector<double> J_values;

    double value;

    vector<vector<double>> dirichlet_coeff;

    eq_type equation_type;

    double tollerance;
    int num_iteration;
    bool projection_status;
    double null_value;

    inline void set_equation_type(eq_type value) { equation_type = value; }

    inline int& get_order(void) { return temp_order; }
    inline void set_order(const int& value) { temp_order = value; }

    inline double& node(const int& pos) { return integration.get_nodes(temp_order)[pos]; }
    inline double& weight(const int& pos) { return integration.get_weights(temp_order)[pos]; }

    inline double x_physical(const int& k, const int& kk, const int& ele_pos) { return mesh->get_mapping()->Fx(node(k), node(kk), ele_pos); }
    inline double y_physical(const int& k, const int& kk, const int& ele_pos) { return mesh->get_mapping()->Fy(node(k), node(kk), ele_pos); }

    inline double& phi(const int& pos, const int& i) { return phi_values[i][integration.necessary_size(temp_order) - 2][pos]; }
    inline double& der_phi(const int& pos, const int& i) { return der_phi_values[i][integration.necessary_size(temp_order) - 2][pos]; }

    inline double& J(const int& idx) { return J_values[idx]; }

    void compute_basis(void);
    void compute_dirichlet_coeff(void);

private:
};

template <class M, class B, class G>
void equation<M, B, G>::compute_basis(void)
{

    unsigned old_size(phi_values.size());
    unsigned new_size(mesh->get_max_deg() + 1);

    if (new_size > old_size) {

        basis.set_order(mesh->get_max_deg());
        integration.set_order(mesh->get_max_deg());
        phi_values.resize(new_size);
        der_phi_values.resize(new_size);
        unsigned old_row_size(0);
        if (old_size > 0)
            old_row_size = phi_values[0].size();
        unsigned new_row_size(integration.necessary_size(mesh->get_max_deg()) - 1);

        for (unsigned k = 0; k != old_size; ++k) {
            phi_values[k].resize(new_row_size);
            der_phi_values[k].resize(new_row_size);
            for (unsigned j = old_row_size; j != new_row_size; ++j) {
                phi_values[k][j].resize(j + 2);
                der_phi_values[k][j].resize(j + 2);
                for (unsigned jj = 0; jj != (j + 2); ++jj) {
                    phi_values[k][j][jj] = basis.eval(integration.get_n(jj, j), k);
                    der_phi_values[k][j][jj] = basis.eval_der(integration.get_n(jj, j), k);
                }
            }
        }

        for (unsigned k = old_size; k != new_size; ++k) {
            phi_values[k].resize(new_row_size);
            der_phi_values[k].resize(new_row_size);
            for (unsigned j = 0; j != new_row_size; ++j) {
                phi_values[k][j].resize(j + 2);
                der_phi_values[k][j].resize(j + 2);
                for (unsigned jj = 0; jj != (j + 2); ++jj) {
                    phi_values[k][j][jj] = basis.eval(integration.get_n(jj, j), k);
                    der_phi_values[k][j][jj] = basis.eval_der(integration.get_n(jj, j), k);
                }
            }
        }
    }
}

template <class M, class B, class G>
void equation<M, B, G>::compute_dirichlet_coeff(void)
{

    set_order(mesh->get_max_deg());
    int num_integration_nodes(integration.necessary_size(get_order()));

    double x0(0);
    double x1(0);
    double y0(0);
    double y1(0);

    vector<double> A(static_cast<unsigned>((mesh->get_max_deg() + 1) * (mesh->get_max_deg() + 1)), 0);
    gsl_matrix_view A_gsl = gsl_matrix_view_array(&(A[0]), mesh->get_max_deg() + 1, mesh->get_max_deg() + 1);

    vector<double> b(static_cast<unsigned>(mesh->get_max_deg() + 1), 0);
    gsl_vector_view b_gsl = gsl_vector_view_array(&(b[0]), mesh->get_max_deg() + 1);

    gsl_vector* x_gsl_ptr = gsl_vector_alloc(mesh->get_max_deg() + 1);
    int s(0);
    gsl_permutation* p_gsl_ptr = gsl_permutation_alloc(mesh->get_max_deg() + 1);

    for (int l = 0; l != mesh->num_lines(); ++l) {
        if (mesh->get_bc_type(l) == DIRICHLET) {

            for (int j = 0; j != (mesh->get_max_deg() + 1); ++j) {
                for (int jj = 0; jj != (mesh->get_max_deg() + 1); ++jj) {
                    A[j * (mesh->get_max_deg() + 1) + jj] = 0;
                    for (int k = 0; k != num_integration_nodes; ++k)
                        A[j * (mesh->get_max_deg() + 1) + jj] += phi(k, j) * phi(k, jj) * weight(k);
                }
            }

            x0 = mesh->x_node(mesh->idx_node_line(l, 0) - mesh->idx_first_node());
            y0 = mesh->y_node(mesh->idx_node_line(l, 0) - mesh->idx_first_node());
            x1 = mesh->x_node(mesh->idx_node_line(l, 1) - mesh->idx_first_node());
            y1 = mesh->y_node(mesh->idx_node_line(l, 1) - mesh->idx_first_node());
            for (int j = 0; j != (mesh->get_max_deg() + 1); ++j) {
                b[j] = 0;
                for (int k = 0; k != num_integration_nodes; ++k)
                    b[j] += phi(k, j) * weight(k) * mesh->eval_bc(l, 0.5 * x0 * (1 - node(k)) + 0.5 * x1 * (1 + node(k)), 0.5 * y0 * (1 - node(k)) + 0.5 * y1 * (1 + node(k)), 0);
                b[j] *= sqrt(pow(0.5 * (x1 - x0), 2) + pow(0.5 * (y1 - y0), 2));
            }

            gsl_linalg_LU_decomp(&A_gsl.matrix, p_gsl_ptr, &s);
            gsl_linalg_LU_solve(&A_gsl.matrix, p_gsl_ptr, &b_gsl.vector, x_gsl_ptr);

            dirichlet_coeff[l].resize(static_cast<unsigned>(mesh->get_max_deg() + 1), 0);
            for (int j = 0; j != (mesh->get_max_deg() + 1); ++j)
                dirichlet_coeff[l][j] = gsl_vector_get(x_gsl_ptr, j);
        }
    }

    gsl_permutation_free(p_gsl_ptr);
    gsl_vector_free(x_gsl_ptr);
}

template <class M, class B, class G>
void equation<M, B, G>::disp_basis(const int& rank)
{
    if (mesh->get_rank() == rank) {
        for (vector<vector<vector<double>>>::iterator iter = phi_values.begin(); iter != phi_values.end(); ++iter) {
            cout << endl
                 << "phi_" << (iter - phi_values.begin()) << endl;
            mat_disp<double>(*iter);
        }
    }
}

template <class M, class B, class G>
void equation<M, B, G>::disp_der_basis(const int& rank)
{
    if (mesh->get_rank() == rank) {
        for (vector<vector<vector<double>>>::iterator iter = der_phi_values.begin(); iter != der_phi_values.end(); ++iter) {
            cout << endl
                 << "phi'_" << (iter - der_phi_values.begin()) << endl;
            mat_disp<double>(*iter);
        }
    }
}

template <class M, class B, class G>
void equation<M, B, G>::disp_dir_coeff(const int& rank)
{
    if (mesh->get_rank() == rank) {
        cout << "Coefficients of the projections of Dirichlet boundary functions on the modal basis" << endl
             << endl;
        for (vector<vector<double>>::iterator iter = dirichlet_coeff.begin(); iter != dirichlet_coeff.end(); ++iter) {
            cout << "line " << mesh->idx_line(iter - dirichlet_coeff.begin()) << ": ";
            for (vector<double>::iterator j = (*iter).begin(); j != (*iter).end(); ++j)
                cout << *j << " ";
            cout << endl;
        }
    }
}

#endif /* EQUATION_H_ */
