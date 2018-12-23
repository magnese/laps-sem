/*
 * nonlinearequation.h
 *
 *  Created on: May 24, 2013
 *      Author: Marco Agnese
 */

#ifndef NONLINEAREQUATION_H_
#define NONLINEAREQUATION_H_

#include "equation.h"

template <class M, class B, class G>
class non_linear_equation : public equation<M, B, G> {

public:
    non_linear_equation(M* mesh_ptr, double toll, double (*fptr)(double, double))
        : equation<M, B, G>(mesh_ptr)
        , initial_guess(fptr)
        , sol_der_x_value(0)
        , sol_der_y_value(0)
        , sol_value(0)
    {
        this->set_equation_type(NON_LINEAR_EQUATION);
        this->tollerance = toll;
        if (fptr == null_f)
            this->set_initial_guess_projected();
    }

    double& LHS(const int& i, const int& ii, const int& j, const int& jj, const int& ele_pos);
    double& RHS(const int& i, const int& ii, const int& ele_pos);
    double& neumann(const int& i, const int& l, const int& ele_pos);

    virtual ~non_linear_equation(void) {}

protected:
    double (*initial_guess)(double, double);

    double sol_der_x_value;
    double sol_der_y_value;
    double sol_value;

    void eval_sol_values(const int& k, const int& kk, const int& ele_pos);

    inline double& sol_x(void) { return sol_der_x_value; }
    inline double& sol_y(void) { return sol_der_y_value; }
    inline double& sol(void) { return sol_value; }

    virtual double eval_Jacobian(const int& k, const int& kk, const int& i, const int& ii, const int& j, const int& jj, const int& ele_pos) = 0;
    virtual double eval_F0(const int& k, const int& kk, const int& i, const int& ii, const int& ele_pos) = 0;
    virtual double eval_neumann(const double fun_value, const int& k, const int& i, const int& ele_pos) = 0;

private:
};

template <class M, class B, class G>
double& non_linear_equation<M, B, G>::LHS(const int& i, const int& ii, const int& j, const int& jj, const int& ele_pos)
{

    // TODO: for edges mode and vertex modes the necessary degree could be smaller and check aliasing problem

    this->set_order(this->mesh->get_deg(ele_pos));
    int num_integration_nodes(this->integration.necessary_size(this->get_order()));
    this->value = 0;

    for (int k = 0; k != num_integration_nodes; ++k) {
        for (int kk = 0; kk != num_integration_nodes; ++kk) {
            this->J_values[0] = this->mesh->get_mapping()->Jxx_1(this->node(k), this->node(kk), ele_pos);
            this->J_values[1] = this->mesh->get_mapping()->Jxy_1(this->node(k), this->node(kk), ele_pos);
            this->J_values[2] = this->mesh->get_mapping()->Jyx_1(this->node(k), this->node(kk), ele_pos);
            this->J_values[3] = this->mesh->get_mapping()->Jyy_1(this->node(k), this->node(kk), ele_pos);
            if (this->projection_status) {
                eval_sol_values(k, kk, ele_pos);
                this->value += eval_Jacobian(k, kk, i, ii, j, jj, ele_pos) * this->weight(k) * this->weight(kk) * fabs(pow(this->J(0) * this->J(3) - this->J(1) * this->J(2), -1));
            } else
                this->value += this->phi(k, j) * this->phi(kk, jj) * this->phi(k, i) * this->phi(kk, ii) * this->weight(k) * this->weight(kk) * fabs(pow(this->J(0) * this->J(3) - this->J(1) * this->J(2), -1));
        }
    }

    return this->value;
}

template <class M, class B, class G>
double& non_linear_equation<M, B, G>::RHS(const int& i, const int& ii, const int& ele_pos)
{

    // TODO: for edges mode and vertex modes the necessary degree could be smaller and check aliasing problem

    this->set_order(this->mesh->get_deg(ele_pos));
    int num_integration_nodes(this->integration.necessary_size(this->get_order()));
    this->value = 0;

    for (int k = 0; k != num_integration_nodes; ++k) {
        for (int kk = 0; kk != num_integration_nodes; ++kk) {
            this->J_values[0] = this->mesh->get_mapping()->Jxx_1(this->node(k), this->node(kk), ele_pos);
            this->J_values[1] = this->mesh->get_mapping()->Jxy_1(this->node(k), this->node(kk), ele_pos);
            this->J_values[2] = this->mesh->get_mapping()->Jyx_1(this->node(k), this->node(kk), ele_pos);
            this->J_values[3] = this->mesh->get_mapping()->Jyy_1(this->node(k), this->node(kk), ele_pos);
            if (this->projection_status) {
                eval_sol_values(k, kk, ele_pos);
                this->value -= eval_F0(k, kk, i, ii, ele_pos) * this->weight(k) * this->weight(kk) * fabs(pow(this->J(0) * this->J(3) - this->J(1) * this->J(2), -1));
            } else
                this->value += initial_guess(this->x_physical(k, kk, ele_pos), this->y_physical(k, kk, ele_pos)) * this->phi(k, i) * this->phi(kk, ii) * this->weight(k) * this->weight(kk) * fabs(pow(this->J(0) * this->J(3) - this->J(1) * this->J(2), -1));
        }
    }

    return this->value;
}

template <class M, class B, class G>
double& non_linear_equation<M, B, G>::neumann(const int& i, const int& l, const int& ele_pos)
{
    //TODO

    this->value = 0;

    if (!(this->projection_status)) {

        this->set_order(this->mesh->get_deg(ele_pos));
        int num_integration_nodes(this->integration.necessary_size(this->get_order()));
        double x0(this->mesh->x_node(this->mesh->idx_node_line(l, 0) - this->mesh->idx_first_node()));
        double x1(this->mesh->y_node(this->mesh->idx_node_line(l, 0) - this->mesh->idx_first_node()));
        double y0(this->mesh->x_node(this->mesh->idx_node_line(l, 1) - this->mesh->idx_first_node()));
        double y1(this->mesh->y_node(this->mesh->idx_node_line(l, 1) - this->mesh->idx_first_node()));

        for (int k = 0; k != num_integration_nodes; ++k)
            this->value += this->weight(k) * eval_neumann(this->mesh->eval_bc(-this->mesh->n(l, ele_pos) - 1, 0.5 * x0 * (1 - this->node(k)) + 0.5 * x1 * (1 + this->node(k)), 0.5 * y0 * (1 - this->node(k)) + 0.5 * y1 * (1 + this->node(k)), 0), k, i, ele_pos);

        this->value *= sqrt(pow(0.5 * (x1 - x0), 2) + pow(0.5 * (y1 - y0), 2));
    }

    return this->value;
}

template <class M, class B, class G>
void non_linear_equation<M, B, G>::eval_sol_values(const int& k, const int& kk, const int& ele_pos)
{

    sol_der_x_value = 0;
    sol_der_y_value = 0;
    sol_value = 0;

    for (int j = 0; j != (this->mesh->get_deg(ele_pos) + 1); ++j) {
        for (int jj = 0; jj != (this->mesh->get_deg(ele_pos) + 1); ++jj) {
            sol_der_x_value += this->map->u(j, jj, ele_pos) * this->der_phi(k, j) * this->phi(kk, jj);
            sol_der_y_value += this->map->u(j, jj, ele_pos) * this->phi(k, j) * this->der_phi(kk, jj);
            sol_value += this->map->u(j, jj, ele_pos) * this->phi(k, j) * this->phi(kk, jj);
        }
    }
}

#endif /* NONLINEAREQUATION_H_ */
