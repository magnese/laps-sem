/*
 * modesmapping.h
 *
 *  Created on: Jun 30, 2013
 *      Author: magnese
 */

#ifndef MODESMAPPING_H_
#define MODESMAPPING_H_

#include "utilities.h"

using namespace std;

template <class M, class B, class P>
class modes_mapping {

public:
    modes_mapping(M* mesh_ptr, B* basis_ptr, P* p_state_ptr)
        : mesh(mesh_ptr)
        , basis(basis_ptr)
        , p_state(p_state_ptr)
        , node_mode_mapping(mesh->num_nodes(), -1)
        , num_dof_nodes(0)
        , num_not_dof_nodes(0)
        , mode_mapping_starting_idx(0)
        , elmdist_dof_edges(mesh->get_size() + 1, 0)
        , num_dof_edges(0)
        , elmdist_not_dof_edges(mesh->get_size() + 1, 0)
        , num_not_dof_edges(0)
        , elmdist_dof_interiors(mesh->get_size() + 1, 0)
        , num_dof_interiors(0)
        , vtx_map(4)
        , offset_inner_modes(0)
        , dim_boundary_modes(0)
        , sol_coeff(0)
        , coeff(0)
        , val(0)
    {
        vtx_map[0] = 3;
        vtx_map[1] = 0;
        vtx_map[2] = 2;
        vtx_map[3] = 1;
    }

    int edge_pos(const int& i, const int& ii);
    inline int& vtx_pos(const int& i, const int& ii) { return vtx_map[2 * i + ii]; }

    void create_mapping(bool static_condensation_status);

    inline bool is_v(const int& i, const int& ii)
    {
        if ((i < 2) && (ii < 2))
            return true;
        else
            return false;
    }
    inline bool is_e(const int& i, const int& ii)
    {
        if (((i < 2) && (ii > 1)) || ((i > 1) && (ii < 2)))
            return true;
        else
            return false;
    }
    inline bool is_i(const int& i, const int& ii)
    {
        if ((i > 1) && (ii > 1))
            return true;
        else
            return false;
    }

    inline int vtx(const int& i, const int& ii, const int& e) { return node_mode_mapping[mesh->v(vtx_pos(i, ii), e)]; }
    inline int& vtx(const int& pos) { return node_mode_mapping[pos]; }
    inline int edge(const int& i, const int& ii, const int& e) { return mode_mapping_starting_idx[5 * e + edge_pos(i, ii)] + (i > ii ? i : ii) - 2; }
    inline int inner(const int& i, const int& ii, const int& e)
    {
        return mode_mapping_starting_idx[5 * e + 4] + (mesh->get_deg(e) - 1) * (i - 2) + ii - 2 - offset_inner_modes;
    }

    inline int& dim_boundary(void) { return dim_boundary_modes; }
    inline int& dim_interiors(void) { return num_dof_interiors; }
    inline int& dim_dof_nodes(void) { return num_dof_nodes; }
    inline int& dim_dof_edges(void) { return num_dof_edges; }
    inline int& dim_ndof_nodes(void) { return num_not_dof_nodes; }
    inline int& dim_ndof_edges(void) { return num_not_dof_edges; }

    inline vector<int>& get_nodes_mapping(void) { return node_mode_mapping; }
    inline void set_ptr_sol_coeff(double** ptr_values) { sol_coeff = ptr_values; }

    double& u(const int& i, const int& ii, const int& e);
    inline double& u(const int& pos) { return (*sol_coeff)[pos]; }
    double& u_value(double& x, double& y, const int& e);

    void disp_mapping(void);

protected:
private:
    M* mesh;
    B* basis;
    P* p_state;

    vector<int> node_mode_mapping;
    int num_dof_nodes;
    int num_not_dof_nodes;

    vector<int> mode_mapping_starting_idx;
    vector<int> elmdist_dof_edges;
    int num_dof_edges;
    vector<int> elmdist_not_dof_edges;
    int num_not_dof_edges;
    vector<int> elmdist_dof_interiors;
    int num_dof_interiors;

    vector<int> vtx_map;
    int offset_inner_modes;
    int dim_boundary_modes;

    double** sol_coeff;
    double coeff;
    double val;
};

template <class M, class B, class P>
void modes_mapping<M, B, P>::create_mapping(bool static_condensation_status)
{

    for (int l = 0; l != mesh->num_lines(); ++l) {
        if (mesh->get_bc_type(l) == DIRICHLET) {
            if (node_mode_mapping[mesh->idx_node_line(l, 0) - mesh->idx_first_node()] < 0) {
                node_mode_mapping[mesh->idx_node_line(l, 0) - mesh->idx_first_node()] = num_not_dof_nodes;
                ++num_not_dof_nodes;
            }
            if (node_mode_mapping[mesh->idx_node_line(l, 1) - mesh->idx_first_node()] < 0) {
                node_mode_mapping[mesh->idx_node_line(l, 1) - mesh->idx_first_node()] = num_not_dof_nodes;
                ++num_not_dof_nodes;
            }
            if (mesh->get_bc_type((l + mesh->num_lines() - 1) % mesh->num_lines()) == PERIODICITY) {
                if (node_mode_mapping[mesh->idx_node_cor_line((l + mesh->num_lines() - 1) % mesh->num_lines(), 1) - mesh->idx_first_node()] < 0)
                    node_mode_mapping[mesh->idx_node_cor_line((l + mesh->num_lines() - 1) % mesh->num_lines(), 1) - mesh->idx_first_node()] = node_mode_mapping[mesh->idx_node_line((l + mesh->num_lines() - 1) % mesh->num_lines(), 1) - mesh->idx_first_node()];
            }
            if (mesh->get_bc_type((l + 1) % mesh->num_lines()) == PERIODICITY) {
                if (node_mode_mapping[mesh->idx_node_cor_line((l + 1) % mesh->num_lines(), 0) - mesh->idx_first_node()] < 0)
                    node_mode_mapping[mesh->idx_node_cor_line((l + 1) % mesh->num_lines(), 0) - mesh->idx_first_node()] = node_mode_mapping[mesh->idx_node_line((l + 1) % mesh->num_lines(), 0) - mesh->idx_first_node()];
            }
        }
    }

    for (int j = 0; j != mesh->num_elements(); ++j) {
        num_dof_interiors += (mesh->get_deg(j) - 1) * (mesh->get_deg(j) - 1);
        for (int jj = 0; jj != 4; ++jj) {
            if (mesh->n(jj, j) > -1) {
                if (mesh->i(j) < mesh->i(mesh->n(jj, j)))
                    num_dof_edges += mesh->get_edge_deg(jj, j) - 1;
            } else {
                if (mesh->get_bc_type(-mesh->n(jj, j) - 1) == DIRICHLET)
                    num_not_dof_edges += mesh->get_deg(j) - 1;
                else
                    num_dof_edges += mesh->get_deg(j) - 1;
            }
        }
    }

    MPI_Gather(&num_not_dof_edges, 1, MPI_INT, &(elmdist_not_dof_edges[1]), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mesh->get_rank() == 0) {
        elmdist_not_dof_edges[0] = num_not_dof_nodes;
        for (vector<int>::iterator iter = (elmdist_not_dof_edges.begin() + 1); iter != elmdist_not_dof_edges.end(); ++iter)
            *iter += *(iter - 1);
    }
    MPI_Bcast(&(elmdist_not_dof_edges[0]), mesh->get_size() + 1, MPI_INT, 0, MPI_COMM_WORLD);
    num_not_dof_edges = *(elmdist_not_dof_edges.end() - 1) - *(elmdist_not_dof_edges.begin());

    num_dof_nodes = num_not_dof_nodes + num_not_dof_edges;

    for (int j = 0; j != mesh->num_lines(); ++j) {
        if (mesh->get_bc_type(j) == PERIODICITY) {
            for (int jj = 0; jj != 2; ++jj) {
                if (((node_mode_mapping[mesh->idx_node_line(j, jj) - mesh->idx_first_node()]) < 0) && ((node_mode_mapping[mesh->idx_node_cor_line(j, jj) - mesh->idx_first_node()]) < 0)) {
                    node_mode_mapping[mesh->idx_node_line(j, jj) - mesh->idx_first_node()] = num_dof_nodes;
                    node_mode_mapping[mesh->idx_node_cor_line(j, jj) - mesh->idx_first_node()] = num_dof_nodes;
                    ++num_dof_nodes;
                } else {
                    if ((node_mode_mapping[mesh->idx_node_line(j, jj) - mesh->idx_first_node()]) > -1)
                        node_mode_mapping[mesh->idx_node_cor_line(j, jj) - mesh->idx_first_node()] = node_mode_mapping[mesh->idx_node_line(j, jj) - mesh->idx_first_node()];
                    else {
                        if ((node_mode_mapping[mesh->idx_node_cor_line(j, jj) - mesh->idx_first_node()]) > -1)
                            node_mode_mapping[mesh->idx_node_line(j, jj) - mesh->idx_first_node()] = node_mode_mapping[mesh->idx_node_cor_line(j, jj) - mesh->idx_first_node()];
                    }
                }
            }
        }
    }

    for (vector<int>::iterator iter = node_mode_mapping.begin(); iter != node_mode_mapping.end(); ++iter) {
        if (*iter < 0) {
            *iter = num_dof_nodes;
            ++num_dof_nodes;
        }
    }
    num_dof_nodes = num_dof_nodes - num_not_dof_nodes - num_not_dof_edges;

    MPI_Gather(&num_dof_edges, 1, MPI_INT, &(elmdist_dof_edges[1]), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mesh->get_rank() == 0) {
        elmdist_dof_edges[0] = num_dof_nodes + num_not_dof_nodes + num_not_dof_edges;
        for (vector<int>::iterator iter = (elmdist_dof_edges.begin() + 1); iter != elmdist_dof_edges.end(); ++iter)
            *iter += *(iter - 1);
    }
    MPI_Bcast(&(elmdist_dof_edges[0]), mesh->get_size() + 1, MPI_INT, 0, MPI_COMM_WORLD);
    num_dof_edges = *(elmdist_dof_edges.end() - 1) - *(elmdist_dof_edges.begin());

    dim_boundary_modes = num_dof_nodes + num_not_dof_nodes + num_dof_edges + num_not_dof_edges;
    if (static_condensation_status)
        offset_inner_modes = dim_boundary_modes;

    MPI_Gather(&num_dof_interiors, 1, MPI_INT, &(elmdist_dof_interiors[1]), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mesh->get_rank() == 0) {
        elmdist_dof_interiors[0] = num_dof_nodes + num_not_dof_nodes + num_dof_edges + num_not_dof_edges;
        for (vector<int>::iterator iter = (elmdist_dof_interiors.begin() + 1); iter != elmdist_dof_interiors.end(); ++iter)
            *iter += *(iter - 1);
    }
    MPI_Bcast(&(elmdist_dof_interiors[0]), mesh->get_size() + 1, MPI_INT, 0, MPI_COMM_WORLD);
    num_dof_interiors = *(elmdist_dof_interiors.end() - 1) - *(elmdist_dof_interiors.begin());

    mode_mapping_starting_idx.resize(static_cast<unsigned>(5 * mesh->num_elements()), -1);
    int cont(elmdist_dof_edges[mesh->get_rank()]);
    int cont_not(elmdist_not_dof_edges[mesh->get_rank()]);
    int cont_int(elmdist_dof_interiors[mesh->get_rank()]);
    for (int j = 0; j != mesh->num_elements(); ++j) {
        mode_mapping_starting_idx[j * 5 + 4] = cont_int;
        cont_int += (mesh->get_deg(j) - 1) * (mesh->get_deg(j) - 1);
        for (int jj = 0; jj != 4; ++jj) {
            if (mesh->n(jj, j) > -1) {
                if (mesh->i(j) < mesh->i(mesh->n(jj, j))) {
                    mode_mapping_starting_idx[j * 5 + jj] = cont;
                    cont += mesh->get_edge_deg(jj, j) - 1;
                } else {
                    if (!mesh->is_ghost(mesh->n(jj, j))) {
                        mode_mapping_starting_idx[j * 5 + jj] = mode_mapping_starting_idx[mesh->n(jj, j) * 5 + mesh->neighbor_orientation(jj, j)];
                    }
                }
            } else {
                if (mesh->get_bc_type(-mesh->n(jj, j) - 1) == DIRICHLET) {
                    mode_mapping_starting_idx[j * 5 + jj] = cont_not;
                    cont_not += mesh->get_deg(j) - 1;
                } else {
                    mode_mapping_starting_idx[j * 5 + jj] = cont;
                    cont += mesh->get_deg(j) - 1;
                }
            }
        }
    }

    mesh->get_ghosts_comm()->set_size(5);
    vector<double> vector_double(mode_mapping_starting_idx.begin(), mode_mapping_starting_idx.end());
    mesh->get_ghosts_comm()->set_global_vec(vector_double);
    for (int j = 0; j != mesh->num_elements(); ++j) {
        for (int jj = 0; jj != 5; ++jj) {
            if (mode_mapping_starting_idx[j * 5 + jj] < 0) {
                mode_mapping_starting_idx[j * 5 + jj] = static_cast<int>(mesh->get_ghosts_comm()->get_ghost_value((mesh->n(jj, j) - mesh->num_elements()) * 5 + mesh->neighbor_orientation(jj, j)));
            }
        }
    }
}

template <class M, class B, class P>
int modes_mapping<M, B, P>::edge_pos(const int& i, const int& ii)
{

    int temp(0);

    if (ii == 1)
        temp = 0;
    if (i == 1)
        temp = 1;
    if (ii == 0)
        temp = 2;
    if (i == 0)
        temp = 3;

    return temp;
}

template <class M, class B, class P>
void modes_mapping<M, B, P>::disp_mapping(void)
{
    if (mesh->get_rank() == 0) {
        cout << "Node modes mapping:" << endl;
        node_vec_disp<int, M>(node_mode_mapping, 1, mesh);
    }
    for (int i = 0; i != mesh->get_size(); ++i) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (mesh->get_rank() == i) {
            cout << endl
                 << "Edge modes mapping on process " << i;
            mesh_vec_disp<int, M>(mode_mapping_starting_idx, 5, mesh);
        }
    }
}

template <class M, class B, class P>
double& modes_mapping<M, B, P>::u(const int& i, const int& ii, const int& e)
{

    if (is_v(i, ii))
        coeff = (*sol_coeff)[vtx(i, ii, e)];
    else {
        if (is_i(i, ii)) {
            if ((i > mesh->get_deg(e)) || (ii > mesh->get_deg(e)))
                coeff = 0;
            else
                coeff = (*sol_coeff)[inner(i, ii, e)];
        } else {
            if (max(i, ii) > mesh->get_edge_deg(edge_pos(i, ii), e))
                coeff = 0;
            else {
                coeff = (*sol_coeff)[edge(i, ii, e)];
                coeff *= pow(mesh->orientation(edge_pos(i, ii), e), max(i, ii));
            }
        }
    }

    return coeff;
}

template <class M, class B, class P>
double& modes_mapping<M, B, P>::u_value(double& x, double& y, const int& e)
{

    val = 0;

    for (int i = 0; i != (mesh->get_deg(e) + 1); ++i) {
        for (int ii = 0; ii != (mesh->get_deg(e) + 1); ++ii)
            val += u(i, ii, e) * basis->eval(x, i) * basis->eval(y, ii);
    }

    return val;
}

#endif /* MODESMAPPING_H_ */
