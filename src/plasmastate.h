/*
 * plasmastate.h
 *
 *  Created on: Feb 19, 2013
 *      Author: Marco Agnese
 */

#ifndef PLASMASTATE_H_
#define PLASMASTATE_H_

#include "modesmapping.h"
#include "utilities.h"

using namespace std;

template <class M, class B>
class plasma_state {

public:
    plasma_state(M* mesh_ptr, B* basis_ptr)
        : mesh(mesh_ptr)
        , basis(basis_ptr)
        , map(0)
        , nodal_values(mesh->num_nodes())
        , num_sample_points(0)
        , non_nodal_values(0)
    {
    }

    inline M* get_mesh(void) { return mesh; }
    inline B* get_basis(void) { return basis; }

    inline void set_ptr_mode_mapping(modes_mapping<M, B, plasma_state<M, B>>* ptr) { map = ptr; }

    void update(const int& num = 0);

    inline double& get_nodal_value(const int& pos) { return nodal_values[pos]; }
    inline vector<double>& get_nodal_values(void) { return nodal_values; }
    void disp_nodes_values(void);

    inline int num_sol_values(void) { return non_nodal_values.size(); }
    inline vector<double> get_sol(const int& j) { return non_nodal_values[j]; }

protected:
private:
    M* mesh;
    B* basis;
    modes_mapping<M, B, plasma_state<M, B>>* map;

    vector<double> nodal_values;

    int num_sample_points;
    vector<vector<double>> non_nodal_values;
};

template <class M, class B>
void plasma_state<M, B>::disp_nodes_values(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (mesh->get_rank() == 0) {
        cout << "Plasma state nodes values:" << endl
             << endl;
        for (int j = 0; j != mesh->num_nodes(); ++j)
            cout << mesh->get_nodes().i(j) << ": " << nodal_values[j] << endl;
        cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

template <class M, class B>
void plasma_state<M, B>::update(const int& num)
{

    for (int i = 0; i != mesh->num_nodes(); ++i)
        nodal_values[i] = map->u(map->vtx(i));

    num_sample_points = num;

    vector<double> null_vec(3, 0);
    non_nodal_values.resize(static_cast<unsigned>(num_sample_points * (num_sample_points + 4) * mesh->num_elements()), null_vec);

    int cont(0);
    vector<double> logical_coor(num_sample_points + 2, 0);
    for (int i = 0; i != (num_sample_points + 2); ++i)
        logical_coor[i] = -1 + (2 / (num_sample_points + 1)) * i;

    double flag(false);

    for (int j = 0; j != mesh->num_elements(); ++j) {

        for (int k = 0; k != (num_sample_points + 2); ++k) {
            for (int kk = 0; kk != (num_sample_points + 2); ++kk) {

                flag = false;

                if ((k > 0) && (k < (num_sample_points + 1)) && (kk > 0) && (kk < (num_sample_points + 1)))
                    flag = true;
                else {
                    if ((k > 0) && (k < (num_sample_points + 1)) && (kk == (num_sample_points + 1))) {
                        if (mesh->n(0, j) < 0)
                            flag = true;
                        else {
                            if (mesh->i(j) > mesh->i(mesh->n(0, j)))
                                flag = true;
                        }
                    } else {
                        if ((k == (num_sample_points + 1)) && (kk > 0) && (kk < (num_sample_points + 1))) {
                            if (mesh->n(1, j) < 0)
                                flag = true;
                            else {
                                if (mesh->i(j) > mesh->i(mesh->n(1, j)))
                                    flag = true;
                            }
                        } else {
                            if ((k > 0) && (k < (num_sample_points + 1)) && (kk == 0)) {
                                if (mesh->n(2, j) < 0)
                                    flag = true;
                                else {
                                    if (mesh->i(j) > mesh->i(mesh->n(2, j)))
                                        flag = true;
                                }
                            } else {
                                if ((k == 0) && (kk > 0) && (kk < (num_sample_points + 1))) {
                                    if (mesh->n(3, j) < 0)
                                        flag = true;
                                    else {
                                        if (mesh->i(j) > mesh->i(mesh->n(3, j)))
                                            flag = true;
                                    }
                                }
                            }
                        }
                    }
                }

                if (flag) {
                    non_nodal_values[cont][0] = mesh->get_mapping()->Fx(logical_coor[k], logical_coor[kk], j);
                    non_nodal_values[cont][1] = mesh->get_mapping()->Fy(logical_coor[k], logical_coor[kk], j);
                    non_nodal_values[cont][2] = map->u_value(logical_coor[k], logical_coor[kk], j);
                    ++cont;
                }
            }
        }
    }

    non_nodal_values.resize(cont);
}

#endif /* PLASMASTATE_H_ */
