/*
 * ghostscommunication.h
 *
 *  Created on: Feb 17, 2013
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef GHOSTSCOMMUNICATION_H_
#define GHOSTSCOMMUNICATION_H_

#include "include.h"

using namespace std;

template <class M>
class ghosts_communication {

public:
    ghosts_communication(M* mesh_ptr)
        : my_mesh(mesh_ptr)
        , data_size(0)
        , global_vec()
        , local_vec()
        , ghosts_idx(0)
        , local_ghosts_positions(0)
        , ghosts_values(0)
    {
        VecCreate(MPI_COMM_WORLD, &local_vec);
    }

    inline void set_size(const int& new_size)
    {
        data_size = new_size;
        set_idxs();
        VecCreateGhost(MPI_COMM_WORLD, data_size * my_mesh->get_elements().num_elements(), PETSC_DECIDE, data_size * my_mesh->get_elements().num_ghosts(), &(ghosts_idx[0]), &global_vec);
    }
    void set_idxs(void);
    void set_global_vec(vector<double>& values, const int& vec_start, const int& vec_end);
    inline void set_global_vec(vector<double>& values) { set_global_vec(values, 0, my_mesh->num_elements()); }

    inline vector<double>& get_ghosts_values(void) { return ghosts_values; }
    inline double& get_ghost_value(const int& j) { return ghosts_values[j]; }

    void disp_ghosts_values(void);
    vector<double>& get_ghosts_global_idx(void);

    inline void free(void)
    {
        PetscFree(global_vec);
        VecDestroy(&local_vec);
    }

private:
    M* my_mesh;
    int data_size;
    Vec global_vec;
    Vec local_vec;
    vector<int> ghosts_idx;
    vector<int> ix;
    vector<int> local_ghosts_positions;
    vector<double> ghosts_values;
};

template <class M>
void ghosts_communication<M>::set_idxs(void)
{
    ghosts_idx.resize(static_cast<unsigned>(data_size * my_mesh->num_ghosts()));
    ghosts_values.resize(static_cast<unsigned>(data_size * my_mesh->num_ghosts()));
    for (int j = my_mesh->num_elements(); j != my_mesh->elements_size(); ++j) {
        for (int k = 0; k != data_size; ++k) {
            ghosts_idx[(j - my_mesh->num_elements()) * data_size + k] = data_size * (my_mesh->get_elements().ele(j).i() - my_mesh->idx_first_global_ele()) + k;
        }
    }
    ix.resize(static_cast<unsigned>(data_size * my_mesh->num_elements()));
    for (int j = 0; j != my_mesh->num_elements(); ++j) {
        for (int k = 0; k != data_size; ++k) {
            ix[j * data_size + k] = (my_mesh->get_elements().ele(j).i() - my_mesh->idx_first_global_ele()) * data_size + k;
        }
    }
    local_ghosts_positions.resize(static_cast<unsigned>(data_size * my_mesh->num_ghosts()));
    for (vector<int>::iterator iter = local_ghosts_positions.begin(); iter != local_ghosts_positions.end(); ++iter) {
        *iter = (iter - local_ghosts_positions.begin()) + data_size * my_mesh->num_elements();
    }
}

template <class M>
void ghosts_communication<M>::set_global_vec(vector<double>& values, const int& vec_start, const int& vec_end)
{
    VecSetValues(global_vec, data_size * (vec_end - vec_start), &(ix[0]), &(values[vec_start]), INSERT_VALUES);
    VecAssemblyBegin(global_vec);
    VecAssemblyEnd(global_vec);
    VecGhostUpdateBegin(global_vec, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(global_vec, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostGetLocalForm(global_vec, &local_vec);
    VecGetValues(local_vec, data_size * my_mesh->num_ghosts(), &(local_ghosts_positions[0]), &(ghosts_values[0]));
    VecGhostRestoreLocalForm(global_vec, &local_vec);
}

template <class M>
void ghosts_communication<M>::disp_ghosts_values(void)
{
    cout << "ghosts values" << endl;
    for (vector<double>::iterator iter = ghosts_values.begin(); iter != ghosts_values.end(); ++iter)
        cout << *iter << " ";
    cout << endl
         << endl;
}

template <class M>
vector<double>& ghosts_communication<M>::get_ghosts_global_idx(void)
{
    vector<double> values;
    set_size(1);
    values.resize(static_cast<unsigned>(my_mesh->num_elements()));
    for (vector<double>::iterator iter = values.begin(); iter != values.end(); ++iter)
        *iter = static_cast<double>(my_mesh->get_elements().ele(iter - values.begin()).i());
    set_global_vec(values);
    return ghosts_values;
}

#endif /* GHOSTSCOMMUNICATION_H_ */
