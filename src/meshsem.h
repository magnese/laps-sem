/*
 * meshsem.h
 *
 *  Created on: Feb 20, 2013
 *      Author: magnese
 */

#ifndef MESHSEM_H_
#define MESHSEM_H_

#include "meshparallel.h"
#include "utilities.h"

class mesh_sem : public mesh_parallel {

public:
    mesh_sem(const int& rank_value, const int& size_value, const int& n_refinements = 0, string file_name = "mesh.msh", flag_read flag1 = READ_MONO,
        flag_write flag2 = WRITE_OFF, const int& deg = 2)
        : mesh_parallel(rank_value, size_value, n_refinements, file_name, flag1, flag2)
        , basis_deg(num_elements())
        , max_deg(0)
    {
        if (deg < 2) {
            if (my_rank == 0)
                cout << "WARNING: Basis degree < 2, will be used degree = 2" << endl
                     << endl;
            set_basis_deg(2);
        } else
            set_basis_deg(deg);
    }

    inline int num_lines(void) { return my_lines.num_lines(); }

    inline int& i(const int& j) { return my_elements.i(j); }
    inline int& v(const int& pos, const int& j) { return my_elements.v(j)[pos]; }
    inline int& n(const int& pos, const int& j) { return my_elements.n(pos, j); }
    inline int& v_idx(const int& pos, const int& j) { return my_nodes.i(v(pos, j)); }
    inline int& idx_line(const int& j) { return my_lines.i(j); }
    inline int& idx_cor_line(const int& j) { return periodicity_lines_mapping[j] > 0 ? idx_line(periodicity_lines_mapping[j] - 1) : idx_line(-periodicity_lines_mapping[j] - 1); }
    inline int& line_ele_pos(const int& j) { return my_lines.elem_pos(j); }

    inline int& idx_first_node(void) { return my_nodes.i(0); }
    inline int& idx_first_line(void) { return my_lines.i(0); }
    inline int& idx_node_line(const int& j, const int& pos = 0) { return my_nodes.i(my_lines.p(j)[pos]); }
    inline int& idx_node_cor_line(const int& j, const int& pos = 0) { return idx_node_line(idx_cor_line(j) - idx_first_line(), (pos + 1) % 2); }

    inline double& x_node(const int& j) { return my_nodes.x(j); }
    inline double& y_node(const int& j) { return my_nodes.y(j); }

    int orientation(const int& value, const int& j);
    inline bool& is_ghost(const int& pos) { return my_elements.get_gf(pos); }

    inline void set_basis_deg(const int& value)
    {
        if (value > max_deg)
            max_deg = value;
        for (int i = 0; i != num_elements(); ++i)
            basis_deg[i] = value;
    }
    inline void set_basis_deg(const int& value, const int& pos)
    {
        basis_deg[pos] = value;
        if (value > max_deg)
            max_deg = value;
    }
    void update_ghosts_deg(void);

    inline vector<int>& get_deg(void) { return basis_deg; }
    inline int& get_deg(const int& j) { return basis_deg[j]; }
    inline int get_edge_deg(const int& pos, const int& j) { return (n(pos, j) < 0 ? basis_deg[j] : min(basis_deg[j], basis_deg[n(pos, j)])); }
    inline int& get_max_deg(void) { return max_deg; }

    inline void disp_line_periodicity(void) { vec_disp<int>(periodicity_lines_mapping); }

private:
    vector<int> basis_deg;
    int max_deg;
};

#endif /* MESHSEM_H_ */
