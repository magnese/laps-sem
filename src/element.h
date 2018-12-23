/*
 * element.h
 *
 *  Created on: Oct 17, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include "include.h"

using namespace std;

class element {

public:
    element()
        : global_index(0)
        , ghost_flag(0)
        , vertex_indices(4, 0)
        , neighboring_elements(0)
        , corner_position(4, 4)
        , surface_tag(0)
        , neighbor_cor_pos(4, -1) {};

    inline int& i(void) { return global_index; }
    inline void set_i(const int& value) { global_index = value; }

    inline vector<int>& v(void) { return vertex_indices; }
    inline void set_v(vector<int>& values) { vertex_indices = values; }
    inline void set_v(const int& values, const int& pos) { vertex_indices[pos] = values; }

    inline vector<int>& n(void) { return neighboring_elements; }
    inline int& n(const int& pos) { return neighboring_elements[pos]; }
    inline void set_n(vector<int>& values) { neighboring_elements = values; }
    inline void set_n(const int& value, const int& pos) { neighboring_elements[pos] = value; }
    inline void add_n(const int& value) { neighboring_elements.push_back(value); }

    inline void set_gf(void) { ghost_flag = 1; }
    inline bool& get_gf(void) { return ghost_flag; }

    inline void set_pos(vector<int>& values) { corner_position = values; }
    inline void set_pos(int& values, const int& j) { corner_position[j] = values; }

    inline void free_pos(void)
    {
        for (int j = 0; j != 4; ++j)
            corner_position[j] = 4;
    }
    int c_n(const int& position, const int& corner_index);

    inline int num_neighbors(void) { return neighboring_elements.size(); }
    inline void resize_neighbors(const int& l, const int& values) { neighboring_elements.resize(l, values); }

    inline int& tag(void) { return surface_tag; }
    inline void set_tag(const int& value) { surface_tag = value; }

    inline int& ncp(const int& pos) { return neighbor_cor_pos[pos]; }
    inline void set_ncp(const int& value, const int& pos) { neighbor_cor_pos[pos] = value; }

    int num_corners(const int& pos);
    inline int num_corners_total(void) { return neighboring_elements.size() - 4; }

    inline void inc_edge_position(const int& l)
    {
        for (int j = l; j != 4; ++j)
            corner_position[j] += 1;
    }

    inline vector<int>& c_pos(void) { return corner_position; }

private:
    int global_index;
    bool ghost_flag;
    vector<int> vertex_indices; // SW vertex, SE, NE, NW
    vector<int> neighboring_elements; // S element, E, N, W
    vector<int> corner_position; // position of the position+1 of the last corner of each vertex followed by the corners in the order SW, SE, NE, NW
    int surface_tag;
    vector<int> neighbor_cor_pos;
};

#endif /* ELEMENT_H_ */
