/*
 * elementslist.h
 *
 *  Created on: Dec 3, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef ELEMENTSLIST_H_
#define ELEMENTSLIST_H_

#include "element.h"
#include "include.h"

using namespace std;

class elements_list {

public:
    elements_list(void)
        : my_elements(0)
        , number_ghosts(0)
    {
    }

    inline int num_elements(void) { return static_cast<int>(my_elements.size()) - number_ghosts; }
    inline int num_ghosts(void) { return number_ghosts; }
    inline vector<int>& v(const int& j) { return my_elements[j].v(); }
    inline vector<int>& n(const int& j) { return my_elements[j].n(); }
    inline int& n(const int& pos, const int& j) { return my_elements[j].n(pos); }
    inline int& i(const int& j) { return my_elements[j].i(); }
    inline element& ele(const int& i) { return my_elements[i]; }

    inline void resize(const int& value) { my_elements.resize(static_cast<unsigned>(value)); }
    inline void set_num_ghosts(const int& value) { number_ghosts = value; }
    inline void set_v(vector<int>& values, const int& j) { my_elements[j].set_v(values); }
    inline void set_v(const int& value, const int& i, const int& j) { (my_elements[j].set_v(value, i)); }
    inline void set_n(vector<int>& values, const int& j) { my_elements[j].set_n(values); }
    inline void set_n(const int& value, const int& pos, const int& j) { my_elements[j].set_n(value, pos); }
    inline void set_i(const int& value, const int& j) { my_elements[j].set_i(value); }
    inline void add_n(const int& value, const int& j) { my_elements[j].add_n(value); }
    inline void set_gf(const int& j) { my_elements[j].set_gf(); }
    inline bool& get_gf(const int& j) { return my_elements[j].get_gf(); }

    inline int& tag(const int& j) { return my_elements[j].tag(); }
    inline void set_tag(const int& value, const int& j) { my_elements[j].set_tag(value); }

    inline void add_element(element& new_element) { my_elements.push_back(new_element); }
    inline void add_ghost(element& new_element)
    {
        my_elements.push_back(new_element);
        ++number_ghosts;
    }
    inline void free_ghosts(void) { my_elements.resize(static_cast<unsigned>(num_elements())); }
    inline void free_neighbors(void)
    {
        for (vector<element>::iterator j = my_elements.begin(); j != my_elements.end(); ++j) {
            j->resize_neighbors(0, 0);
            j->free_pos();
        }
    }
    inline void free(void)
    {
        my_elements.resize(0);
        number_ghosts = 0;
    }

    inline void set_pos(vector<int>& values, const int& j) { my_elements[j].set_pos(values); }
    inline void set_pos(int& values, const int& j, const int& i) { my_elements[j].set_pos(values, i); }
    inline void free_pos(const int& j) { my_elements[j].free_pos(); }
    inline int c_n(const int& position, const int& corner_index, const int& j) { return my_elements[j].c_n(position, corner_index); }
    inline int num_neighbors(const int& j) { return my_elements[j].num_neighbors(); }
    inline void resize_neighbors(const int& l, const int& values, const int& j) { my_elements[j].resize_neighbors(l, values); }
    inline int num_corners(const int& pos, const int& j) { return my_elements[j].num_corners(pos); }
    inline int num_corners_total(const int& j) { return my_elements[j].num_corners_total(); }
    inline void inc_edge_position(const int& l, const int& j) { my_elements[j].inc_edge_position(l); }

    inline int& ncp(const int& pos, const int& j) { return my_elements[j].ncp(pos); }
    inline void set_ncp(const int& value, const int& pos, const int& j) { my_elements[j].set_ncp(value, pos); }

    inline vector<int>& c_pos(const int& j) { return my_elements[j].c_pos(); }

    void add_from_file(const char* file_name);
    void copy(const int& l, const int& j);
    void print(void);

private:
    vector<element> my_elements;
    int number_ghosts;
};

#endif /* ELEMENSTLIST_H_ */
