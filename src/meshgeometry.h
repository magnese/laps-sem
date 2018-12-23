/*
 * meshgeometry.h
 *
 *  Created on: Dec 12, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef MESHGEOMETRY_H_
#define MESHGEOMETRY_H_

#include "include.h"
#include "mesh.h"

using namespace std;

class mesh_geometry : public mesh {

public:
    mesh_geometry()
        : mesh()
        , structure_param(0) {};
    void read_nodes_list(string file_name);
    void write_msh(string file_name);
    void create_lines(void);
    void write_geo(string file_name);
    void split_mesh(string file_name);
    void read_mesh_struct(string file_name, const int& lines_of_param);
    void create_elements(const int& tag_value, const loop& flag = OPENED);

protected:
private:
    vector<int> structure_param;
};

#endif /* MESHGEOMETRY_H_ */
