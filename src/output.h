/*
 * output.h
 *
 *  Created on: Mar 15, 2013
 *      Author: Marco Agnese
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "utilities.h"

using namespace std;

template <class M, class PS>
class output {

public:
    output(M* mesh_ptr, PS* ps_ptr)
        : mesh(mesh_ptr)
        , p_state(ps_ptr) {};

    inline M* get_mesh(void) { return mesh; }
    inline PS* get_ps(void) { return p_state; }

    void gmsh_dump(string file_name);
    void dump(string file_name, bool all_values = true);
    void dump_binary(string file_name, bool all_values = true);

protected:
private:
    M* mesh;
    PS* p_state;
};

template <class M, class PS>
void output<M, PS>::gmsh_dump(string file_name)
{

    stringstream my_stringstream;
    string rank_n;

    file_name.resize(file_name.size() - 4);
    my_stringstream << mesh->get_rank();
    rank_n = my_stringstream.str();
    file_name = file_name + "_" + rank_n + ".msh";

    if (mesh->get_rank() == 0)
        mesh->write_file(file_name, WRITE_GHOSTS_OFF, p_state->get_nodal_values());
    else
        mesh->write_file(file_name, WRITE_GHOSTS_OFF);
}

template <class M, class PS>
void output<M, PS>::dump(string file_name, bool all_values)
{

    stringstream my_stringstream;
    string rank;
    my_stringstream << (mesh->get_rank());
    rank = my_stringstream.str();

    string file_name_sol(file_name + "_" + rank + ".dat");
    ofstream file;
    file.open(file_name_sol.c_str());

    if (file.is_open()) {

        if (mesh->get_rank() == 0) {
            //TODO: each process insert a portion!
            for (int j = 0; j != mesh->num_nodes(); ++j)
                file << mesh->get_nodes().x(j) << " " << mesh->get_nodes().y(j) << " " << p_state->get_nodal_value(j) << endl;
        }

        if (all_values) {
            for (int j = 0; j != p_state->num_sol_values(); ++j)
                file << p_state->get_sol(j)[0] << " " << p_state->get_sol(j)[1] << " " << p_state->get_sol(j)[2] << endl;
        }

        file.close();
    }
}

template <class M, class PS>
void output<M, PS>::dump_binary(string file_name, bool all_values)
{

    stringstream my_stringstream;
    string rank;
    my_stringstream << (mesh->get_rank());
    rank = my_stringstream.str();

    string file_name_x(file_name + "_" + rank + ".x");
    ofstream file_x;
    file_x.open(file_name_x.c_str());

    string file_name_y(file_name + "_" + rank + ".y");
    ofstream file_y;
    file_y.open(file_name_y.c_str());

    string file_name_sol(file_name + "_" + rank + ".z");
    ofstream file_sol;
    file_sol.open(file_name_sol.c_str());

    if (file_x.is_open() && file_y.is_open() && file_sol.is_open()) {

        double temp(0);

        if (mesh->get_rank() == 0) {
            //TODO: each process insert a portion!
            for (int j = 0; j != mesh->num_nodes(); ++j) {
                temp = mesh->get_nodes().x(j);
                file_x.write((char*)&temp, sizeof(double));
                temp = mesh->get_nodes().y(j);
                file_y.write((char*)&temp, sizeof(double));
                temp = p_state->get_nodal_value(j);
                file_sol.write((char*)&temp, sizeof(double));
            }
        }

        if (all_values) {
            for (int j = 0; j != p_state->num_sol_values(); ++j) {
                temp = p_state->get_sol(j)[0];
                file_x.write((char*)&temp, sizeof(double));
                temp = p_state->get_sol(j)[1];
                file_y.write((char*)&temp, sizeof(double));
                temp = p_state->get_sol(j)[2];
                file_sol.write((char*)&temp, sizeof(double));
            }
        }

        file_x.close();
        file_y.close();
        file_sol.close();
    }
}

#endif /* OUTPUT_H_ */
