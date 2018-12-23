/*
 * mesh.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#include "mesh.h"

typedef enum { NODE,
    ELEMENT,
    GHOST,
    NONE } flag_type;
typedef enum { FIRST,
    SECOND,
    OTHER } flag_temp;
vector<double> null_vec(0);

void mesh::read_file(string file_name)
{

    fstream file;
    string my_str;
    istringstream my_str_converted;

    double x_temp, y_temp;
    vector<int> element_index(4);
    int cont2(0), i_temp, i_temp2, i_temp_idx, i_temp3, i_temp4;
    int cont_ghosts(0);
    flag_type my_flag(NONE);
    flag_temp my_flag2(FIRST);

    line new_line;

    file.open(file_name.c_str());

    if (file.is_open()) {
        while (getline(file, my_str)) {
            if (my_str.compare("$Nodes") == 0)
                my_flag = NODE;
            if (my_str.compare("$EndNodes") == 0) {
                my_flag = NONE;
                my_flag2 = FIRST;
            }
            if (my_str.compare("$Elements") == 0)
                my_flag = ELEMENT;
            if (my_str.compare("$EndElements") == 0) {
                my_flag = NONE;
                my_flag2 = FIRST;
            }

            // add nodes to my_nodes
            if (my_flag == NODE) {
                if (my_flag2 == OTHER) {
                    my_str_converted.str(my_str);
                    my_str_converted >> i_temp;
                    my_nodes.set_i(i_temp, i_temp - 1);
                    my_str_converted >> x_temp >> y_temp;
                    my_nodes.set_x(x_temp, i_temp - 1);
                    my_nodes.set_y(y_temp, i_temp - 1);
                }
                if (my_flag2 == SECOND) {
                    my_str_converted.str(my_str);
                    my_str_converted >> i_temp;
                    my_nodes.resize(i_temp);
                    my_flag2 = OTHER;
                }
                if (my_flag2 == FIRST)
                    my_flag2 = SECOND;
            }

            // add elements to my_elements
            if (my_flag == ELEMENT) {
                if (my_flag2 == OTHER) {
                    my_str_converted.str(my_str);
                    my_str_converted >> i_temp_idx >> i_temp3;
                    if (i_temp3 == 1) {
                        my_str_converted >> i_temp >> i_temp >> i_temp2 >> i_temp;
                        new_line.set_lt(i_temp);
                        new_line.set_elem_idx(i_temp2);
                        my_str_converted >> i_temp >> i_temp2;
                        new_line.set_p(i_temp - 1, i_temp2 - 1);
                        new_line.set_i(i_temp_idx);
                        my_lines.add_line(new_line);
                    }
                    if (i_temp3 == 3) {
                        my_str_converted >> i_temp;
                        for (int i = 0; i != i_temp; ++i) {
                            if (i == 2)
                                my_str_converted >> i_temp4;
                            else
                                my_str_converted >> i_temp2; // put tags in the garbage
                        }
                        for (vector<int>::iterator i = element_index.begin(); i != element_index.end(); ++i) {
                            my_str_converted >> (*i);
                            --*i;
                        }
                        if (i_temp4 < 0) {
                            ++cont_ghosts;
                            i_temp4 = -i_temp4 - 1;
                        }
                        my_elements.set_v(element_index, cont2);
                        my_elements.set_i(i_temp_idx, cont2);
                        my_elements.set_tag(i_temp4, cont2);
                        ++cont2;
                    }
                }
                if (my_flag2 == SECOND) {
                    my_str_converted.str(my_str);
                    my_str_converted >> i_temp;
                    my_elements.resize(i_temp);
                    my_flag2 = OTHER;
                }
                if (my_flag2 == FIRST)
                    my_flag2 = SECOND;
            }
            my_elements.set_num_ghosts(cont_ghosts);
            my_str_converted.clear();
        }
        file.close();
    }
    my_elements.resize(cont2);
}

void mesh::input_metis_conversion(void)
{

    eptr.resize(static_cast<unsigned>(my_elements.num_elements() + my_elements.num_ghosts()) + 1);
    eind.resize(static_cast<unsigned>(my_elements.num_elements() + my_elements.num_ghosts()) * 4);
    for (vector<long int>::iterator i = eptr.begin(); i != eptr.end(); ++i) {
        int pos(i - eptr.begin());
        *i = pos * 4;
        if (pos < (my_elements.num_elements() + my_elements.num_ghosts())) {
            for (int j = 0; j != 4; ++j)
                eind[pos * 4 + j] = static_cast<long int>(my_elements.v(pos)[j]);
        }
    }
}

void mesh::eval_topology()
{

    long int ne(my_elements.num_elements() + my_elements.num_ghosts());
    long int nn(my_nodes.num_nodes());
    long int ncommon(1);
    long int numflag(0);

    for (int j = 0; j != my_elements.num_elements(); ++j) {
        my_elements.resize_neighbors(0, 0, j);
        my_elements.free_pos(j);
    }

    input_metis_conversion();

    METIS_MeshToDual(&ne, &nn, &(eptr[0]), &(eind[0]), &ncommon, &numflag, &xadj_ptr, &adjncy_ptr);

    for (int i = 0; i != my_elements.num_elements(); ++i) {
        for (int j = *(xadj_ptr + i); j != *(xadj_ptr + i + 1); ++j)
            my_elements.add_n(static_cast<int>(*(adjncy_ptr + j)), i);
    }
    METIS_Free(xadj_ptr);
    METIS_Free(adjncy_ptr);
}

void mesh::reorder_neighbors(void)
{
    vector<int> map(0);
    int temp_n;
    vector<int> inc_counter(4);
    vector<int> missing_position(4);
    vector<int> temp_neigh;
    for (int j = 0; j != my_elements.num_elements(); ++j) {
        my_elements.free_pos(j);
        temp_n = 4;
        for (vector<int>::iterator i = inc_counter.begin(); i != inc_counter.end(); ++i) {
            *i = 0;
            missing_position[i - inc_counter.begin()] = 0;
        }
        map.resize(static_cast<unsigned>(my_elements.num_neighbors(j)));
        for (int k = 0; k != my_elements.num_neighbors(j); ++k) {
            map[k] = find_positions_neighbors(j, my_elements.n(j)[k]);
            if (map[k] < 4) {
                --temp_n;
                missing_position[map[k]] = 1;
            } else
                my_elements.inc_edge_position(map[k] - 3, j);
        }
        my_elements.resize_neighbors(my_elements.num_neighbors(j) + temp_n, -1, j);
        temp_neigh.resize(static_cast<int>(my_elements.num_neighbors(j)));
        for (int i = 0; i != 4; ++i) {
            if (missing_position[i] == 0)
                map.push_back(i);
        }
        for (int m = 0; m != (my_elements.num_neighbors(j)); ++m) {
            if (map[m] > 3) {
                ++inc_counter[map[m] - 4];
                map[m] = my_elements.c_pos(j)[map[m] - 4] + inc_counter[map[m] - 4] - 1;
            }
        }
        for (int k = 0; k != (my_elements.num_neighbors(j)); ++k)
            temp_neigh[map[k]] = my_elements.n(j)[k];
        my_elements.set_n(temp_neigh, j);
    }
}

int mesh::find_positions_neighbors(const int& j, const int& k)
{

    vector<int> position(2, 0);
    int temp(0);

    for (int s = 0; (s != 4) && (temp < 2); ++s) {
        for (int t = 0; (t != 4) && (temp < 2); ++t) {
            if ((my_elements.v(j)[s] == my_elements.v(k)[t])) {
                position[temp] = s;
                ++temp;
            }
        }
    }
    if (temp == 1)
        position[0] += 4;
    else {
        if ((position[1] - position[0]) > 1)
            position[0] = position[1];
    }
    return position[0];
}

void mesh::eval_border(void)
{

    int my_rank(0), my_size(0);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &my_size);

    vector<int> temp_values(2 * my_lines.num_lines(), 0);
    vector<int> recv_buffer(0);
    if (my_rank == 0)
        recv_buffer.resize(static_cast<unsigned>(2 * my_lines.num_lines() * my_size), 0);

    for (int j = 0; j != my_elements.num_elements(); ++j) {
        for (int k = 0; k != 4; ++k) {
            if (my_elements.n(j)[k] == -1) {
                for (int l = 0; l != my_lines.num_lines(); ++l) {
                    if (my_elements.v(j)[k] == my_lines.p(l)[0]) {
                        my_elements.set_n((-l - 1), k, j);
                        temp_values[2 * l] = my_elements.i(j);
                        temp_values[2 * l + 1] = k;
                    }
                }
            }
        }
    }

    MPI_Gather(&(temp_values[0]), 2 * my_lines.num_lines(), MPI_INT, &(recv_buffer[0]), 2 * my_lines.num_lines(), MPI_INT, 0, MPI_COMM_WORLD);
    if (my_rank == 0) {
        for (int j = 0; j != (2 * my_lines.num_lines()); ++j) {
            for (int jj = 1; jj != my_size; ++jj)
                temp_values[j] += recv_buffer[j + 2 * my_lines.num_lines() * jj];
        }
    }
    MPI_Bcast(&(temp_values[0]), 2 * my_lines.num_lines(), MPI_INT, 0, MPI_COMM_WORLD);
    for (int j = 0; j != my_lines.num_lines(); ++j) {
        my_lines.set_elem_idx(temp_values[2 * j], j);
        my_lines.set_elem_pos(temp_values[2 * j + 1], j);
    }
}

void mesh::set_border_function(double (*fptr)(double, double, double), const int& starting_line, const int& ending_line, const flag_bc& type, const int& dof)
{
    int idx(border_function_ptr[dof].size());
    if (ending_line < starting_line) {
        for (int j = (starting_line - my_lines.i(0)); j != (my_lines.num_lines() - my_lines.i(0)); ++j)
            my_lines.set_lt(idx, j);
        for (int j = 0; j != (ending_line + 1 - my_lines.i(0)); ++j)
            my_lines.set_lt(idx, j);
    } else {
        for (int j = (starting_line - my_lines.i(0)); j != (ending_line + 1 - my_lines.i(0)); ++j)
            my_lines.set_lt(idx, j);
    }
    border_function_ptr[dof].push_back(fptr);
    lines_bc_tags.push_back(type);
}

void mesh::set_periodicity_mesh(const int& starting_line, const int& ending_line, const int& cor_starting_line, const int& cor_ending_line)
{

    int temp(1);
    int temp_idx(0);
    int temp_cor_idx(0);
    int temp_pos(0);
    int temp_cor_pos(0);
    int num_lines_periodic(ending_line - starting_line + 1);
    int start_pos(starting_line - my_lines.i(0));
    int start_cor_pos(cor_starting_line - my_lines.i(0));
    element new_ghost;

    new_ghost.set_tag(0);
    new_ghost.set_gf();

    set_border_function(&null_fun, starting_line, ending_line, PERIODICITY);
    if (cor_starting_line > cor_ending_line) {
        set_border_function(&null_fun, cor_ending_line, cor_starting_line, PERIODICITY);
        temp = -1;

    } else
        set_border_function(&null_fun, cor_starting_line, cor_ending_line, PERIODICITY);

    periodicity_lines_mapping.resize(static_cast<unsigned>(my_lines.num_lines()), 0);
    for (int j = 0; j != num_lines_periodic; ++j) {
        periodicity_lines_mapping[start_pos + j] = ((start_cor_pos + j * temp) + 1) * temp;
        periodicity_lines_mapping[start_cor_pos + j * temp] = ((start_pos + j) + 1) * temp;
    }

    for (int j = 0; j != num_lines_periodic; ++j) {
        temp_idx = -1;
        temp_cor_idx = -1;
        for (int jj = 0; jj != num_elements(); ++jj) {
            if (my_elements.i(jj) == my_lines.elem_idx(start_pos + j)) {
                temp_idx = my_elements.i(jj);
                temp_pos = jj;
            }
            if (my_elements.i(jj) == my_lines.elem_idx(start_cor_pos + j * temp)) {
                temp_cor_idx = my_lines.elem_idx(start_cor_pos + j * temp);
                temp_cor_pos = jj;
            }
        }
        if (temp_idx > -1 && temp_cor_idx > -1) {
            for (int l = 0; l != 4; ++l) {
                if (-my_elements.n(l, temp_pos) - 1 == (start_pos + j)) {
                    my_elements.set_ncp(my_lines.elem_pos(start_cor_pos + j * temp), l, temp_pos);
                    my_elements.set_n(temp_cor_pos, l, temp_pos);
                }
                if (-my_elements.n(l, temp_cor_pos) - 1 == (start_cor_pos + j * temp)) {
                    my_elements.set_ncp(my_lines.elem_pos(start_pos + j), l, temp_cor_pos);
                    my_elements.set_n(temp_pos, l, temp_cor_pos);
                }
            }
        }
        if (temp_idx > -1 && temp_cor_idx < 0) {
            for (int l = 0; l != 4; ++l) {
                if (-my_elements.n(l, temp_pos) - 1 == (start_pos + j)) {
                    my_elements.set_ncp(my_lines.elem_pos(start_cor_pos + j * temp), l, temp_pos);
                    my_elements.set_n(elements_size(), l, temp_pos);
                }
            }
            new_ghost.set_i(my_lines.elem_idx(start_cor_pos + j * temp));
            my_elements.add_ghost(new_ghost);
        }
    }

    for (int j = 0; j != num_lines_periodic; ++j) {
        temp_idx = -1;
        temp_cor_idx = -1;
        for (int jj = 0; jj != num_elements(); ++jj) {
            if (my_elements.i(jj) == my_lines.elem_idx(start_pos + j)) {
                temp_idx = my_elements.i(jj);
                temp_pos = jj;
            }
            if (my_elements.i(jj) == my_lines.elem_idx(start_cor_pos + j * temp)) {
                temp_cor_idx = my_lines.elem_idx(start_cor_pos + j * temp);
                temp_cor_pos = jj;
            }
        }
        if (temp_idx < 0 && temp_cor_idx > -1) {
            for (int l = 0; l != 4; ++l) {
                if (-my_elements.n(l, temp_cor_pos) - 1 == (start_cor_pos + j * temp)) {
                    my_elements.set_ncp(my_lines.elem_pos(start_pos + j), l, temp_cor_pos);
                    my_elements.set_n(elements_size(), l, temp_cor_pos);
                }
            }
            new_ghost.set_i(my_lines.elem_idx(start_pos + j));
            my_elements.add_ghost(new_ghost);
        }
    }
}

void mesh::write_file(string file_name, const flag_ghost& flag, vector<double>& nodal_values)
{

    ofstream file;
    stringstream my_stringstream;

    file.open(file_name.c_str());

    if (file.is_open()) {
        file << "$MeshFormat\n2 0 8 \n$EndMeshFormat\n$Nodes\n";
        file << my_nodes.num_nodes() << endl;
        for (int i = 0; i != my_nodes.num_nodes(); ++i)
            file << my_nodes.i(i) << " " << my_nodes.x(i) << " " << my_nodes.y(i) << " 0" << endl;
        file << "$EndNodes\n$Elements" << endl;
        if (flag == WRITE_GHOSTS_ON)
            file << my_elements.num_elements() + my_elements.num_ghosts() + my_lines.num_lines() << endl;
        else
            file << my_elements.num_elements() + my_lines.num_lines() << endl;
        for (int j = 0; j != my_lines.num_lines(); ++j)
            file << my_lines.i(j) << " 1 3 0 " << my_lines.elem_idx(j) << " " << my_lines.lt(j) << " " << my_lines.p(j)[0] + 1 << " " << my_lines.p(j)[1] + 1 << endl;
        for (int i = 0; i != my_elements.num_elements(); ++i) {
            file << my_elements.i(i) << " 3 3 0 0 " << my_elements.tag(i) << " " << my_elements.v(i)[0] + 1
                 << " " << my_elements.v(i)[1] + 1 << " " << my_elements.v(i)[2] + 1 << " " << my_elements.v(i)[3] + 1 << endl;
        }
        if (flag == WRITE_GHOSTS_ON) {
            for (int i = my_elements.num_elements(); i != (my_elements.num_elements() + my_elements.num_ghosts()); ++i) {
                file << my_elements.i(i) << " 3 3 0 0 " << (-my_elements.tag(i) - 1) << " " << my_elements.v(i)[0] + 1
                     << " " << my_elements.v(i)[1] + 1 << " " << my_elements.v(i)[2] + 1 << " " << my_elements.v(i)[3] + 1 << endl;
            }
        }
        file << "$EndElements" << endl;
        if (static_cast<int>(nodal_values.size()) == my_nodes.num_nodes()) {
            file << "$NodeData" << endl;
            file << "1" << endl;
            file << "\"Scalar node values\"" << endl;
            file << "1" << endl;
            file << "0.0" << endl;
            file << "3" << endl;
            file << "0" << endl;
            file << "1" << endl;
            file << my_nodes.num_nodes() << endl;
            for (int j = 0; j != my_nodes.num_nodes(); ++j)
                file << my_nodes.i(j) << " " << nodal_values[j] << endl;
            file << "$EndNodeData" << endl;
        }
        file.close();
    }
}

int mesh::corresponding_edge(const int& pos, const int& j)
{
    int temp(-1);
    if (my_elements.n(pos, j) > -1) {
        for (int i = 0; i != 4; ++i) {
            if (my_elements.v(j)[pos] == my_elements.v(my_elements.n(pos, j))[i])
                temp = (i + 3) % 4;
        }
    }
    return temp;
}

void mesh::disp_elements_neighbors(void)
{
    cout << "Elements:" << endl
         << endl;
    for (int j = 0; j != num_elements(); ++j)
        cout << my_elements.i(j) << " " << (my_elements.n(0, j) > -1 ? my_elements.i(my_elements.n(0, j)) : my_elements.n(0, j)) << " "
             << (my_elements.n(1, j) > -1 ? my_elements.i(my_elements.n(1, j)) : my_elements.n(1, j)) << " "
             << (my_elements.n(2, j) > -1 ? my_elements.i(my_elements.n(2, j)) : my_elements.n(2, j)) << " "
             << (my_elements.n(3, j) > -1 ? my_elements.i(my_elements.n(3, j)) : my_elements.n(3, j)) << endl;
    cout << endl;
}

void mesh::disp_ghosts_nodes(void)
{
    cout << "Ghosts:" << endl
         << endl;
    for (int j = num_elements(); j != elements_size(); ++j)
        cout << my_elements.i(j) << " " << my_nodes.i(my_elements.v(j)[0]) << " " << my_nodes.i(my_elements.v(j)[1]) << " " << my_nodes.i(my_elements.v(j)[2]) << " " << my_nodes.i(my_elements.v(j)[3]) << endl;
}

void mesh::set_neighbors_cor_pos(void)
{
    for (int j = 0; j != num_elements(); ++j) {
        for (int jj = 0; jj != 4; ++jj)
            my_elements.set_ncp(corresponding_edge(jj, j), jj, j);
    }
}
