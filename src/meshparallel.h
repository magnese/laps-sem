/*
 * meshparallel.h
 *
 *  Created on: Dec 3, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef MESHPARALLEL_H_
#define MESHPARALLEL_H_

#include "ghostscommunication.h"
#include "include.h"
#include "mesh.h"

using namespace std;

class mesh_parallel : public mesh {

public:
    mesh_parallel(const int& rank_value, const int& size_value, const int& ncon_ = 1, const int& ncommonnodes_ = 2, const float& ubvecc_ = 1.05)
        : my_rank(rank_value)
        , num_CPUs(size_value)
        , elmdist(0)
        , i_buffer(0)
        , d_buffer(0)
        , ghosts(0)
        , elm_proc(0)
        , temp_dim(0)
        , ubvec(0)
        , tpwgts(0)
        , ncon(ncon_)
        , ncommonnodes(ncommonnodes_)
        , ubvecc(ubvecc_)
        , options(3, 0)
        , partition(0)
        , edgecut(0)
        , temp_ele()
        , send_elements_buffer(0)
        , receive_elements_buffer(0)
        , send_elements_count(0)
        , send_ghosts_buffer(0)
        , receive_ghosts_buffer(0)
        , send_ghosts_count(0)
        , weights(0)
        , ghosts_mapping(0)
        , my_ghosts_comm() {};

    mesh_parallel(const int& rank_value, const int& size_value, const int& n_refinements = 0, string file_name = "mesh.msh", flag_read flag1 = READ_MONO, flag_write flag2 = WRITE_OFF,
        const int& ncon_ = 1, const int& ncommonnodes_ = 2, const float& ubvecc_ = 1.05)
        : my_rank(rank_value)
        , num_CPUs(size_value)
        , elmdist(0)
        , i_buffer(0)
        , d_buffer(0)
        , ghosts(0)
        , elm_proc(0)
        , temp_dim(0)
        , ubvec(0)
        , tpwgts(0)
        , ncon(ncon_)
        , ncommonnodes(ncommonnodes_)
        , ubvecc(ubvecc_)
        , options(3, 0)
        , partition(0)
        , edgecut(0)
        , temp_ele()
        , send_elements_buffer(0)
        , receive_elements_buffer(0)
        , send_elements_count(0)
        , send_ghosts_buffer(0)
        , receive_ghosts_buffer(0)
        , send_ghosts_count(0)
        , weights(0)
        , ghosts_mapping(0)
        , my_ghosts_comm()
    {
        if (flag1 != READ_PARALLEL_PAR) {
            if (flag1 == READ_MONO) {
                if (my_rank == 0) {
                    cout << "START: Single processor reading from file" << file_name << endl;
                    read_file(file_name);
                    cout << "FINISH: Single processor reading from file" << file_name << endl
                         << endl;
                    cout << "START: Single processor evaluation of topology" << endl;
                    eval_topology();
                    cout << "FINISH: Single processor evaluation of topology" << endl
                         << endl;
                }
                if (my_rank == 0)
                    cout << "START: Prepartitioning" << endl;
                pre_partitioning();
                if (my_rank == 0) {
                    cout << "FINISH: Prepartitioning" << endl
                         << endl;
                    cout << "START: Prepartitioning evaluation of topology" << endl;
                }
                eval_topology();
                if (my_rank == 0)
                    cout << "FINISH: Prepartitioning evaluation of topology" << endl
                         << endl;
                if (flag2 == WRITE_PREPARTITIONING || flag2 == WRITE_ALL) {
                    if (my_rank == 0)
                        cout << "START: writing prepartitioning" << endl;
                    write_file(file_name_pre(file_name));
                    if (my_rank == 0)
                        cout << "FINISH: writing prepartitioning" << endl
                             << endl;
                }
            }
            if (flag1 == READ_PARALLEL_PRE) {
                if (my_rank == 0)
                    cout << "START: Prepartitioning reading from file" << endl;
                read_file(file_name_pre(file_name));
                if (my_rank == 0)
                    cout << "FINISH: Prepartitioning reading from file" << endl
                         << endl;
                if (my_rank == 0)
                    cout << "START: Prepartitioning evaluation of topology" << endl;
                eval_topology();
                if (my_rank == 0)
                    cout << "FINISH: Prepartitioning evaluation of topology" << endl
                         << endl;
                compute_elmdist();
            }
            if (my_rank == 0)
                cout << "START: Partitioning" << endl;
            uniform_weights();
            partitioning();
            if (my_rank == 0)
                cout << "FINISH: Partitioning" << endl
                     << endl;
            if (flag2 == WRITE_PARTITIONING || flag2 == WRITE_ALL) {
                if (my_rank == 0)
                    cout << "START: Writing partitioning" << endl;
                write_file(file_name_par(file_name));
                if (my_rank == 0)
                    cout << "FINISH: Writing partitioning" << endl
                         << endl;
            }
        } else {
            if (my_rank == 0)
                cout << "START: Partitioning reading from file" << endl;
            read_file(file_name_par(file_name));
            if (my_rank == 0)
                cout << "FINISH: Partitioning reading from file" << endl
                     << endl;
        }
        compute_elmdist();
        if (my_rank == 0)
            cout << "START: Partitioning evaluation of topology" << endl;
        eval_topology();
        if (my_rank == 0)
            cout << "FINISH: Partitioning evaluation of topology" << endl
                 << endl;
        if (my_rank == 0)
            cout << "START: Partitioning refinements" << endl;
        for (int i = 0; i != n_refinements; ++i) {
            if (my_rank == 0)
                cout << "Refinement number " << i + 1 << endl;
            uniform_weights();
            refinement();
        }
        if (my_rank == 0)
            cout << "FINISH: Partitioning refinements" << endl
                 << endl;
        if (my_rank == 0)
            cout << "START: Partitioning reordering" << endl;
        reorder_neighbors();
        if (my_rank == 0)
            cout << "FINISH: Partitioning reordering" << endl
                 << endl;
        if (flag2 == WRITE_REFINEMENT || flag2 == WRITE_ALL) {
            if (my_rank == 0)
                cout << "START: Writing partitioning refined" << endl;
            write_file(file_name_ref(file_name));
            if (my_rank == 0)
                cout << "FINISH: Writing partitioning refined" << endl
                     << endl;
        }
        my_ghosts_comm = new ghosts_communication<mesh_parallel>(this);
        eval_border();
        set_gfs();
        set_neighbors_cor_pos();
    }

    inline int& get_size(void) { return num_CPUs; }
    inline int& get_rank(void) { return my_rank; }

    string file_name_pre(string file_name);
    string file_name_par(string file_name);
    string file_name_ref(string file_name);

    void pre_partitioning(void);

    inline void free_unused_pre(void)
    {
        i_buffer.resize(0);
        d_buffer.resize(0);
        ghosts.resize(0);
    }
    inline void set_gfs(void)
    {
        for (int j = num_elements(); j != num_elements() + num_ghosts(); ++j)
            my_elements.set_gf(j);
    }
    inline bool& get_gf(const int& j) { return my_elements.get_gf(j); }

    void set_options(const int& ncon_, const int& ncommonnodes_, const float& ubvecc_);
    void partitioning(void);
    void refinement(void);
    inline void free_unused_par(void)
    {
        send_elements_buffer.resize(0);
        receive_elements_buffer.resize(0);
        send_elements_count.resize(0);
        send_ghosts_buffer.resize(0);
        receive_ghosts_buffer.resize(0);
        send_ghosts_count.resize(0);
        ghosts.resize(0);
        elm_proc.resize(0);
    }

    void uniform_weights(void) { weights.resize(static_cast<unsigned>(my_elements.num_elements() * ncon), 1); }

    inline int& get_ghosts_rank(const int& j) { return ghosts_mapping[j * 2]; }
    inline int& get_ghosts_pos(const int& j) { return ghosts_mapping[j * 2 + 1]; }

    inline vector<long int>& get_elmdist(void) { return elmdist; }
    inline ghosts_communication<mesh_parallel>*& get_ghosts_comm(void) { return my_ghosts_comm; }

    void set_periodicity(const int& starting_line, const int& ending_line, const int& cor_starting_line, const int& cor_ending_line);

    ~mesh_parallel(void) { delete my_ghosts_comm; }

protected:
    inline int find_pos_elmdist(const int& value)
    {
        vector<long int>::iterator iter(elmdist.begin() + 1);
        for (; iter != elmdist.end(); ++iter) {
            if (*iter > value)
                break;
        }
        return iter - elmdist.begin() - 1;
    }

    void fill_buffer_elements(const int& index);
    void fill_buffer_ghosts(const int& index);
    void unpack_buffer_elements(const int& pos);
    void fill_buffer_nodes(void);
    void unpack_buffer_nodes(void);
    void fill_buffer_lines(void);
    void unpack_buffer_lines(void);
    void search_ghosts(void);

    void create_ghosts_mapping(void);

    void compute_elmdist(void);
    void input_parmetis_conversion(void);
    void set_ghosts_elm_proc(void);

    void partitioning_move_elements(void);
    void create_par_send_buffers(void);
    void unpack_par_elements_buffer(const int& num_new_elements);
    void unpack_par_ghosts_buffer(void);

    int my_rank;
    int num_CPUs;

    vector<long int> elmdist;

    vector<int> i_buffer;
    vector<double> d_buffer;
    vector<vector<int>> ghosts;
    vector<vector<int>> elm_proc;
    int temp_dim;

    vector<float> ubvec;
    vector<float> tpwgts;
    long int ncon;
    long int ncommonnodes;
    float ubvecc;
    vector<long int> options;
    vector<long int> partition;
    long int edgecut;

    element temp_ele;

    vector<int> send_elements_buffer;
    vector<int> receive_elements_buffer;
    vector<int> send_elements_count;
    vector<int> send_ghosts_buffer;
    vector<int> receive_ghosts_buffer;
    vector<int> send_ghosts_count;
    vector<long int> weights;

    vector<int> ghosts_mapping;

    ghosts_communication<mesh_parallel>* my_ghosts_comm;
};

#endif /* MESHPREPARTITIONING_H_ */
