/*
 * mesh.h
 *
 *  Created on: Dec 12, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef MESH_H_
#define MESH_H_

#include "nodeslist.h"
#include "elementslist.h"
#include "lineslist.h"
#include "mapping.h"
#include "include.h"

extern vector<double> null_vec;

using namespace std;

class mesh{

public:

  mesh():my_nodes(),my_elements(),my_lines(),eptr(0),eind(0),xadj_ptr(NULL),adjncy_ptr(NULL),border_function_ptr(1),lines_bc_tags(0),
  periodicity_lines_mapping(0),my_mapping_ptr(){
    my_mapping_ptr=new mapping<mesh>(this);
  }

  inline void print(void){my_nodes.print();my_elements.print();my_lines.print();}

  void read_file(string file_name);
  void eval_topology();
  void reorder_neighbors(void);
  void eval_border(void);
  void write_file(string file_name,const flag_ghost& flag=WRITE_GHOSTS_ON,vector<double>& nodal_values=null_vec);
  void set_border_function(double (*fptr)(double,double,double),const int& starting_line,const int& ending_line, const flag_bc& type,const int& dof=0);
    void set_periodicity_mesh(const int& starting_line,const int& ending_line,const int& cor_starting_line,const int& cor_ending_line);
    inline void set_simmetry(const int& starting_line,const int& ending_line){set_border_function(&null_fun,starting_line,ending_line,NEUMANN);}

  int corresponding_edge(const int& pos,const int& j);
  int corresponding_node(const int& pos,const int& j);
  void set_neighbors_cor_pos(void);

  inline int num_elements(void){return my_elements.num_elements();}
  inline int num_ghosts(void){return my_elements.num_ghosts();}
  inline int num_nodes(void){return my_nodes.num_nodes();}
  inline int elements_size(void){return num_elements()+num_ghosts();}
  inline int idx_first_global_ele(void){return my_lines.lin(my_lines.num_lines()-1).i()+1;}
  inline double& x(const int& vertex, const int& j){return my_nodes.nod(my_elements.ele(j).v()[vertex]).x();}
  inline double& y(const int& vertex, const int& j){return my_nodes.nod(my_elements.ele(j).v()[vertex]).y();}
  inline flag_bc& get_bc_type(const int& i){return lines_bc_tags[my_lines.lt(i)];}
  inline void resize_bc(const int& ndof){border_function_ptr.resize(static_cast<unsigned>(ndof));}
  inline int& neighbor_orientation(const int& pos,const int& j){return my_elements.ncp(pos,j);}
  inline vector<int>& periodicity_lines_map(void){return periodicity_lines_mapping;}

  inline nodes_list& get_nodes(void){return my_nodes;}
  inline elements_list& get_elements(void){return my_elements;}
  inline lines_list& get_lines(void){return my_lines;}
  inline vector<vector<double (*)(double,double,double)> >& get_bc(void){return border_function_ptr;}
  inline double eval_bc(const int& j, const double& x_value, const double& y_value, const double& t_value,const int& dof=0){return border_function_ptr[dof][my_lines.lt(j)](x_value,y_value,t_value);}
  inline mapping<mesh>*& get_mapping(void){return my_mapping_ptr;}

  void disp_elements_neighbors(void);
  void disp_ghosts_nodes(void);

  ~mesh(void){delete my_mapping_ptr;}

protected:

  void input_metis_conversion(void);
  int find_positions_neighbors(const int& j,const int& k);

  nodes_list my_nodes;
  elements_list my_elements;
  lines_list my_lines;

  vector<long int> eptr;
  vector<long int> eind;
  long int* xadj_ptr;
  long int* adjncy_ptr;

  vector<vector<double (*)(double,double,double)> > border_function_ptr;
  vector<flag_bc> lines_bc_tags;
  vector<int> periodicity_lines_mapping;

  mapping<mesh>* my_mapping_ptr;

private:

};
#endif /* MESH_H_ */
