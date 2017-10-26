/*
 * meshsem.cpp
 *
 *  Created on: Feb 20, 2013
 *      Author: magnese
 */

#include "meshsem.h"

void mesh_sem::update_ghosts_deg(void){

  vector<double> deg_elements(basis_deg.begin(),basis_deg.begin()+num_elements());
  get_ghosts_comm()->set_size(1);
  get_ghosts_comm()->set_global_vec(deg_elements);

  basis_deg.resize(static_cast<unsigned>(elements_size()));
  for(int j=0;j!=num_ghosts();++j) basis_deg[j+num_elements()]=static_cast<int>(get_ghosts_comm()->get_ghost_value(j));
}

int mesh_sem::orientation(const int& pos,const int& j){

  int temp(1);
  if(((pos==0)||(pos==3))&&((neighbor_orientation(pos,j)==0)||(neighbor_orientation(pos,j)==1))) temp=-1;
  else{
    if(((pos==1)||(pos==2))&&((neighbor_orientation(pos,j)==2)||(neighbor_orientation(pos,j)==3))) temp=-1;
  }
  return temp;

}
