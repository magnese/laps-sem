/*
 * element.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#include "element.h"
int element::c_n(const int& position, const int& corner_index){
  int temp;
  if (corner_index<3){
    if ((corner_position[(corner_index+1)]-corner_position[corner_index]+1)>position)
      temp=neighboring_elements[corner_position[corner_index]+position];
    else temp=-1;
  }
  else {
    if ((num_neighbors()-corner_position[corner_index])>position)
      temp=neighboring_elements[corner_position[corner_index]+position];
    else temp=-1;

  }
  return temp;
}

int element::num_corners(const int& pos){
  int temp;
  if(pos<3) temp=corner_position[pos+1]-corner_position[pos];
  else temp=neighboring_elements.size()-corner_position[3];
  return temp;
}
