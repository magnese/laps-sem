/*
 * legmodalbasis.cpp
 *
 *  Created on: Oct 15, 2012
 *      Author: Marco Agnese
 */

#include "legmodalbasis.h"

void leg_modal_basis::set_basis_order(const int& j){
  if(max_order()<j){
    leg_pol.set_order(j);
    int old_size(max_order()+1);
    double temp(0);
    basis_coeff.resize(static_cast<unsigned>(j+1));
    for(int i=old_size;i!=(j+1);++i){
      basis_coeff[i]=leg_pol.c(i);
      temp=-pow(static_cast<double>(4*i-2),-0.5);
      for(vector<double>::iterator k=basis_coeff[i].begin();k!=basis_coeff[i].end();++k){
        (*k)*=temp;
        int position(k-basis_coeff[i].begin()-2);
        if(position>=0) (*k)-=temp*leg_pol.c(i-2)[position];
      }
    }
  }
}
