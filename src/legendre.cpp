/*
 * legendre.cpp
 *
 *  Created on: Oct 15, 2012
 *      Author: Marco Agnese
 */

#include "legendre.h"

using namespace std;

void legendre::set_order(const int& new_leg_order){
  int old_leg_order(max_order());
  if(old_leg_order<new_leg_order){
    double temp1(0),temp2(0);
    leg_coeff.resize(static_cast<unsigned>(new_leg_order+1));
    leg_der_coeff.resize(static_cast<unsigned>(new_leg_order+1));
    for(int i=(old_leg_order+1);i!=(new_leg_order+1);++i){
      leg_coeff[i]=leg_coeff[i-1];
      leg_coeff[i].resize(static_cast<unsigned>(i+1));
      temp1=(2*i-1)/(static_cast<double>(i));
      temp2=(i-1)/(static_cast<double>(i));
      for(vector<double>::iterator j=leg_coeff[i].begin();j!=leg_coeff[i].end();++j){
        (*j)*=temp1;
        int position=(j-leg_coeff[i].begin())-2;
        if(position>=0) (*j)-=temp2*leg_coeff[i-2][position];
      }
      poly_der(leg_coeff[i],leg_der_coeff[i]);
    }
  }
  set_reverse_order();
}

void legendre::set_reverse_order(void){
  int temp_dim(leg_coeff_reverse.size());
  if(temp_dim<static_cast<int>(leg_coeff.size())){
    leg_coeff_reverse.resize(leg_coeff.size());
    leg_der_coeff_reverse.resize(leg_der_coeff.size());
    for(int i=temp_dim;i!=static_cast<int>(leg_coeff.size());++i){
      leg_coeff_reverse[i]=poly_reverse_coeff(leg_coeff[i]);
      leg_der_coeff_reverse[i]=poly_reverse_coeff(leg_der_coeff[i]);
    }
  }
}
