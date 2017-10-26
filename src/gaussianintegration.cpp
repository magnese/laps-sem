/*
 * gaussianintegration.cpp
 *
 *  Created on: Oct 23, 2012
 *      Author: Marco Agnese
 */

#include "gaussianintegration.h"

void gaussian_integration::set_order(const int& order){
  int old_size(nodes.size());
  int new_size(necessary_size(order)-1);
  if(new_size>old_size){
    weights.resize(static_cast<unsigned>(new_size));
    nodes.resize(static_cast<unsigned>(new_size));
    leg_pol.set_order(new_size+1);
    double *dr=NULL;
    for(int i=(old_size+1);i!=(new_size+1);++i){
      nodes[i-1].resize(i+1);
      *(nodes[i-1].begin())=-1;
      *(nodes[i-1].end()-1)=1;
      dr=new double[2*i-2];
      gsl_poly_complex_workspace *w=gsl_poly_complex_workspace_alloc(i);
      gsl_poly_complex_solve(&(leg_pol.c_der_reverse(i)[0]),i,w,dr);
      gsl_poly_complex_workspace_free(w);
      for(int j=0;j!=(i-1);++j) nodes[i-1][j+1]=dr[2*j];
      delete []dr;
      dr=NULL;
      sort(nodes[i-1].begin(),nodes[i-1].end());
      weights[i-1].resize(i+1);
      double temp=2.0/(i*(i+1));
      for(int j=0;j!=(i+1);++j) weights[i-1][j]=temp/pow(leg_pol.eval(nodes[i-1][j],i),2);
    }
  }
}

void gaussian_integration::disp(void){;
  cout<<"Nodes and weights of Gauss-Legendre-Lobatto quadrature formula"<<endl<<endl;
  for(int i=0;i!=static_cast<int>(nodes.size());++i){
    cout<<"degree : "<<i+1<<endl;
    for(int j=0;j!=(i+2);++j){
      cout<<"\t"<<nodes[i][j]<<" \t"<<weights[i][j]<<endl;
    }
    cout<<endl;
  }
}

