/*
 * utilities.h
 *
 *  Created on: Nov 13, 2012
 *      Author: Marco Agnese
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "include.h"

#include "/nh/research/zehuag/GSL/include/gsl/gsl_vector.h"
#include "/nh/research/zehuag/GSL/include/gsl/gsl_matrix.h"
#include "/nh/research/zehuag/GSL/include/gsl/gsl_poly.h"
#include "/nh/research/zehuag/GSL/include/gsl/gsl_linalg.h"
#include "/nh/research/zehuag/GSL/include/gsl/gsl_mode.h"
#include "/nh/research/zehuag/GSL/include/gsl/gsl_permutation.h"

/*
#include "/nh/research/zehuag/petsc-dev/include/petscpc.h"
#include "/nh/research/zehuag/petsc-dev/include/petscvec.h"
#include "/nh/research/zehuag/petsc-dev/include/petscis.h"
#include "/nh/research/zehuag/petsc-dev/include/petscmat.h"
#include "/nh/research/zehuag/petsc-dev/include/petscksp.h"
#include "/nh/research/zehuag/petsc-dev/include/petsc-private/vecimpl.h"
#include "/nh/research/zehuag/petsc-dev/include/petscdm.h"
#include "/nh/research/zehuag/petsc-dev/include/petscsys.h"
#include "/nh/research/zehuag/petsc-dev/include/petscsf.h"
#include "/nh/research/zehuag/petsc-dev/include/petscbag.h"
*/

#include "petscpc.h"
#include "petscvec.h"
#include "petscis.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petsc-private/vecimpl.h"
#include "petscdm.h"
#include "petscsys.h"
#include "petscsf.h"
#include "petscbag.h"

typedef enum{STATIC_CONDENSATION_ON,STATIC_CONDENSATION_OFF} static_condesation_flag;
typedef enum{LINEAR_EQUATION,NON_LINEAR_EQUATION} eq_type;

using namespace std;

// routines for polynomials (coefficients decreasing order)

vector<double> poly_reverse_coeff(const vector<double>& coeff);
double poly_val(const vector<double>& coeff,const double& x);
void poly_der(const vector<double>& coeff,vector<double>& der_coeff);
int poly_deg(const vector<double>& coeff);

// display routines

template<typename T>
void vec_disp(const vector<T>& vec){
  typename vector<T>::const_iterator j;
  for(j=vec.begin();j!=vec.end();++j) cout<<*j<<" ";
  cout<<endl;
}

template<typename T>
void vec_disp_formatted(const vector<T>& vec,const int& c){
  typename vector<T>::const_iterator j;
  for(j=vec.begin();j!=vec.end();++j){
    if((j-vec.begin())%c==0) cout<<endl<<((j-vec.begin())/c)<<": ";
    cout<<*j<<" ";
  }
  cout<<endl;
}

template<typename T,class M>
void mesh_vec_disp(const vector<T>& vec,const int& c,M* mesh_ptr){
  typename vector<T>::const_iterator j;
  for(j=vec.begin();j!=vec.end();++j){
    if((j-vec.begin())%c==0) cout<<endl<<mesh_ptr->i((j-vec.begin())/c)<<": ";
    cout<<*j<<" ";
  }
  cout<<endl;
}

template<typename T,class M>
void ghost_vec_disp(const vector<T>& vec,const int& c,M* mesh_ptr){
  typename vector<T>::const_iterator j;
  for(j=vec.begin();j!=vec.end();++j){
    if((j-vec.begin())%c==0) cout<<endl<<mesh_ptr->i((j-vec.begin())/c+mesh_ptr->num_elements())<<": ";
    cout<<*j<<" ";
  }
  cout<<endl;
}

template<typename T,class M>
void node_vec_disp(const vector<T>& vec,const int& c,M* mesh_ptr){
  typename vector<T>::const_iterator j;
  for(j=vec.begin();j!=vec.end();++j){
    if((j-vec.begin())%c==0) cout<<" || "<<mesh_ptr->get_nodes().i((j-vec.begin())/c)<<": ";
    cout<<*j<<" ";
  }
  cout<<endl;
}

template<typename T>
void mat_disp(const vector<vector<T> >& mat){
  typename vector<vector<T> >::const_iterator j;
  for(j=mat.begin();j!=mat.end();++j) vec_disp<T>(*j);
}

// miscellaneous

int min_positive(const int& x,const int& y,const int& z);

// debugging

template<class M>
void breakpoint(const int& number,M* mesh_ptr){
  if(mesh_ptr->get_rank()==0) cout<<endl<<"Breakpoint "<<number<<":"<<endl;
  for(int j=0;j!=mesh_ptr->get_size();++j){
    if(mesh_ptr->get_rank()==j) cout<<"status process "<<j<<" OK"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(mesh_ptr->get_rank()==0) cout<<endl;
}

template<typename T,class M>
void parallel_value_disp(const T& value,M* mesh_ptr,const int& number=0){
  if(mesh_ptr->get_rank()==0) cout<<endl<<"Output "<<number<<":"<<endl;
  for(int j=0;j!=mesh_ptr->get_size();++j){
    if(mesh_ptr->get_rank()==j){ cout<<"value = "<<value<<" on process "<<j<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(mesh_ptr->get_rank()==0) cout<<endl;
}

template<typename T,class M>
void parallel_vec_disp(const vector<T>& vec,M* mesh_ptr,const int& number=0){
  if(mesh_ptr->get_rank()==0) cout<<endl<<"Output "<<number<<":"<<endl;
  for(int j=0;j!=mesh_ptr->get_size();++j){
    if(mesh_ptr->get_rank()==j){ cout<<"process: "<<j<<endl;
    vec_disp<T>(vec);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if(mesh_ptr->get_rank()==0) cout<<endl;
}

#endif /* UTILITIES_H_ */



