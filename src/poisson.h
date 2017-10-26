/*
 * poisson.h
 *
 *  Created on: Mar 8, 2013
 *      Author: magnese
 */

#ifndef POISSON_H_
#define POISSON_H_

#include "linearequation.h"

using namespace std;

template<class M,class B, class G>
class poisson:public linear_equation<M,B,G> {

public:

  poisson(M* mesh_ptr,double (*fptr)(double,double)):linear_equation<M,B,G>(mesh_ptr),source(fptr){}

  ~poisson(void){};

private:

  double (*source)(double,double);

  double eval_LHS(const int& k,const int& kk,const int& i,const int& ii,const int& j,const int& jj,const int& ele_pos);
  double eval_RHS(const int& k,const int& kk,const int& i,const int& ii,const int& ele_pos);
  double eval_neumann(const double fun_value,const int& k,const int& i,const int& ele_pos);
};

template<class M,class B, class G>
double poisson<M,B,G>::eval_LHS(const int& k,const int& kk,const int& i,const int& ii,const int& j,const int& jj,const int& ele_pos){

  return this->der_phi(k,j)*this->phi(kk,jj)*((pow(this->J(0),2)+pow(this->J(1),2))*this->der_phi(k,i)*this->phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->phi(k,i)*this->der_phi(kk,ii))+
      this->phi(k,j)*this->der_phi(kk,jj)*((pow(this->J(2),2)+pow(this->J(3),2))*this->phi(k,i)*this->der_phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->der_phi(k,i)*this->phi(kk,ii));

}

template<class M,class B, class G>
double poisson<M,B,G>::eval_RHS(const int& k,const int& kk,const int& i,const int& ii,const int& ele_pos){

  return source(this->x_physical(k,kk,ele_pos),this->y_physical(k,kk,ele_pos))*this->phi(k,i)*this->phi(kk,ii);
}

template<class M,class B, class G>
double poisson<M,B,G>::eval_neumann(const double fun_value,const int& k,const int& i,const int& ele_pos){
  return fun_value*this->phi(k,i);
}

#endif /* POISSON_H_ */
