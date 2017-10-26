/*
 * lineargradshafranov.h
 *
 *  Created on: July 2, 2013
 *      Author: Marco Agnese
 */

#ifndef LINEARGRADSHAFRANOVMOD_H_
#define LINEARGRADSHAFRANOVMOD_H_

#include "linearequation.h"

using namespace std;

template<class M,class B, class G>
class linear_grad_shafranov_mod:public linear_equation<M,B,G> {

public:

  linear_grad_shafranov_mod(M* mesh_ptr,double p_0_value,double u_0_value,double (*fptr)(double,double)=null_f):
    linear_equation<M,B,G>(mesh_ptr),u_0(u_0_value),p_0(p_0_value),source(fptr){}

  ~linear_grad_shafranov_mod(void){};

private:

  double u_0;
  double p_0;
  double (*source)(double,double);

  double eval_LHS(const int& k,const int& kk,const int& i,const int& ii,const int& j,const int& jj,const int& ele_pos);
  double eval_RHS(const int& k,const int& kk,const int& i,const int& ii,const int& ele_pos);
  double eval_neumann(const double fun_value,const int& k,const int& i,const int& ele_pos);

};

template<class M,class B, class G>
double linear_grad_shafranov_mod<M,B,G>::eval_LHS(const int& k,const int& kk,const int& i,const int& ii,const int& j,const int& jj,const int& ele_pos){

  return (this->der_phi(k,j)*this->phi(kk,jj)*((pow(this->J(0),2)+pow(this->J(1),2))*this->der_phi(k,i)*this->phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->phi(k,i)*this->der_phi(kk,ii))+
      this->phi(k,j)*this->der_phi(kk,jj)*((pow(this->J(2),2)+pow(this->J(3),2))*this->phi(k,i)*this->der_phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->der_phi(k,i)*this->phi(kk,ii)))*this->y_physical(k,kk,ele_pos)
      +2*this->phi(k,i)*this->phi(kk,ii)*(this->der_phi(k,j)*this->phi(kk,jj)*this->J(1)+this->phi(k,j)*this->der_phi(kk,jj)*this->J(3))-
      2*p_0/u_0*this->phi(k,i)*this->phi(kk,ii)*this->phi(k,j)*this->phi(kk,jj);
}

template<class M,class B, class G>
double linear_grad_shafranov_mod<M,B,G>::eval_RHS(const int& k,const int& kk,const int& i,const int& ii,const int& ele_pos){

  return (source(this->x_physical(k,kk,ele_pos),this->y_physical(k,kk,ele_pos))*this->phi(k,i)*this->phi(kk,ii))*pow(this->y_physical(k,kk,ele_pos),3);
}

template<class M,class B, class G>
double linear_grad_shafranov_mod<M,B,G>::eval_neumann(const double fun_value,const int& k,const int& i,const int& ele_pos){
  return 0;
}

#endif /* LINEARGRADSHAFRANOV_H_ */
