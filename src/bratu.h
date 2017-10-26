/*
 * bratu.h
 *
 *  Created on: Jun 27, 2013
 *      Author: Marco Agnese
 */

#ifndef BRATU_H_
#define BRATU_H_

#include "nonlinearequation.h"

using namespace std;

template<class M,class B,class G>
class bratu:public non_linear_equation<M,B,G> {

public:

  bratu(M* mesh_ptr,double toll,double lambda_value,double (*fptr)(double,double)=null_f):non_linear_equation<M,B,G>(mesh_ptr,
    toll,fptr),lambda(lambda_value){}

  ~bratu(void){};

private:

  double lambda;

  double eval_Jacobian(const int& k,const int& kk,const int& i,const int& ii,const int& j,const int& jj,const int& ele_pos);
  double eval_F0(const int& k,const int& kk,const int& i,const int& ii,const int& ele_pos);
  double eval_neumann(const double fun_value,const int& k,const int& i,const int& ele_pos);

};

template<class M,class B, class G>
double bratu<M,B,G>::eval_Jacobian(const int& k,const int& kk,const int& i,const int& ii,const int& j,const int& jj,const int& ele_pos){

  return (this->der_phi(k,j)*this->phi(kk,jj)*((pow(this->J(0),2)+pow(this->J(1),2))*this->der_phi(k,i)*this->phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->phi(k,i)*this->der_phi(kk,ii))+
      this->phi(k,j)*this->der_phi(kk,jj)*((pow(this->J(2),2)+pow(this->J(3),2))*this->phi(k,i)*this->der_phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->der_phi(k,i)*this->phi(kk,ii)))
      -lambda*this->phi(k,j)*this->phi(kk,jj)*this->phi(k,i)*this->phi(kk,ii)*exp(this->sol());
}

template<class M,class B,class G>
double bratu<M,B,G>::eval_F0(const int& k,const int& kk,const int& i,const int& ii,const int& ele_pos){

  return (this->sol_x()*((pow(this->J(0),2)+pow(this->J(1),2))*this->der_phi(k,i)*this->phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->phi(k,i)*this->der_phi(kk,ii))+
      this->sol_y()*((pow(this->J(2),2)+pow(this->J(3),2))*this->phi(k,i)*this->der_phi(kk,ii)+
      (this->J(0)*this->J(2)+this->J(1)*this->J(3))*this->der_phi(k,i)*this->phi(kk,ii)))
      -lambda*this->phi(k,i)*this->phi(kk,ii)*exp(this->sol());
}

template<class M,class B, class G>
double bratu<M,B,G>::eval_neumann(const double fun_value,const int& k,const int& i,const int& ele_pos){
  return 0; //TODO
}

#endif /* BRATU_H_ */



