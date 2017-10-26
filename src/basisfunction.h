/*
 * basisfunction.h
 *
 *  Created on: Nov 30, 2012
 *      Author: Marco Agnese
 */

#ifndef BASISFUNCTION_H_
#define BASISFUNCTION_H_

#include "utilities.h"

using namespace std;

class basis_function {

public:

  basis_function(void):basis_coeff(0),der_basis_coeff(0){}

  inline void set_order(const int& j){set_basis_order(j);set_der_order();}
  inline int max_order(void){return static_cast<int>(basis_coeff.size()-1);}


  inline double eval(double& x,const int& j){return poly_val(basis_coeff[j],x);}
  inline double eval_der(double& x,const int& j){return poly_val(der_basis_coeff[j],x);}

  inline void disp(void){
    cout<<endl<<"Coefficients of the basis in decreasing order"<<endl<<endl;mat_disp<double>(basis_coeff);cout<<endl;}
  inline void disp_der(void){
    cout<<endl<<"Coefficients of the first derivative of the basis in decreasing order"<<endl<<endl;mat_disp<double>(der_basis_coeff);cout<<endl;}

  virtual ~basis_function(){}

protected:

    vector<vector<double> > basis_coeff;
  vector<vector<double> > der_basis_coeff;

  virtual void set_basis_order(const int& j)=0;

private:

  void set_der_order(void);

};

#endif /* BASISFUNCTION_H_ */
