/*
 * mappingbase.h
 *
 *  Created on: Feb 26, 2013
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef MAPPINGBASE_H_
#define MAPPINGBASE_H_

#include "include.h"

using namespace std;

class mapping_base{

public:
  mapping_base(void):temp(0){};

  inline double& Fxb(const double& x_logic,const double& y_logic,const double& x0,const double& x1,const double& x2,const double& x3){
    temp=0.25*((x0*(1-x_logic)*(1+y_logic))+x1*(1+x_logic)*(1+y_logic)+x2*(1+x_logic)*(1-y_logic)+x3*(1-x_logic)*(1-y_logic));
    return temp;}
  inline double& Fyb(const double& x_logic,const double& y_logic,const double& y0,const double& y1,const double& y2,const double& y3){
    temp=0.25*((y0*(1-x_logic)*(1+y_logic))+y1*(1+x_logic)*(1+y_logic)+y2*(1+x_logic)*(1-y_logic)+y3*(1-x_logic)*(1-y_logic));
    return temp;}

  inline double& Jxxb(const double& x_logic,const double& y_logic,const double& x0,const double& x1,const double& x2,const double& x3){
    temp=0.25*((-x0*(1+y_logic))+x1*(1+y_logic)+x2*(1-y_logic)-x3*(1-y_logic));
    return temp;}
  inline double& Jxyb(const double& x_logic,const double& y_logic,const double& x0,const double& x1,const double& x2,const double& x3){
    temp=0.25*((x0*(1-x_logic))+x1*(1+x_logic)-x2*(1+x_logic)-x3*(1-x_logic));
    return temp;}
  inline double& Jyxb(const double& x_logic,const double& y_logic,const double& y0,const double& y1,const double& y2,const double& y3){
    temp=0.25*((-y0*(1+y_logic))+y1*(1+y_logic)+y2*(1-y_logic)-y3*(1-y_logic));
    return temp;}
  inline double& Jyyb(const double& x_logic,const double& y_logic,const double& y0,const double& y1,const double& y2,const double& y3){
    temp=0.25*((y0*(1-x_logic))+y1*(1+x_logic)-y2*(1+x_logic)-y3*(1-x_logic));
    return temp;}

protected:

  double temp;

private:

};

#endif /* MAPPINGBASE_H_ */
