/*
 * node.h
 *
 *  Created on: Oct 16, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef NODE_H_
#define NODE_H_

#include "include.h"

using namespace std;

class node{

public:

  node(void):coordinates(2,0),global_index(0){}

  inline int& i(void){return global_index;}
  inline double& x(void){return coordinates[0];}
  inline double& y(void){return coordinates[1];}

  inline void set_x(double& value){coordinates[0]=value;}
  inline void set_y(double& value){coordinates[1]=value;}
  inline void set_i(const int& value){global_index=value;}

  void print(void);

private:

  vector<double> coordinates;
  int global_index;

};

#endif /* NODE_H_ */
