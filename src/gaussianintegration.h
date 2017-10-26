/*
 * gaussianintegration.h
 *
 *  Created on: Oct 23, 2012
 *      Author: Marco Agnese
 */

#ifndef GAUSSIANINTEGRATION_H_
#define GAUSSIANINTEGRATION_H_

#include "legendre.h"
#include "utilities.h"

using namespace std;

class gaussian_integration{

public:

  gaussian_integration(const int& order=3):leg_pol(2),nodes(2),weights(2){
    nodes[0].resize(2);
    nodes[0][0]=-1;
    nodes[0][1]=1;
    nodes[1].resize(3);
    nodes[1][0]=-1;
    nodes[1][1]=0;
    nodes[1][2]=1;
    weights[0].resize(2);
    weights[0][0]=1;
    weights[0][1]=1;
    weights[1].resize(3);
    weights[1][0]=1.0/3.0;
    weights[1][1]=4.0/3.0;
    weights[1][2]=1.0/3.0;
    set_order(order);
  }

  void set_order(const int& order);
  inline int max_order(void){return static_cast<int>((nodes.size()+1)*2-3);}
  inline int necessary_size(const int& order){return ceil((order+3)/2.0);}

  inline vector<double>& get_nodes(const int& order){return nodes[necessary_size(order)-2];}
  inline vector<double>& get_weights(const int& order){return weights[necessary_size(order)-2];}

  inline double& get_n(const int& idx, const int& pos){return nodes[pos][idx];}
  inline double& get_w(const int& idx, const int& pos){return weights[pos][idx];}

  void disp(void);

private:

  legendre leg_pol;
  vector<vector<double> > nodes;
  vector<vector<double> > weights;

};

#endif /* GAUSSIANINTEGRATION_H_ */
