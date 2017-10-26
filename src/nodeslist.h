/*
 * nodeslist.h
 *
 *  Created on: Sep 27, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef NODESLIST_H_
#define NODESLIST_H_

#include "node.h"
#include "include.h"

using namespace std;

class nodes_list{

public:

  nodes_list(void):my_nodes(0){}

  inline int long num_nodes(void){return static_cast<int>(my_nodes.size());}
  inline double& x(const int& j){return my_nodes[j].x();}
  inline double& y(const int& j){return my_nodes[j].y();}
  inline int& i(const int& j){return my_nodes[j].i();}
  inline node& nod(const int& i){return my_nodes[i];}

  inline void set_x(double& value,const int& j){my_nodes[j].set_x(value);}
  inline void set_y(double& value,const int& j){my_nodes[j].set_y(value);}
  inline void set_i(const int& value,const int& j){my_nodes[j].set_i(value);}

  inline void add_node(node& new_node){my_nodes.push_back(new_node);}
  inline void resize(const int& value){my_nodes.resize(static_cast<unsigned>(value));}

  void print(void);

private:

  vector<node> my_nodes;

};


#endif /* NODESLIST_H_ */
