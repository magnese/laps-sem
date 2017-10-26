/*
 * nodeslist.cpp
 *
 *  Created on: Sep 27, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#include "nodeslist.h"

void nodes_list::print(void){
  cout<<endl<<"node list"<<endl<<endl;
  for(vector<node>::iterator i=my_nodes.begin();i!=my_nodes.end();++i) i->print();
  cout<<endl;
}
