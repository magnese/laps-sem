/*
 * elementslist.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#include "elementslist.h"

void elements_list::copy(const int& l,const int& j){
  if(l!=j){
    set_v(v(l),j);
    set_n(n(l),j);
    set_i(i(l),j);
  }
}

void elements_list::print(void){
  cout<<endl<<"element list"<<endl<<endl;
  for(vector<element>::iterator j=my_elements.begin();j!=my_elements.end();++j){
    cout<<j->i()<<"\t"<<j->v()[0]+1<<" "<<j->v()[1]+1<<" "<<j->v()[2]+1<<" "<<j->v()[3]+1<<"\t\t";
    for(int l=0;l!=j->num_neighbors();++l) cout<<(j->n()[l]>-1?i(j->n()[l]):-1)<<" ";
    cout<<"\t\t";
    for(int l=0;l!=4;++l) cout<<j->c_pos()[l]<<" ";
    cout<<endl;
  }
  cout<<endl;
}


