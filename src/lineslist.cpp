/*
 * lineslist.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */
#include "lineslist.h"

void lines_list::print(void){
  cout<<endl<<"lines list"<<endl<<endl;
  for(vector<line>::iterator j=my_lines.begin();j!=my_lines.end();++j)
    cout<<j->i()<<" "<<j->p()[0]+1<<" "<<j->p()[1]+1<<" "<<j->lt()<<endl;
  cout<<endl;
}
