/*
 * lineslist.h
 *
 *  Created on: Oct 17, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef LINESLIST_H_
#define LINESLIST_H_

#include "line.h"
#include "include.h"

using namespace std;

class lines_list{

public:

lines_list(void):my_lines(0){}

inline int num_lines(void){return static_cast<int>(my_lines.size());}
inline vector<int>& p(const int& j){return my_lines[j].p();}

inline int& i(const int& j){return my_lines[j].i();}
inline int& lt(const int& j){return my_lines[j].lt();}
inline line& lin(const int& i){return my_lines[i];}
inline int& elem_idx(const int& j){return my_lines[j].elem_idx();}
inline int& elem_pos(const int& j){return my_lines[j].elem_pos();}

inline void set_p(const int& point_a, const int& point_b,const  int& j){my_lines[j].set_p(point_a,point_b);}
inline void set_i(const int& value,const int& j){my_lines[j].set_i(value);}
inline void set_lt(const int& value,const int& j){my_lines[j].set_lt(value);}
inline void set_elem_idx(const int& value,const int& j){my_lines[j].set_elem_idx(value);}
inline void set_elem_pos(const int& value,const int& j){my_lines[j].set_elem_pos(value);}

inline void resize(const int& value){my_lines.resize(static_cast<unsigned>(value));}
void print(void);

inline void add_line(line& new_line){my_lines.push_back(new_line);}

private:

vector<line> my_lines;

};

#endif /* LINESLIST_H_ */
