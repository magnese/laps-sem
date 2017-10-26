/*
 * line.h
 *
 *  Created on: Oct 17, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef LINE_H_
#define LINE_H_

#include "include.h"

using namespace std;

class line{

public:

line(void):points(2,0),global_index(0),line_tag(0),element_idx(0),element_pos(0){}

inline vector<int>& p(void){return points;}
inline int& i(void){return global_index;}
inline int& lt(void){return line_tag;}
inline int& elem_idx(void){return element_idx;}
inline int& elem_pos(void){return element_pos;}

inline void set_p(const int& point_a, const int& point_b){points[0]=point_a;points[1]=point_b;}
inline void set_i(const int& value){global_index=value;}
inline void set_lt(const int& value){line_tag=value;}
inline void set_elem_idx(const int& value){element_idx=value;}
inline void set_elem_pos(const int& value){element_pos=value;}

private:

vector <int> points;
int global_index;
int line_tag;
int element_idx;
int element_pos;

};

#endif /* LINE_H_ */
