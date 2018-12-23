/*
 * mapping.h
 *
 *  Created on: Feb 20, 2013
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef MAPPING_H_
#define MAPPING_H_

#include "include.h"
#include "mappingbase.h"

using namespace std;
template <class M>
class mapping : public mapping_base {
public:
    mapping(M* mesh_ptr)
        : my_mesh(mesh_ptr) {};

    inline double& Fx(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Fxb(x_logic, y_logic, my_mesh->x(0, j), my_mesh->x(1, j), my_mesh->x(2, j), my_mesh->x(3, j));
        return temp;
    }
    inline double& Fy(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Fyb(x_logic, y_logic, my_mesh->y(0, j), my_mesh->y(1, j), my_mesh->y(2, j), my_mesh->y(3, j));
        return temp;
    }

    inline double& Jxx(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jxxb(x_logic, y_logic, my_mesh->x(0, j), my_mesh->x(1, j), my_mesh->x(2, j), my_mesh->x(3, j));
        return temp;
    }
    inline double& Jxy(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jxyb(x_logic, y_logic, my_mesh->x(0, j), my_mesh->x(1, j), my_mesh->x(2, j), my_mesh->x(3, j));
        return temp;
    }
    inline double& Jyx(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jyxb(x_logic, y_logic, my_mesh->y(0, j), my_mesh->y(1, j), my_mesh->y(2, j), my_mesh->y(3, j));
        return temp;
    }
    inline double& Jyy(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jyyb(x_logic, y_logic, my_mesh->y(0, j), my_mesh->y(1, j), my_mesh->y(2, j), my_mesh->y(3, j));
        return temp;
    }

    inline double& det_J(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jxx(x_logic, y_logic, j) * Jyy(x_logic, y_logic, j) - Jxy(x_logic, y_logic, j) * Jyx(x_logic, y_logic, j);
        return temp;
    }
    inline double& Jxx_1(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jyy(x_logic, y_logic, j) / det_J(x_logic, y_logic, j);
        return temp;
    }
    inline double& Jxy_1(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = -Jxy(x_logic, y_logic, j) / det_J(x_logic, y_logic, j);
        return temp;
    }
    inline double& Jyx_1(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = -Jyx(x_logic, y_logic, j) / det_J(x_logic, y_logic, j);
        return temp;
    }
    inline double& Jyy_1(const double& x_logic, const double& y_logic, const int& j)
    {
        temp = Jxx(x_logic, y_logic, j) / det_J(x_logic, y_logic, j);
        return temp;
    }

private:
    M* my_mesh;
};

#endif /* MAPPING_H_ */
