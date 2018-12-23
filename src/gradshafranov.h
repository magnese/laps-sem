/*
 * gradshafranov.h
 *
 *  Created on: Jun 4, 2013
 *      Author: Marco Agnese
 */

#ifndef GRADSHAFRANOV_H_
#define GRADSHAFRANOV_H_

#include "nonlinearequation.h"

using namespace std;

template <class M, class B, class G>
class grad_shafranov : public non_linear_equation<M, B, G> {

public:
    grad_shafranov(M* mesh_ptr, double toll, double u_0_value, double p_0_value, double (*fptr)(double, double) = null_f)
        : non_linear_equation<M, B, G>(mesh_ptr,
              toll, fptr)
        , u_0(u_0_value)
        , p_0(p_0_value)
    {
    }

    ~grad_shafranov(void) {};

private:
    double u_0;
    double p_0;

    double eval_Jacobian(const int& k, const int& kk, const int& i, const int& ii, const int& j, const int& jj, const int& ele_pos);
    double eval_F0(const int& k, const int& kk, const int& i, const int& ii, const int& ele_pos);
    double eval_neumann(const double fun_value, const int& k, const int& i, const int& ele_pos);
};

template <class M, class B, class G>
double grad_shafranov<M, B, G>::eval_Jacobian(const int& k, const int& kk, const int& i, const int& ii, const int& j, const int& jj, const int& ele_pos)
{

    return (this->der_phi(k, j) * this->phi(kk, jj) * ((pow(this->J(0), 2) + pow(this->J(1), 2)) * this->der_phi(k, i) * this->phi(kk, ii) + (this->J(0) * this->J(2) + this->J(1) * this->J(3)) * this->phi(k, i) * this->der_phi(kk, ii)) + this->phi(k, j) * this->der_phi(kk, jj) * ((pow(this->J(2), 2) + pow(this->J(3), 2)) * this->phi(k, i) * this->der_phi(kk, ii) + (this->J(0) * this->J(2) + this->J(1) * this->J(3)) * this->der_phi(k, i) * this->phi(kk, ii))) * this->y_physical(k, kk, ele_pos)
        + 2 * this->phi(k, i) * this->phi(kk, ii) * (this->der_phi(k, j) * this->phi(kk, jj) * this->J(1) + this->phi(k, j) * this->der_phi(kk, jj) * this->J(3)) + 2 * (p_0 / pow(u_0, 2)) * pow(this->y_physical(k, kk, ele_pos), 3) * this->phi(k, j) * this->phi(kk, jj) * this->phi(k, i) * this->phi(kk, ii) * this->sol();
}

template <class M, class B, class G>
double grad_shafranov<M, B, G>::eval_F0(const int& k, const int& kk, const int& i, const int& ii, const int& ele_pos)
{

    return (this->sol_x() * ((pow(this->J(0), 2) + pow(this->J(1), 2)) * this->der_phi(k, i) * this->phi(kk, ii) + (this->J(0) * this->J(2) + this->J(1) * this->J(3)) * this->phi(k, i) * this->der_phi(kk, ii)) + this->sol_y() * ((pow(this->J(2), 2) + pow(this->J(3), 2)) * this->phi(k, i) * this->der_phi(kk, ii) + (this->J(0) * this->J(2) + this->J(1) * this->J(3)) * this->der_phi(k, i) * this->phi(kk, ii))) * this->y_physical(k, kk, ele_pos)
        + 2 * this->phi(k, i) * this->phi(kk, ii) * (this->sol_x() * this->J(1) + this->sol_y() * this->J(3)) - p_0 * pow(this->y_physical(k, kk, ele_pos), 3) * (1 - pow(this->sol() / u_0, 2)) * this->phi(k, i) * this->phi(kk, ii);
}

template <class M, class B, class G>
double grad_shafranov<M, B, G>::eval_neumann(const double fun_value, const int& k, const int& i, const int& ele_pos)
{
    return 0;
}

#endif /* GRADSHAFRANOV_H_ */
