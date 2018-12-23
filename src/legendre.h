/*
 * legendre.h
 *
 *  Created on: Oct 15, 2012
 *      Author: Marco Agnese
 */

#ifndef LEGENDRE_H_
#define LEGENDRE_H_

#include "utilities.h"

using namespace std;

class legendre {

public:
    legendre(const int& order = 1)
        : leg_coeff(2)
        , leg_der_coeff(2)
        , leg_coeff_reverse(0)
        , leg_der_coeff_reverse(0)
    {
        leg_coeff[0].resize(1);
        leg_coeff[0][0] = 1;
        leg_coeff[1].resize(2);
        leg_coeff[1][0] = 1;
        leg_coeff[1][1] = 0;
        leg_der_coeff[0].resize(1);
        leg_der_coeff[0][0] = 0;
        leg_der_coeff[1].resize(1);
        leg_der_coeff[1][0] = 1;
        set_order(order);
    }

    void set_order(const int& new_leg_order);
    inline int max_order(void) { return static_cast<int>(leg_coeff.size() - 1); }

    inline vector<double>& c(const int& j) { return leg_coeff[j]; }
    inline vector<double>& c_der(const int& j) { return leg_der_coeff[j]; }
    inline vector<double>& c_reverse(const int& j) { return leg_coeff_reverse[j]; }
    inline vector<double>& c_der_reverse(const int& j) { return leg_der_coeff_reverse[j]; }

    inline double eval(const double& x, const int& j) { return poly_val(c(j), x); }
    inline double eval_der(const double& x, const int& j) { return poly_val(c_der(j), x); }

    inline void disp(void)
    {
        cout << endl
             << "Legendre coefficients in decreasing order" << endl
             << endl;
        mat_disp<double>(leg_coeff);
        cout << endl;
    }
    inline void disp_der(void)
    {
        cout << endl
             << "Legendre derivative coefficients in decreasing order" << endl
             << endl;
        mat_disp<double>(leg_der_coeff);
        cout << endl;
    }
    inline void disp_reverse(void)
    {
        cout << endl
             << "Legendre coefficients in increasing order" << endl
             << endl;
        mat_disp<double>(leg_coeff_reverse);
        cout << endl;
    }
    inline void disp_der_reverse(void)
    {
        cout << endl
             << "Legendre derivative coefficients in increasing order" << endl
             << endl;
        mat_disp<double>(leg_der_coeff_reverse);
        cout << endl;
    }

private:
    vector<vector<double>> leg_coeff;
    vector<vector<double>> leg_der_coeff;
    vector<vector<double>> leg_coeff_reverse;
    vector<vector<double>> leg_der_coeff_reverse;

    void set_reverse_order(void);
};

#endif /* LEGENDRE_H_ */
