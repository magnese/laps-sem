/*
 * utilities.cpp
 *
 *  Created on: Feb 11, 2013
 *      Author: Marco Agnese
 */

#include "utilities.h"

vector<double> poly_reverse_coeff(const vector<double>& coeff)
{
    vector<double> reverse_coeff(coeff);
    reverse(reverse_coeff.begin(), reverse_coeff.end());
    return reverse_coeff;
}

double poly_val(const vector<double>& coeff, const double& x)
{
    double value(0);
    for (vector<double>::const_reverse_iterator j = coeff.rbegin(); j != coeff.rend(); ++j)
        value += *j * pow(x, static_cast<int>(j - coeff.rbegin()));
    return value;
}

void poly_der(const vector<double>& coeff, vector<double>& der_coeff)
{
    der_coeff = coeff;
    for (vector<double>::reverse_iterator j = der_coeff.rbegin(); j != der_coeff.rend(); ++j)
        (*j) *= static_cast<int>(j - der_coeff.rbegin());
    der_coeff.resize(coeff.size() - 1);
}

int poly_deg(const vector<double>& coeff)
{
    int degree(coeff.size() - 1);
    vector<double>::const_iterator j = coeff.begin();
    while ((*j == 0) && (j != coeff.end()))
        ++j;
    degree -= static_cast<int>(j - coeff.begin());
    return degree;
}

int min_positive(const int& x, const int& y, const int& z)
{
    if (x < 0) {
        if (y < 0) {
            if (z < 0)
                return -1;
            else
                return z;
        } else {
            if (z < 0)
                return y;
            else
                return min(y, z);
        }
    } else {
        if (y < 0) {
            if (z < 0)
                return x;
            else
                return min(x, z);
        } else {
            if (z < 0)
                return min(x, y);
            else
                return min(min(x, y), min(x, z));
        }
    }
}
