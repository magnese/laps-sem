/*
 * basisfunction.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: Marco Agnese
 */

#include "basisfunction.h"

void basis_function::set_der_order(void)
{
    unsigned db_size(der_basis_coeff.size());
    unsigned b_size(basis_coeff.size());
    if (b_size > db_size) {
        der_basis_coeff.resize(b_size);
        for (unsigned i = db_size; i != (b_size); ++i)
            poly_der(basis_coeff[i], der_basis_coeff[i]);
    }
}
