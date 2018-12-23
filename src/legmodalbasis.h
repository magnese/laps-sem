/*
 * legmodalbasis.h
 *
 *  Created on: Oct 15, 2012
 *      Author: Marco Agnese
 */

#ifndef LEGMODALBASIS_H_
#define LEGMODALBASIS_H_

#include "basisfunction.h"
#include "legendre.h"
#include "utilities.h"

class leg_modal_basis : public basis_function {

public:
    leg_modal_basis(const int& j = 1)
        : leg_pol(j)
    {
        basis_coeff.resize(2);
        der_basis_coeff.resize(2);
        basis_coeff[0].resize(2);
        basis_coeff[0][0] = -0.5;
        basis_coeff[0][1] = 0.5;
        der_basis_coeff[0].resize(1);
        der_basis_coeff[0][0] = -0.5;
        basis_coeff[1].resize(2);
        basis_coeff[1][0] = 0.5;
        basis_coeff[1][1] = 0.5;
        der_basis_coeff[1].resize(1);
        der_basis_coeff[1][0] = 0.5;
        set_order(j);
    }

private:
    legendre leg_pol;

    virtual void set_basis_order(const int& j);
};

#endif /* LEGMODALBASIS_H_ */
