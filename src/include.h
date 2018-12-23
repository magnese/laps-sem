/*
 * include.h
 *
 *  Created on: May 13, 2013
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#ifndef INCLUDE_H_
#define INCLUDE_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "math.h"
#include "metis.h"
#include "mpi.h"
#include "parmetis.h"
#include "petscvec.h"

typedef enum { DIRICHLET,
    NEUMANN,
    PERIODICITY,
    WALL,
    NOCONDITION } flag_bc;
typedef enum { WRITE_GHOSTS_ON,
    WRITE_GHOSTS_OFF } flag_ghost;

typedef enum { READ_MONO,
    READ_PARALLEL_PRE,
    READ_PARALLEL_PAR } flag_read;
typedef enum { WRITE_OFF,
    WRITE_PREPARTITIONING,
    WRITE_PARTITIONING,
    WRITE_REFINEMENT,
    WRITE_ALL } flag_write;

typedef enum { CLOSED,
    OPENED } loop;

typedef enum { OUTPUT_SEQ,
    OUTPUT_MPI,
    OUTPUT_NODAL } flag_output_type;

inline double null_fun(double x, double y, double t) { return 0; }
inline double null_f(double x, double y) { return 0; }

#endif /* INCLUDE_H_ */
