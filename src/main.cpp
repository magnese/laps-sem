/*
 * main.h

 *
 *  Created on: Dec 4, 2012
 *      Author: Marco Agnese
 */

#include "meshsem.h"

#include "gaussianintegration.h"
#include "legmodalbasis.h"

#include "bratu.h"
#include "gradshafranov.h"
#include "lineargradshafranov.h"
#include "lineargradshafranovmod.h"
#include "poisson.h"

#include "SEM.h"

#include "output.h"

#include <time.h>

using namespace std;

static char help[] = "SEM.\n\n";

inline double source_term(double x, double y) { return 0; }
//inline double source_term(double x,double y){return -2*x*x-2*y*y+4;}
//inline double source_term(double x,double y){return 2*pow(3.1415926,2)*sin(3.1415926*x)*sin(3.1415926*y);}

inline double dir_bc(double x, double y, double t) { return 2; }

inline double initial_guess(double x, double y) { return 0; }

int main(int argc, char** argv)
{

    /*initialization*/

    clock_t t(0);
    int my_rank;
    int my_size;
    int order(2);

    PetscInitialize(&argc, &argv, (char*)0, help);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &my_size);

    if (my_rank == 0)
        t = clock();

    /*mesh*/

    /*
  string folder("square_4k");
  mesh_sem my_mesh(my_rank,my_size,4,"../mesh/"+folder+"/square_4k.msh",READ_MONO,WRITE_OFF,order);
  my_mesh.set_border_function(dir_bc,5,260,DIRICHLET);
  */

    string folder("mirror_2k");
    mesh_sem my_mesh(my_rank, my_size, 4, "../mesh/" + folder + "/mirror_field_lines_frontal.msh", READ_MONO, WRITE_OFF, order);
    my_mesh.set_periodicity(93, 106, 198, 185);
    my_mesh.set_border_function(dir_bc, 107, 184, DIRICHLET);
    my_mesh.set_border_function(null_fun, 199, 276, NEUMANN);

    /*
  string folder("tokamak_1.7k");
  mesh_sem my_mesh(my_rank,my_size,4,"../mesh/"+folder+"/tokamak.msh",READ_MONO,WRITE_OFF,order);
  my_mesh.set_border_function(dir_bc,80,237,DIRICHLET);
  */

    /*problem*/

    //poisson<mesh_sem,leg_modal_basis,gaussian_integration> my_equation(&my_mesh,source_term);
    //SEM<mesh_sem,poisson<mesh_sem,leg_modal_basis,gaussian_integration>,leg_modal_basis> my_simulator(&my_mesh,&my_equation,STATIC_CONDENSATION_OFF);

    //linear_grad_shafranov<mesh_sem,leg_modal_basis,gaussian_integration> my_equation(&my_mesh,source_term);
    //SEM<mesh_sem,linear_grad_shafranov<mesh_sem,leg_modal_basis,gaussian_integration>,leg_modal_basis > my_simulator(&my_mesh,&my_equation,STATIC_CONDENSATION_OFF);

    //linear_grad_shafranov_mod<mesh_sem,leg_modal_basis,gaussian_integration> my_equation(&my_mesh,1,1,source_term);
    //SEM<mesh_sem,linear_grad_shafranov_mod<mesh_sem,leg_modal_basis,gaussian_integration>,leg_modal_basis > my_simulator(&my_mesh,&my_equation,STATIC_CONDENSATION_OFF);

    //bratu<mesh_sem,leg_modal_basis,gaussian_integration> my_equation(&my_mesh,1.e-7,1,initial_guess);
    //SEM<mesh_sem,bratu<mesh_sem,leg_modal_basis,gaussian_integration>,leg_modal_basis> my_simulator(&my_mesh,&my_equation,STATIC_CONDENSATION_OFF);

    grad_shafranov<mesh_sem, leg_modal_basis, gaussian_integration> my_equation(&my_mesh, 1.e-10, 1, 1, initial_guess);
    SEM<mesh_sem, grad_shafranov<mesh_sem, leg_modal_basis, gaussian_integration>, leg_modal_basis> my_simulator(&my_mesh, &my_equation, STATIC_CONDENSATION_OFF);

    /*output*/

    output<mesh_sem, plasma_state<mesh_sem, leg_modal_basis>> dump(&my_mesh, my_simulator.get_ps());
    dump.gmsh_dump("../mesh/" + folder + "/sol.msh");
    dump.dump_binary("../mesh/" + folder + "/values", true);
    //dump.dump("../mesh/"+folder+"/values",true);

    /*finalization*/

    my_simulator.free_petsc_variables();

    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank == 0) {
        t = clock() - t;
        cout << endl
             << "Time: " << static_cast<double>(t) / CLOCKS_PER_SEC << " s" << endl
             << endl;
    }

    PetscFinalize();

    return 0;
}
