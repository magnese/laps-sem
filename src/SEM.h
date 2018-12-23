/*
 * SEM.h
 *
 *  Created on: Feb 22, 2013
 *      Author: Marco Agnese
 */

#ifndef SEM_H_
#define SEM_H_

#include "modesmapping.h"
#include "plasmastate.h"
#include "utilities.h"

using namespace std;

static PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void* dummy)
{

    PetscPrintf(PETSC_COMM_WORLD, "iteration %D solution vector:\n", n);
    PetscPrintf(PETSC_COMM_WORLD, "iteration %D KSP Residual norm %14.12e \n", n, rnorm);
    return 0;
}

template <class M, class E, class B>
class SEM {

public:
    SEM(M* mesh_ptr, E* equation_ptr, static_condesation_flag flag = STATIC_CONDENSATION_OFF)
        : mesh(mesh_ptr)
        , equation(equation_ptr)
        , p_state(mesh, equation->get_basis())
        , map(mesh, equation->get_basis(), &p_state)
        , static_condensation_status(false)
        , values(0)
        , num_max_iteration(50)
    {
        if (flag == STATIC_CONDENSATION_ON)
            static_condensation_status = true;
        mesh->update_ghosts_deg();
        map.create_mapping(static_condensation_status);
        equation->set_ptr_mode_mapping(&map);
        p_state.set_ptr_mode_mapping(&map);
        create_petsc_variables();
        create_scatter();
        if (equation->get_equation_type() == LINEAR_EQUATION) {
            if (static_condensation_status) {
                if (mesh->get_rank() == 0)
                    cout << "START: Assembling linear system with static condensation" << endl;
                assemble();
                assemble_dirichlet_condition_LHS();
                assemble_dirichlet_condition_RHS();
                if (mesh->get_rank() == 0)
                    cout << "END: Assembling linear system with static condensation" << endl
                         << endl;
                if (mesh->get_rank() == 0)
                    cout << "START: Solving linear system with static condensation" << endl;
                solve();
                if (mesh->get_rank() == 0)
                    cout << "END: Solving linear system with static condensation" << endl
                         << endl;
            } else {
                if (mesh->get_rank() == 0)
                    cout << "START: Assembling linear system without static condensation" << endl;
                assemble_no_static_condensation();
                assemble_dirichlet_condition_LHS();
                assemble_dirichlet_condition_RHS();
                if (mesh->get_rank() == 0)
                    cout << "END: Assembling linear system without static condensation" << endl
                         << endl;
                if (mesh->get_rank() == 0)
                    cout << "START: Solving linear system without static condensation" << endl;
                solve_no_static_condensation();
                if (mesh->get_rank() == 0)
                    cout << "END: Solving linear system without static condensation" << endl
                         << endl;
            }
            scatter_solution();
        } else {
            if (static_condensation_status)
                nl_solve();
            else
                nl_solve_no_static_condensation();
        }
        if (mesh->get_rank() == 0)
            cout << "START: Storing solution in plasma state" << endl;
        p_state.update((mesh->get_max_deg() - 1));
        if (mesh->get_rank() == 0)
            cout << "END: Storing solution in plasma state" << endl
                 << endl;
    }

    inline plasma_state<M, B>* get_ps(void) { return &p_state; }

    inline void disp_nodes_values(void) { p_state.disp_nodes_values(); }
    inline void disp_Abb(void) { MatView(Abb, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Abi(void) { MatView(Abi, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Aii(void) { MatView(Aii, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Fb(void) { VecView(Fb, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Fi(void) { VecView(Fi, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Xb(void) { VecView(Xb, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Xi(void) { VecView(Xi, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Xb_sol(void) { VecView(Xb_sol, PETSC_VIEWER_STDOUT_WORLD); }
    inline void disp_Xi_sol(void) { VecView(Xi_sol, PETSC_VIEWER_STDOUT_WORLD); }

    void free_petsc_variables(void);

    ~SEM() {}

private:
    M* mesh;
    E* equation;
    plasma_state<M, B> p_state;
    modes_mapping<M, B, plasma_state<M, B>> map;

    bool static_condensation_status;

    Mat Abb;
    Mat Abi;
    Mat Aii;

    Vec Fb;
    Vec Fi;

    Vec Xb;
    Vec Xi;

    Vec Xb_sol;
    Vec Xi_sol;

    Vec X_seq;
    VecScatter scatter_X;
    double* values;

    int num_max_iteration;

    void create_petsc_variables(void);
    void initialize_petsc_variables(void);
    void create_scatter(void);
    void scatter_solution(void);

    void assemble(void);
    void solve(void);
    void nl_solve(void);

    void assemble_no_static_condensation(void);
    void solve_no_static_condensation(void);
    void nl_solve_no_static_condensation(void);

    void assemble_dirichlet_condition_LHS(void);
    void assemble_dirichlet_condition_RHS(void);
};

template <class M, class E, class B>
void SEM<M, E, B>::create_petsc_variables(void)
{

    if (static_condensation_status) {
        MatCreate(MPI_COMM_WORLD, &Abb);
        MatSetSizes(Abb, PETSC_DECIDE, PETSC_DECIDE, map.dim_boundary(), map.dim_boundary());
        MatSetType(Abb, MATMPIAIJ);
        MatSetUp(Abb);
        //MatMPIAIJSetPreallocation(Abb,nz,NULL,onz,NULL); //TODO: preallocation!!

        MatCreate(MPI_COMM_WORLD, &Abi);
        MatSetSizes(Abi, PETSC_DECIDE, PETSC_DECIDE, map.dim_boundary(), map.dim_interiors());
        MatSetType(Abi, MATMPIAIJ);
        MatSetUp(Abi);
        //MatMPIAIJSetPreallocation(Abb,nz,NULL,onz,NULL); //TODO: preallocation!!

        MatCreate(MPI_COMM_WORLD, &Aii);
        MatSetSizes(Aii, PETSC_DECIDE, PETSC_DECIDE, map.dim_interiors(), map.dim_interiors());
        MatSetType(Aii, MATMPIAIJ);
        MatMPIAIJSetPreallocation(Aii, (mesh->get_max_deg() - 1) * (mesh->get_max_deg() - 1), NULL, 0, NULL);

        VecCreate(MPI_COMM_WORLD, &Fb);
        VecSetSizes(Fb, PETSC_DECIDE, map.dim_boundary());
        VecSetType(Fb, VECMPI);

        VecCreate(MPI_COMM_WORLD, &Fi);
        VecSetSizes(Fi, PETSC_DECIDE, map.dim_interiors());
        VecSetType(Fi, VECMPI);

        VecDuplicate(Fb, &Xb);
        VecDuplicate(Fi, &Xi);

        if (equation->get_equation_type() == NON_LINEAR_EQUATION) {
            VecDuplicate(Xb, &Xb_sol);
            VecDuplicate(Xi, &Xi_sol);
        }
    } else {
        MatCreate(MPI_COMM_WORLD, &Abb);
        MatSetSizes(Abb, PETSC_DECIDE, PETSC_DECIDE, map.dim_boundary() + map.dim_interiors(), map.dim_boundary() + map.dim_interiors());
        MatSetType(Abb, MATMPIAIJ);
        MatSetUp(Abb);
        //MatMPIAIJSetPreallocation(Abb,nz,NULL,onz,NULL); //TODO: preallocation!!

        VecCreate(MPI_COMM_WORLD, &Fb);
        VecSetType(Fb, VECMPI);
        VecSetSizes(Fb, PETSC_DECIDE, map.dim_boundary() + map.dim_interiors());

        VecDuplicate(Fb, &Xb);

        if (equation->get_equation_type() == NON_LINEAR_EQUATION)
            VecDuplicate(Xb, &Xb_sol);
    }
}

template <class M, class E, class B>
void SEM<M, E, B>::create_scatter(void)
{

    if (equation->get_equation_type() == NON_LINEAR_EQUATION)
        VecScatterCreateToAll(Xb_sol, &scatter_X, &X_seq);
    else
        VecScatterCreateToAll(Xb, &scatter_X, &X_seq);
}

template <class M, class E, class B>
void SEM<M, E, B>::scatter_solution(void)
{

    if (equation->get_equation_type() == NON_LINEAR_EQUATION) {
        VecScatterBegin(scatter_X, Xb_sol, X_seq, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter_X, Xb_sol, X_seq, INSERT_VALUES, SCATTER_FORWARD);
    } else {
        VecScatterBegin(scatter_X, Xb, X_seq, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter_X, Xb, X_seq, INSERT_VALUES, SCATTER_FORWARD);
    }

    VecGetArray(X_seq, &values);
    map.set_ptr_sol_coeff(&values);
}

template <class M, class E, class B>
void SEM<M, E, B>::initialize_petsc_variables(void)
{

    MatZeroEntries(Abb);
    MatAssemblyBegin(Abb, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Abb, MAT_FINAL_ASSEMBLY);

    VecZeroEntries(Fb);
    VecAssemblyBegin(Fb);
    VecAssemblyEnd(Fb);

    if (static_condensation_status) {

        MatZeroEntries(Abi);
        MatAssemblyBegin(Abi, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Abi, MAT_FINAL_ASSEMBLY);

        MatZeroEntries(Aii);
        MatAssemblyBegin(Aii, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Aii, MAT_FINAL_ASSEMBLY);

        VecZeroEntries(Fi);
        VecAssemblyBegin(Fi);
        VecAssemblyEnd(Fi);
    }
}

template <class M, class E, class B>
void SEM<M, E, B>::free_petsc_variables(void)
{

    MatDestroy(&Abb);
    VecDestroy(&Fb);
    VecDestroy(&Xb);

    if (equation->get_equation_type() == NON_LINEAR_EQUATION)
        VecDestroy(&Xb_sol);

    if (static_condensation_status) {
        MatDestroy(&Abi);
        MatDestroy(&Aii);
        VecDestroy(&Fi);
        VecDestroy(&Xi);

        if (equation->get_equation_type() == NON_LINEAR_EQUATION)
            VecDestroy(&Xi_sol);
    }

    VecScatterDestroy(&scatter_X);
    VecDestroy(&X_seq);
}

template <class M, class E, class B>
void SEM<M, E, B>::assemble(void)
{

    // TODO: optimization when matrices are symmetric

    for (int e = 0; e != mesh->num_elements(); ++e) {
        for (int i = 0; i != (mesh->get_deg(e) + 1); ++i) {
            for (int ii = 0; ii != (mesh->get_deg(e) + 1); ++ii) {
                if (map.is_i(i, ii))
                    VecSetValue(Fi, map.inner(i, ii, e), equation->RHS(i, ii, e), ADD_VALUES);
                else {
                    if (map.is_v(i, ii)) {
                        if ((mesh->n(map.vtx_pos(i, ii), e) > -1) && (mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) > -1))
                            VecSetValue(Fb, map.vtx(i, ii, e), equation->RHS(i, ii, e), ADD_VALUES);
                        else {
                            if ((mesh->n(map.vtx_pos(i, ii), e) < 0) && (mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) > -1)) {
                                if (mesh->get_bc_type(-mesh->n(map.vtx_pos(i, ii), e) - 1) != DIRICHLET)
                                    VecSetValue(Fb, map.vtx(i, ii, e), equation->neumann(i, map.vtx_pos(i, ii), e) + equation->RHS(i, ii, e), ADD_VALUES);
                            } else {
                                if ((mesh->n(map.vtx_pos(i, ii), e) > -1) && (mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) < 0)) {
                                    if (mesh->get_bc_type(-mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) - 1) != DIRICHLET)
                                        VecSetValue(Fb, map.vtx(i, ii, e), equation->neumann(ii, (map.vtx_pos(i, ii) + 3) % 4, e) + equation->RHS(i, ii, e), ADD_VALUES);
                                } else {
                                    if ((mesh->get_bc_type(-mesh->n(map.vtx_pos(i, ii), e) - 1) != DIRICHLET) && (mesh->get_bc_type(-mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) - 1) != DIRICHLET))
                                        VecSetValue(Fb, map.vtx(i, ii, e), equation->neumann(i, map.vtx_pos(i, ii), e) + equation->RHS(i, ii, e), ADD_VALUES);
                                }
                            }
                        }
                    } else {
                        if (mesh->n(map.edge_pos(i, ii), e) < 0) {
                            if (mesh->get_bc_type(-mesh->n(map.edge_pos(i, ii), e) - 1) != DIRICHLET)
                                VecSetValue(Fb, map.edge(i, ii, e), equation->neumann(max(i, ii), map.edge_pos(i, ii), e) + equation->RHS(i, ii, e), ADD_VALUES);
                        } else
                            VecSetValue(Fb, map.edge(i, ii, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->RHS(i, ii, e), ADD_VALUES);
                    }
                }
                for (int j = 0; j != (mesh->get_deg(e) + 1); ++j) {
                    for (int jj = 0; jj != (mesh->get_deg(e) + 1); ++jj) {
                        if (map.is_i(j, jj)) {
                            if (map.is_i(i, ii))
                                MatSetValue(Aii, map.inner(i, ii, e), map.inner(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                            else {
                                if (map.is_e(i, ii)) {
                                    if (map.edge(i, ii, e) > (map.dim_ndof_nodes() + map.dim_ndof_edges() - 1))
                                        MatSetValue(Abi, map.edge(i, ii, e), map.inner(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                } else {
                                    if (map.vtx(i, ii, e) > (map.dim_ndof_nodes() - 1))
                                        MatSetValue(Abi, map.vtx(i, ii, e), map.inner(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                }
                            }
                        } else {
                            if (map.is_e(j, jj)) {
                                if (map.is_e(i, ii)) {
                                    if (map.edge(i, ii, e) > (map.dim_ndof_nodes() + map.dim_ndof_edges() - 1))
                                        MatSetValue(Abb, map.edge(i, ii, e), map.edge(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * pow(mesh->orientation(map.edge_pos(j, jj), e), max(j, jj)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                } else {
                                    if (map.is_v(i, ii)) {
                                        if (map.vtx(i, ii, e) > (map.dim_ndof_nodes() - 1))
                                            MatSetValue(Abb, map.vtx(i, ii, e), map.edge(j, jj, e), pow(mesh->orientation(map.edge_pos(j, jj), e), max(j, jj)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                    }
                                }
                            } else {
                                if (map.is_e(i, ii)) {
                                    if (map.edge(i, ii, e) > (map.dim_ndof_nodes() + map.dim_ndof_edges() - 1))
                                        MatSetValue(Abb, map.edge(i, ii, e), map.vtx(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                } else {
                                    if (map.is_v(i, ii)) {
                                        if (map.vtx(i, ii, e) > (map.dim_ndof_nodes() - 1))
                                            MatSetValue(Abb, map.vtx(i, ii, e), map.vtx(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    MatAssemblyBegin(Abb, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Abb, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(Abi, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Abi, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(Aii, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Aii, MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(Fb);
    VecAssemblyEnd(Fb);

    VecAssemblyBegin(Fi);
    VecAssemblyEnd(Fi);
}

template <class M, class E, class B>
void SEM<M, E, B>::solve(void)
{

    //TODO
}

template <class M, class E, class B>
void SEM<M, E, B>::assemble_no_static_condensation(void)
{

    // TODO: optimization when matrices are symmetric

    for (int e = 0; e != mesh->num_elements(); ++e) {
        for (int i = 0; i != (mesh->get_deg(e) + 1); ++i) {
            for (int ii = 0; ii != (mesh->get_deg(e) + 1); ++ii) {
                if (map.is_i(i, ii))
                    VecSetValue(Fb, map.inner(i, ii, e), equation->RHS(i, ii, e), ADD_VALUES);
                else {
                    if (map.is_v(i, ii)) {
                        if ((mesh->n(map.vtx_pos(i, ii), e) > -1) && (mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) > -1))
                            VecSetValue(Fb, map.vtx(i, ii, e), equation->RHS(i, ii, e), ADD_VALUES);
                        else {
                            if ((mesh->n(map.vtx_pos(i, ii), e) < 0) && (mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) > -1)) {
                                if (mesh->get_bc_type(-mesh->n(map.vtx_pos(i, ii), e) - 1) != DIRICHLET)
                                    VecSetValue(Fb, map.vtx(i, ii, e), equation->neumann(i, map.vtx_pos(i, ii), e) + equation->RHS(i, ii, e), ADD_VALUES);
                            } else {
                                if ((mesh->n(map.vtx_pos(i, ii), e) > -1) && (mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) < 0)) {
                                    if (mesh->get_bc_type(-mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) - 1) != DIRICHLET)
                                        VecSetValue(Fb, map.vtx(i, ii, e), equation->neumann(ii, (map.vtx_pos(i, ii) + 3) % 4, e) + equation->RHS(i, ii, e), ADD_VALUES);
                                } else {
                                    if ((mesh->get_bc_type(-mesh->n(map.vtx_pos(i, ii), e) - 1) != DIRICHLET) && (mesh->get_bc_type(-mesh->n((map.vtx_pos(i, ii) + 3) % 4, e) - 1) != DIRICHLET))
                                        VecSetValue(Fb, map.vtx(i, ii, e), equation->neumann(i, map.vtx_pos(i, ii), e) + equation->RHS(i, ii, e), ADD_VALUES);
                                }
                            }
                        }
                    } else {
                        if (mesh->n(map.edge_pos(i, ii), e) < 0) {
                            if (mesh->get_bc_type(-mesh->n(map.edge_pos(i, ii), e) - 1) != DIRICHLET)
                                VecSetValue(Fb, map.edge(i, ii, e), equation->neumann(max(i, ii), map.edge_pos(i, ii), e) + equation->RHS(i, ii, e), ADD_VALUES);
                        } else
                            VecSetValue(Fb, map.edge(i, ii, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->RHS(i, ii, e), ADD_VALUES);
                    }
                }
                for (int j = 0; j != (mesh->get_deg(e) + 1); ++j) {
                    for (int jj = 0; jj != (mesh->get_deg(e) + 1); ++jj) {
                        if (map.is_i(j, jj)) {
                            if (map.is_i(i, ii))
                                MatSetValue(Abb, map.inner(i, ii, e), map.inner(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                            else {
                                if (map.is_e(i, ii)) {
                                    if (map.edge(i, ii, e) > (map.dim_ndof_nodes() + map.dim_ndof_edges() - 1))
                                        MatSetValue(Abb, map.edge(i, ii, e), map.inner(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                } else {
                                    if (map.vtx(i, ii, e) > (map.dim_ndof_nodes() - 1))
                                        MatSetValue(Abb, map.vtx(i, ii, e), map.inner(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                }
                            }
                        } else {
                            if (map.is_e(j, jj)) {
                                if (map.is_i(i, ii))
                                    MatSetValue(Abb, map.inner(i, ii, e), map.edge(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                else {
                                    if (map.is_e(i, ii)) {
                                        if (map.edge(i, ii, e) > (map.dim_ndof_nodes() + map.dim_ndof_edges() - 1))
                                            MatSetValue(Abb, map.edge(i, ii, e), map.edge(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * pow(mesh->orientation(map.edge_pos(j, jj), e), max(j, jj)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                    } else {
                                        if (map.vtx(i, ii, e) > (map.dim_ndof_nodes() - 1))
                                            MatSetValue(Abb, map.vtx(i, ii, e), map.edge(j, jj, e), pow(mesh->orientation(map.edge_pos(j, jj), e), max(j, jj)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                    }
                                }
                            } else {
                                if (map.is_i(i, ii))
                                    MatSetValue(Abb, map.inner(i, ii, e), map.vtx(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                else {
                                    if (map.is_e(i, ii)) {
                                        if (map.edge(i, ii, e) > (map.dim_ndof_nodes() + map.dim_ndof_edges() - 1))
                                            MatSetValue(Abb, map.edge(i, ii, e), map.vtx(j, jj, e), pow(mesh->orientation(map.edge_pos(i, ii), e), max(i, ii)) * equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                    } else {
                                        if (map.vtx(i, ii, e) > (map.dim_ndof_nodes() - 1))
                                            MatSetValue(Abb, map.vtx(i, ii, e), map.vtx(j, jj, e), equation->LHS(i, ii, j, jj, e), ADD_VALUES);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    MatAssemblyBegin(Abb, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Abb, MAT_FINAL_ASSEMBLY);

    VecAssemblyBegin(Fb);
    VecAssemblyEnd(Fb);
}

template <class M, class E, class B>
void SEM<M, E, B>::solve_no_static_condensation(void)
{

    KSP ksp;

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, Abb, Abb, DIFFERENT_NONZERO_PATTERN);
    KSPSetTolerances(ksp, 1.e-16, 1.e-16, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    //KSPMonitorSet(ksp,MyKSPMonitor,NULL,0);
    KSPSolve(ksp, Fb, Xb);
    KSPDestroy(&ksp);
}

template <class M, class E, class B>
void SEM<M, E, B>::nl_solve_no_static_condensation(void)
{

    double vec_norm(0);
    bool flag(true);

    if (mesh->get_rank() == 0)
        cout << "START: Solving non linear system without static condensation" << endl;

    if (!(equation->is_initial_guess_projected())) {
        assemble_no_static_condensation();
        assemble_dirichlet_condition_LHS();
        assemble_dirichlet_condition_RHS();
        solve_no_static_condensation();
        VecCopy(Xb, Xb_sol);
        equation->set_initial_guess_projected();
    } else {
        assemble_dirichlet_condition_RHS();
        VecCopy(Fb, Xb_sol);
    }
    scatter_solution();

    while (flag && equation->get_iter() < num_max_iteration) {
        initialize_petsc_variables();
        assemble_no_static_condensation();
        assemble_dirichlet_condition_LHS();
        VecNorm(Fb, NORM_2, &vec_norm);
        if (vec_norm < equation->get_toll())
            flag = false;
        else {
            solve_no_static_condensation();
            VecAYPX(Xb_sol, 1, Xb);
            scatter_solution();
        }
        equation->inc_iter();
        if (mesh->get_rank() == 0)
            cout << "Iteration number " << equation->get_iter() << ", L^2 norm of the residual " << vec_norm << endl;
    }

    if (equation->get_iter() == num_max_iteration) {
        if (mesh->get_rank() == 0)
            cout << endl
                 << "ERROR: Reached maximum number of iterations before convergence" << endl
                 << endl;
    }

    if (mesh->get_rank() == 0)
        cout << "END: Solving non linear system without static condensation" << endl
             << endl;
}

template <class M, class E, class B>
void SEM<M, E, B>::nl_solve(void)
{

    double vec_norm_1(0);
    double vec_norm_2(0);
    bool flag(true);

    if (mesh->get_rank() == 0)
        cout << "START: Solving non linear system with static condensation" << endl;

    if (!(equation->is_initial_guess_projected())) {
        assemble();
        assemble_dirichlet_condition_LHS();
        assemble_dirichlet_condition_RHS();
        solve();
        VecCopy(Xb, Xb_sol);
        VecCopy(Xi, Xi_sol);
        equation->set_initial_guess_projected();
    } else {
        assemble_dirichlet_condition_RHS();
        VecCopy(Fb, Xb_sol);
    }
    scatter_solution();

    while (flag && equation->get_iter() < num_max_iteration) {
        initialize_petsc_variables();
        assemble();
        assemble_dirichlet_condition_LHS();
        VecNorm(Fb, NORM_2, &vec_norm_1);
        VecNorm(Fi, NORM_2, &vec_norm_2);
        if ((vec_norm_1 + vec_norm_2) < equation->get_toll())
            flag = false;
        else {
            solve();
            VecAYPX(Xb_sol, 1, Xb);
            VecAYPX(Xi_sol, 1, Xi);
            scatter_solution();
            equation->inc_iter();
        }

        if (mesh->get_rank() == 0)
            cout << "Iteration number " << equation->get_iter() << ", L^2 norm of the residual " << vec_norm_1 + vec_norm_2 << endl;
    }

    if (equation->get_iter() == num_max_iteration) {
        if (mesh->get_rank() == 0)
            cout << endl
                 << "ERROR: Reached maximum number of iterations before convergence" << endl
                 << endl;
    }

    if (mesh->get_rank() == 0)
        cout << "END: Solving non linear system with static condensation" << endl
             << endl;
}

template <class M, class E, class B>
void SEM<M, E, B>::assemble_dirichlet_condition_LHS(void)
{

    //TODO: each process insert a part instead of only the first
    if (mesh->get_rank() == 0) {
        for (int j = 0; j != (map.dim_ndof_nodes() + map.dim_ndof_edges()); ++j)
            MatSetValue(Abb, j, j, 1, INSERT_VALUES);
    }
    MatAssemblyBegin(Abb, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Abb, MAT_FINAL_ASSEMBLY);
}

template <class M, class E, class B>
void SEM<M, E, B>::assemble_dirichlet_condition_RHS(void)
{

    //TODO: each process insert a part instead of only the first

    if (mesh->get_rank() == 0) {
        for (int l = 0; l != mesh->num_lines(); ++l) {
            if (mesh->get_bc_type(l) == DIRICHLET) {
                VecSetValue(Fb, map.get_nodes_mapping()[mesh->idx_node_line(l) - mesh->idx_first_node()], mesh->eval_bc(l, mesh->x_node(mesh->idx_node_line(l) - mesh->idx_first_node()), mesh->y_node(mesh->idx_node_line(l) - mesh->idx_first_node()), 0), INSERT_VALUES);
                VecSetValue(Fb, map.get_nodes_mapping()[mesh->idx_node_line(l, 1) - mesh->idx_first_node()], mesh->eval_bc(l, mesh->x_node(mesh->idx_node_line(l, 1) - mesh->idx_first_node()), mesh->y_node(mesh->idx_node_line(l, 1) - mesh->idx_first_node()), 0), INSERT_VALUES);
            }
        }
    }
    for (int e = 0; e != mesh->num_elements(); ++e) {
        for (int i = 0; i != (mesh->get_deg(e) + 1); ++i) {
            for (int ii = 0; ii != (mesh->get_deg(e) + 1); ++ii) {
                if (map.is_e(i, ii)) {
                    if (mesh->n(map.edge_pos(i, ii), e) < 0) {
                        if (mesh->get_bc_type(-mesh->n(map.edge_pos(i, ii), e) - 1) == DIRICHLET)
                            VecSetValue(Fb, map.edge(i, ii, e), equation->dirichlet(max(i, ii), map.edge_pos(i, ii), e), INSERT_VALUES);
                    }
                }
            }
        }
    }
    VecAssemblyBegin(Fb);
    VecAssemblyEnd(Fb);
}

#endif /* SEM_H_ */
