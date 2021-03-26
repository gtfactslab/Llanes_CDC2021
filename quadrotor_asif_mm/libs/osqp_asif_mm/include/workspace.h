#ifndef WORKSPACE_H
#define WORKSPACE_H

/*
 * This file was autogenerated by OSQP-Matlab on March 05, 2021 at 19:56:47.
 * 
 * This file contains the prototypes for all the workspace variables needed
 * by OSQP. The actual data is contained inside workspace.c.
 */

#include "types.h"
#include "qdldl_interface.h"

// Data structure prototypes
extern csc Pdata;
extern csc Adata;
extern c_float qdata[4];
extern c_float ldata[132];
extern c_float udata[132];
extern OSQPData data;

// Settings structure prototype
extern OSQPSettings settings;

// Scaling structure prototypes
extern c_float Dscaling[4];
extern c_float Dinvscaling[4];
extern c_float Escaling[132];
extern c_float Einvscaling[132];
extern OSQPScaling scaling;

// Prototypes for linsys_solver structure
extern csc linsys_solver_L;
extern c_float linsys_solver_Dinv[136];
extern c_int linsys_solver_P[136];
extern c_float linsys_solver_bp[136];
extern c_float linsys_solver_sol[136];
extern c_float linsys_solver_rho_inv_vec[132];
extern c_int linsys_solver_Pdiag_idx[4];
extern csc linsys_solver_KKT;
extern c_int linsys_solver_PtoKKT[4];
extern c_int linsys_solver_AtoKKT[520];
extern c_int linsys_solver_rhotoKKT[132];
extern QDLDL_float linsys_solver_D[136];
extern QDLDL_int linsys_solver_etree[136];
extern QDLDL_int linsys_solver_Lnz[136];
extern QDLDL_int   linsys_solver_iwork[408];
extern QDLDL_bool  linsys_solver_bwork[136];
extern QDLDL_float linsys_solver_fwork[136];
extern qdldl_solver linsys_solver;

// Prototypes for solution
extern c_float xsolution[4];
extern c_float ysolution[132];

extern OSQPSolution solution;

// Prototype for info structure
extern OSQPInfo info;

// Prototypes for the workspace
extern c_float work_rho_vec[132];
extern c_float work_rho_inv_vec[132];
extern c_int work_constr_type[132];
extern c_float work_x[4];
extern c_float work_y[132];
extern c_float work_z[132];
extern c_float work_xz_tilde[136];
extern c_float work_x_prev[4];
extern c_float work_z_prev[132];
extern c_float work_Ax[132];
extern c_float work_Px[4];
extern c_float work_Aty[4];
extern c_float work_delta_y[132];
extern c_float work_Atdelta_y[4];
extern c_float work_delta_x[4];
extern c_float work_Pdelta_x[4];
extern c_float work_Adelta_x[132];
extern c_float work_D_temp[4];
extern c_float work_D_temp_A[4];
extern c_float work_E_temp[132];

extern OSQPWorkspace workspace;

#endif // ifndef WORKSPACE_H