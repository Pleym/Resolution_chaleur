#ifndef _LIB_POISSON1D_H_
#define _LIB_POISSON1D_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <cblas.h>

// Matrix and vector operations
void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, int* kv);
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int* la, int* kv);
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);
void set_grid_points_1D(double* x, int* la);
double relative_forward_error(double* x, double* y, int* la);

// Richardson method
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, 
                     int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite);
double richardson_alpha_opt(int *la);
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                  int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite);

// Jacobi and Gauss-Seidel methods
void jacobi_iteration(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                     int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite);
void gauss_seidel_iteration(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                           int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite);
void plot_iterative_convergence(double* error_vec, int nb_iter, char* method_name);

// Eigenvalues
void eig_poisson1D(double* eigval, int *la);
double eigmax_poisson1D(int *la);
double eigmin_poisson1D(int *la);

// Matrix extraction functions
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv);
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv);

// Writers
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename);
void write_vec(double* vec, int* la, char* filename);
void write_xy(double* vec, double* x, int* la, char* filename);

#endif
