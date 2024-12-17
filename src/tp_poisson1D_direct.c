/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "lapacke.h"
#include <cblas.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double *AB;
  double *Y;  
  clock_t start, end;
  double time_used = 0.0;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  // Save a copy of RHS before it gets modified by the solver
  double* RHS_copy = (double*) malloc(sizeof(double)*la);
  memcpy(RHS_copy, RHS, la*sizeof(double));
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);
  Y = (double *) malloc(sizeof(double)*la);  // Allocation de Y

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  // Multiplication matrice-vecteur avec dgbmv
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB, lab, EX_SOL, 1, 0.0, Y, 1);
  write_vec(Y, &la, "Mvec.dat");

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if (IMPLEM == TRF) {
    clock_t start = clock();
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    clock_t end = clock();
    double time_fact = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for factorization = %lf\n", time_fact);
    
    if (info == 0) {
      start = clock();
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      end = clock();
      double time_solve = (double)(end - start) / CLOCKS_PER_SEC;
      printf("Time for solve = %lf\n", time_solve);
      printf("Total time = %lf\n", time_fact + time_solve);
    }
  }

  /* LU for tridiagonal matrix */
  if (IMPLEM == TRI) {
    clock_t start = clock();
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    clock_t end = clock();
    double time_fact = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time for factorization = %lf\n", time_fact);
    
    if (info == 0) {
      start = clock();
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      end = clock();
      double time_solve = (double)(end - start) / CLOCKS_PER_SEC;
      printf("Time for solve = %lf\n", time_solve);
      printf("Total time = %lf\n", time_fact + time_solve);
    }
  }

  /* Direct solve */
  if (IMPLEM == SV) {
    clock_t start = clock();
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    clock_t end = clock();
    double time_total = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total time = %lf\n", time_total);
  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Validation avec produit matrice-vecteur connu */
  double* b_exact = (double*) malloc(sizeof(double)*la);
  double* x_copy = (double*) malloc(sizeof(double)*la);
  
  // Sauvegarder une copie de la matrice AB avant qu'elle ne soit modifiée par dgbsv
  double* AB_copy = (double*) malloc(sizeof(double)*lab*la);
  memcpy(AB_copy, AB, sizeof(double)*lab*la);
  
  // Calculer b = A*x_exact
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB_copy, lab, EX_SOL, 1, 0.0, b_exact, 1);
  write_vec(b_exact, &la, "b_exact.dat");
  
  // Copier b_exact car il sera modifié par dgbsv
  memcpy(x_copy, b_exact, sizeof(double)*la);
  
  // Résoudre le système avec b_exact
  dgbsv_(&la, &kl, &ku, &NRHS, AB_copy, &lab, ipiv, x_copy, &la, &info);
  
  // Comparer la solution obtenue avec EX_SOL
  double validation_error = relative_forward_error(x_copy, EX_SOL, &la);
  printf("\nValidation avec solution exacte: ||x - x_exact|| / ||x_exact|| = %e\n", validation_error);
  
  free(b_exact);
  free(x_copy);
  free(AB_copy);

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);  
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(Y);  
  printf("\n\n--------- End -----------\n");
}
