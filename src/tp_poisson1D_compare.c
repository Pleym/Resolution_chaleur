#include "lib_poisson1D.h"
#include <string.h>

int main(int argc, char *argv[]) {
    int ierr;
    int jj;
    int nbpoints, la;
    int ku, kl, lab, kv;
    int *ipiv;
    int info;
    int NRHS;
    double T0, T1;
    double *RHS, *SOL, *EX_SOL, *X;
    double *AB;
    double *MB;
    double temp, relres;
    double tol;
    int maxit;
    double *resvec;
    int nbite;

    // Set problem parameters
    nbpoints = 102;  // Increased problem size
    la = nbpoints - 2;
    T0 = -5.0;
    T1 = 5.0;
    tol = 1e-10;    // Increased tolerance
    maxit = 50;

    printf("Taille du problème : la = %d\n", la);
    printf("Conditions aux limites : T0 = %f, T1 = %f\n", T0, T1);
    printf("Tolérance : %e, Nombre maximum d'itérations : %d\n", tol, maxit);

    // Allocate arrays
    RHS = (double *)calloc(la, sizeof(double));
    SOL = (double *)calloc(la, sizeof(double));
    EX_SOL = (double *)calloc(la, sizeof(double));
    X = (double *)calloc(la, sizeof(double));
    resvec = (double *)calloc(maxit, sizeof(double));

    // Set grid points
    set_grid_points_1D(X, &la);
    printf("\nPremiers points de la grille : %f %f %f\n", X[0], X[1], X[2]);

    // Set RHS and exact solution
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    printf("\nPremières valeurs du second membre : %f %f %f\n", RHS[0], RHS[1], RHS[2]);
    printf("Premières valeurs de la solution exacte : %f %f %f\n", EX_SOL[0], EX_SOL[1], EX_SOL[2]);

    // Define the band matrix
    kv = 1;
    ku = 1;
    kl = 1;
    lab = kv + kl + ku + 1;
    AB = (double *)calloc(lab * la, sizeof(double));
    MB = (double *)calloc(lab * la, sizeof(double));

    // Set the band matrix
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    printf("\nDimensions de la matrice bande : lab = %d, la = %d\n", lab, la);
    printf("Premiers éléments diagonaux : %f %f %f\n", 
           AB[lab*0 + ku], AB[lab*1 + ku], AB[lab*2 + ku]);

    // Test Richardson method
    printf("\nTest de la méthode de Richardson...\n");
    memset(SOL, 0, la * sizeof(double));  // Initialize solution to 0
    memset(resvec, 0, maxit * sizeof(double));
    double alpha = richardson_alpha_opt(&la);
    printf("Alpha optimal pour Richardson : %f\n", alpha);
    richardson_alpha(AB, RHS, SOL, &alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("Méthode de Richardson : nombre d'itérations = %d\n", nbite);
    printf("Méthode de Richardson : erreur avant relative = %e\n", relres);
    printf("Premières valeurs de la solution : %f %f %f\n", SOL[0], SOL[1], SOL[2]);
    plot_iterative_convergence(resvec, nbite, "Richardson");

    // Test Jacobi method
    printf("\nTest de la méthode de Jacobi...\n");
    memset(SOL, 0, la * sizeof(double));
    memset(resvec, 0, maxit * sizeof(double));
    extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    jacobi_iteration(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("Méthode de Jacobi : nombre d'itérations = %d\n", nbite);
    printf("Méthode de Jacobi : erreur avant relative = %e\n", relres);
    printf("Premières valeurs de la solution : %f %f %f\n", SOL[0], SOL[1], SOL[2]);
    plot_iterative_convergence(resvec, nbite, "Jacobi");

    // Test Gauss-Seidel method
    printf("\nTest de la méthode de Gauss-Seidel...\n");
    memset(SOL, 0, la * sizeof(double));
    memset(resvec, 0, maxit * sizeof(double));
    extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    gauss_seidel_iteration(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    relres = relative_forward_error(SOL, EX_SOL, &la);
    printf("Méthode de Gauss-Seidel : nombre d'itérations = %d\n", nbite);
    printf("Méthode de Gauss-Seidel : erreur avant relative = %e\n", relres);
    printf("Premières valeurs de la solution : %f %f %f\n", SOL[0], SOL[1], SOL[2]);
    plot_iterative_convergence(resvec, nbite, "GaussSeidel");

    // Free memory
    free(RHS);
    free(SOL);
    free(EX_SOL);
    free(X);
    free(resvec);
    free(AB);
    free(MB);

    return 0;
}
