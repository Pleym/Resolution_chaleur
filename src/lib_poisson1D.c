/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "lapacke.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// construit la matrice trigonale en format GB 
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv) {
    int n = *la;
    int k = *kv;
    
    for (int i = 0; i < (*lab) * n; i++) {
        AB[i] = 0.0;
    }
    
    for (int j = 0; j < n; j++) {
        AB[k + j*(*lab)] = 2.0;
        
        if (j > 0) {
            AB[(k-1) + j*(*lab)] = -1.0;
        }
        
        if (j < n-1) {
            AB[(k+1) + j*(*lab)] = -1.0;
        }
    }
}

// créer une matrice identité en bande 
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv) {
    int n = *la;
    for (int i = 0; i < n; ++i) {
        AB[*lab + i] = 1.0;
        *lab += *kv; 
    }
}

// définit le second membre du vecteur aux conditions limites de Dirichlet 
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1) {
    int n = *la;
    
    for (int i = 0; i < n; i++) {
        RHS[i] = 0.0;
    }
    
    RHS[0] = *BC0;
    RHS[n-1] = *BC1;
}

// donne la solution analytique
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1) {
    int n = *la;
    double h = 1.0 / (n + 1);
    
    for (int i = 0; i < n; i++) {
        double x = (i + 1) * h;
        EX_SOL[i] = *BC0 * (1.0 - x) + *BC1 * x;
    }
}

// définit les points de la grille
void set_grid_points_1D(double* x, int* la) {
    int n = *la; 
    double h = 1.0 / (n + 1); 
    
    for (int i = 0; i < n; ++i) {
        x[i] = (i + 1) * h; 
    }
}

double relative_forward_error(double* x, double* y, int* la) {
    double norm_diff = 0.0, norm_y = 0.0; 
    
    for(int i = 0; i < *la; i++) {
        norm_diff += (x[i] - y[i]) * (x[i] - y[i]);
        norm_y += y[i] * y[i];
    }
    
    return sqrt(norm_diff) / sqrt(norm_y);
}

int indexABCol(int i, int j, int *lab) {
    return 0;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) {
    // Vérification des conditions pour une matrice tridiag
    if (*lab != 4 || *kl != 1 || *ku != 1 || *la < *n) {
        *info = -1;
        return *info;
    }
    
    *info = 0;
    
    // Factorisation LU pour matrice tridiagonale
    for (int k = 0; k < *n - 1; k++) {
        double diag = fabs(AB[indexABCol(k, k, lab)]);
        double subdiag = fabs(AB[indexABCol(k + 1, k, lab)]);
        
        if (subdiag > diag) {
            // Échange des lignes k et k+1
            for (int j = k; j <= (k + 2 < *n ? k + 2 : *n - 1); j++) {
                double temp = AB[indexABCol(k, j, lab)];
                AB[indexABCol(k, j, lab)] = AB[indexABCol(k + 1, j, lab)];
                AB[indexABCol(k + 1, j, lab)] = temp;
            }
            ipiv[k] = k + 2;  
        } else {
            ipiv[k] = k + 1;  
        }
        
        // Vérification pivot nul
        if (AB[indexABCol(k, k, lab)] == 0.0) {
            *info = k + 1;
            return *info;
        }
        
        // Élimination de Gauss
        double multiplier = AB[indexABCol(k + 1, k, lab)] / AB[indexABCol(k, k, lab)];
        AB[indexABCol(k + 1, k, lab)] = multiplier;  // Stockage du multiplicateur dans L
        AB[indexABCol(k + 1, k + 1, lab)] -= multiplier * AB[indexABCol(k, k + 1, lab)];
    }
    
    // Dernier pivot
    ipiv[*n - 1] = *n;
    if (AB[indexABCol(*n - 1, *n - 1, lab)] == 0.0) {
        *info = *n;
    }
    
    return *info;
}