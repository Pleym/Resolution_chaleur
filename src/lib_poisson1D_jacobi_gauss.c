#include "lib_poisson1D.h"
#include <string.h>

double jacobi_alpha_opt(int *la) {
    double lambda_min = eigmin_poisson1D(la);
    double lambda_max = eigmax_poisson1D(la);
    return 2.0 / (lambda_min + lambda_max);
}

double gauss_seidel_alpha_opt(int *la) {
    double lambda_min = eigmin_poisson1D(la);
    double lambda_max = eigmax_poisson1D(la);
    return 2.0 / (lambda_min + lambda_max);
}

void jacobi_iteration(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                     int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int it;
    double *resvec_temp = (double *)calloc(*la, sizeof(double));
    double *X_temp = (double *)calloc(*la, sizeof(double));
    double *EX_SOL = (double *)calloc(*la, sizeof(double));
    double T0 = -5.0, T1 = 5.0;
    double *grid = (double *)calloc(*la, sizeof(double));
    double alpha_opt = jacobi_alpha_opt(la);
    
    printf("Alpha optimal pour Jacobi : %f\n", alpha_opt);
    
    // Setup grid and exact solution for error computation
    set_grid_points_1D(grid, la);
    set_analytical_solution_DBC_1D(EX_SOL, grid, la, &T0, &T1);
    
    for (it = 0; it < *maxit; it++) {
        // Save current X
        cblas_dcopy(*la, X, 1, X_temp, 1);
        
        // Jacobi iteration
        for (int i = 0; i < *la; i++) {
            double diag = AB[(*kl+1) + i*(*lab)];  // Diagonal element
            double sum = RHS[i];
            
            // Lower diagonal contribution
            if (i > 0) {
                sum -= AB[*kl + i*(*lab)] * X_temp[i-1];
            }
            
            // Upper diagonal contribution
            if (i < *la-1) {
                sum -= AB[(*kl+2) + i*(*lab)] * X_temp[i+1];
            }
            
            // Update X using Jacobi avec relaxation
            X[i] = (1.0 - alpha_opt) * X_temp[i] + (alpha_opt / diag) * sum;
        }
        
        // Compute relative forward error for this iteration
        cblas_dcopy(*la, X, 1, resvec_temp, 1);
        resvec[it] = relative_forward_error(resvec_temp, EX_SOL, la);
        
        // Check convergence
        if (resvec[it] < *tol) {
            it++;
            break;
        }
    }
    
    *nbite = it;
    free(resvec_temp);
    free(X_temp);
    free(EX_SOL);
    free(grid);
}

void gauss_seidel_iteration(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                           int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) {
    int it;
    double *resvec_temp = (double *)calloc(*la, sizeof(double));
    double *X_old = (double *)calloc(*la, sizeof(double));
    double *EX_SOL = (double *)calloc(*la, sizeof(double));
    double T0 = -5.0, T1 = 5.0;
    double *grid = (double *)calloc(*la, sizeof(double));
    double alpha_opt = gauss_seidel_alpha_opt(la);
    
    printf("Alpha optimal pour Gauss-Seidel : %f\n", alpha_opt);
    
    // Setup grid and exact solution for error computation
    set_grid_points_1D(grid, la);
    set_analytical_solution_DBC_1D(EX_SOL, grid, la, &T0, &T1);
    
    for (it = 0; it < *maxit; it++) {
        // Save old X for relaxation
        cblas_dcopy(*la, X, 1, X_old, 1);
        
        // Gauss-Seidel iteration
        for (int i = 0; i < *la; i++) {
            double diag = AB[(*kl+1) + i*(*lab)];  // Diagonal element
            double sum = RHS[i];
            
            // Lower diagonal contribution (using updated values)
            if (i > 0) {
                sum -= AB[*kl + i*(*lab)] * X[i-1];
            }
            
            // Upper diagonal contribution (using old values)
            if (i < *la-1) {
                sum -= AB[(*kl+2) + i*(*lab)] * X[i+1];
            }
            
            // Update X using Gauss-Seidel avec relaxation
            X[i] = (1.0 - alpha_opt) * X_old[i] + (alpha_opt / diag) * sum;
        }
        
        // Compute relative forward error for this iteration
        cblas_dcopy(*la, X, 1, resvec_temp, 1);
        resvec[it] = relative_forward_error(resvec_temp, EX_SOL, la);
        
        // Check convergence
        if (resvec[it] < *tol) {
            it++;
            break;
        }
    }
    
    *nbite = it;
    free(resvec_temp);
    free(X_old);
    free(EX_SOL);
    free(grid);
}

void plot_iterative_convergence(double* error_vec, int nb_iter, char* method_name) {
    char data_filename[256];
    char gnuplot_filename[256];
    char plot_filename[256];
    
    // Create filenames
    sprintf(data_filename, "convergence_%s.dat", method_name);
    sprintf(gnuplot_filename, "plot_%s.gnu", method_name);
    sprintf(plot_filename, "convergence_%s.png", method_name);
    
    // Write data to file
    FILE* data_file = fopen(data_filename, "w");
    if (data_file != NULL) {
        for (int i = 0; i < nb_iter; i++) {
            fprintf(data_file, "%d %e\n", i, error_vec[i]);
        }
        fclose(data_file);
        
        // Create gnuplot script
        FILE* gnuplot_file = fopen(gnuplot_filename, "w");
        if (gnuplot_file != NULL) {
            fprintf(gnuplot_file, "set terminal png\n");
            fprintf(gnuplot_file, "set output '%s'\n", plot_filename);
            fprintf(gnuplot_file, "set title 'Convergence de la méthode %s'\n", method_name);
            fprintf(gnuplot_file, "set xlabel 'Nombre d''itérations'\n");
            fprintf(gnuplot_file, "set ylabel 'Erreur avant relative'\n");
            if (strcmp(method_name, "Jacobi") == 0) {
                fprintf(gnuplot_file, "plot '%s' using 1:2 with lines linecolor rgb '#0000FF' linewidth 2 title 'Méthode %s'\n", data_filename, method_name);
            } else {
                fprintf(gnuplot_file, "plot '%s' using 1:2 with lines linecolor rgb '#FF0000' linewidth 2 title 'Méthode %s'\n", data_filename, method_name);
            }
            fclose(gnuplot_file);
            
            // Execute gnuplot script
            char command[512];
            sprintf(command, "gnuplot %s", gnuplot_filename);
            system(command);
        }
    }
}
