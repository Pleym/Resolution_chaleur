/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"


extern void dgbtrf_(int* m, int* n, int* kl, int* ku, double* ab, int* ldab, int* ipiv, int* info);
extern void dgbtrs_(char* trans, int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);

void eig_poisson1D(double* eigval, int *la){

  return -2.0*cos(M_PI/(*la+1)) +2.0;
}

double eigmax_poisson1D(int *la){

  double h = 1.0 / ((*la) + 1.0);
  double sin_theta = sin(*la * M_PI * h / 2.0);
  return 4.0 * sin_theta * sin_theta;
}

double eigmin_poisson1D(int *la){

  double h = 1.0 / (*la + 1.0);
  double sin_theta = sin(M_PI * h / 2.0);
  return 4.0 * sin_theta * sin_theta;
}

double richardson_alpha_opt(int *la){
  double min = -2.0*cos(M_PI/(*la+1)) + 2.0;
  double max = -2.0*cos((*la)*M_PI/(*la+1)) + 2.0;
  return 2.0/(min + max);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, 
                     int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
  int nb_iterations = 0;
  double norm;
  double *temp_vec;
  int i;
  
  temp_vec = (double *) malloc(sizeof(double)*(*la));
  
  norm = cblas_dnrm2(*la, RHS, 1);
  norm = 1.0/norm;

  memcpy(temp_vec, RHS, sizeof(double)*(*la));

  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp_vec, 1);

  while (nb_iterations < *maxit && cblas_dnrm2(*la, temp_vec, 1)*norm > *tol){
    cblas_daxpy(*la, *alpha_rich, temp_vec, 1, X, 1);
    memcpy(temp_vec, RHS, sizeof(double)*(*la));
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp_vec, 1);
    resvec[nb_iterations] = cblas_dnrm2(*la, temp_vec, 1)*norm;
    nb_iterations++;
  }

  *nbite = nb_iterations;
  free(temp_vec);
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
                  int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite){
  int nb_iterations = 0;
  double *temp_vec;
  double *temp_x;
  int i;
  int j;
  double temp_sum;
  double temp_value;
  double temp_alpha;
  
  temp_vec = (double *) malloc(sizeof(double)*(*la));
  temp_x = (double *) malloc(sizeof(double)*(*la));
  
  while (nb_iterations < *maxit){
    memcpy(temp_x, X, sizeof(double)*(*la));
    
    for (i=0; i<*la; i++){
      temp_sum = 0.0;
      for (j=0; j<*la; j++){
        if (j != i){
          temp_value = AB[(*kl+j-i)+i*(*lab)];
          temp_sum += temp_value*X[j];
        }
      }
      temp_alpha = MB[i];
      X[i] = (1.0-temp_alpha)*X[i] + temp_alpha*(RHS[i]-temp_sum)/AB[*kl+i*(*lab)];
    }
    
    cblas_dcopy(*la, RHS, 1, temp_vec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, temp_vec, 1);
    resvec[nb_iterations] = cblas_dnrm2(*la, temp_vec, 1);
    
    if (resvec[nb_iterations] < *tol){
      break;
    }
    nb_iterations++;
  }
  
  *nbite = nb_iterations + 1;
  free(temp_vec);
  free(temp_x);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv){
  for (int i = 0; i < *la; i++) {
    for (int j = 0; j < *lab; j++) {
      int ind = (*lab)*i;
      MB[ind + j] = 0.0;
    }
    int ind = ((*lab)*i) + (*ku);
    MB[ind] = AB[ind];
  }
}


void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv){
  for (int i=0; i<*la; i++){
    MB[*lab * i + *ku] = AB[i * (*lab) + *ku];
  }
}