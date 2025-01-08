/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    // Pour une matrice tridiagonale en format bande colonne majeure :
    // - la diagonale principale est stockée dans la ligne kl+ku
    // - la sous-diagonale est stockée dans la ligne kl+ku-1
    // - la sur-diagonale est stockée dans la ligne kl+ku+1
    int i, j;
    int kl = 1;  // nombre de sous-diagonales
    int ku = 1;  // nombre de sur-diagonales
    
    // Initialisation à zéro
    for(i = 0; i < *lab * (*la); i++){
        AB[i] = 0.0;
    }
    
    // Remplissage de la diagonale principale (2.0)
    for(i = 0; i < *la; i++){
        AB[(*kv+1) + i*(*lab)] = 2.0;
    }
    
    // Remplissage de la sous-diagonale (-1.0)
    for(i = 1; i < *la; i++){
        AB[*kv + i*(*lab)] = -1.0;
    }
    
    // Remplissage de la sur-diagonale (-1.0)
    for(i = 0; i < *la-1; i++){
        AB[(*kv+2) + i*(*lab)] = -1.0;
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
	double norme_x = sqrt(cblas_ddot(*la, x, 1, x, 1));
	cblas_daxpy(*la, -1, y, 1, x, 1);
	double res_norme = sqrt(cblas_ddot(*la, x, 1, x, 1));
	return res_norme / norme_x;
}
int indexABCol(int i, int j, int *lab){
  return (j + 1) * (*lab - 1) + i - 1;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
