/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}

void plot_richardson_convergence(double* error_vec, int nb_iter, char* method_name) {
    char data_filename[256];
    snprintf(data_filename, sizeof(data_filename), "richardson_%s_data.dat", method_name);
    
    // Écriture des données
    FILE* data_file = fopen(data_filename, "w");
    if (!data_file) {
        perror("Erreur lors de l'ouverture du fichier de données");
        return;
    }
    
    for (int i = 0; i < nb_iter; i++) {
        fprintf(data_file, "%d %.10e\n", i, error_vec[i]);
    }
    fclose(data_file);
    
    // Création du script gnuplot
    char script_filename[256];
    snprintf(script_filename, sizeof(script_filename), "richardson_%s.gnuplot", method_name);
    
    FILE* script_file = fopen(script_filename, "w");
    if (!script_file) {
        perror("Erreur lors de l'ouverture du fichier script");
        return;
    }

    fprintf(script_file, 
        "set terminal png size 1200,800\n"
        "set output 'richardson_%s.png'\n"
        "set title 'Convergence de la methode de Richardson'\n"
        "set xlabel 'Iterations'\n"
        "set ylabel 'Erreur Avant'\n"
        "set grid\n"
        "set format y '%%g'\n"
        "plot '%s' using 1:2 with linespoints title 'Erreur'\n",
        method_name, data_filename);
    
    fclose(script_file);
    
    // Exécution de gnuplot
    char command[512];
    snprintf(command, sizeof(command), "gnuplot %s", script_filename);
    int res = system(command);
    if (res) {
        perror("Erreur lors de l'exécution de gnuplot");
    } else {
        printf("Graphique généré : richardson_%s.png\n", method_name);
    }
    
    // Nettoyage
    remove(data_filename);
    remove(script_filename);
}