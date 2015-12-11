#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <fstream>
#include "solver.hpp"



// This is the function which solves the Jacobi method 
double jacobi_solve(double **matrix_first, double *vector_guess,double* vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int indx,k;
 double *guess_old, *temp_b, *residual_vector;
 double sum_squares, root_sum_squares;
 guess_old = (double *)malloc(nelem_in_array*sizeof(double));
 temp_b = (double *)malloc(nelem_in_array*sizeof(double));
 residual_vector = (double *)malloc(nelem_in_array*sizeof(double));

 
 for (indx = 0; indx < nelem_in_array; indx++)
 {
    guess_old[indx] = vector_guess[indx];
 }   

 for (indx = 0; indx < nelem_in_array; indx++)
 {
    double product_sum= 0.;
    for (k = 0; k < nelem_in_array; k++) 
    {
        if(k!=indx)
        {    
            product_sum = product_sum + matrix_first[indx][k]*guess_old[k];   // Implementing the equation of jacobi method
        }    
    }
    vector_guess[indx] =  1./ matrix_first[indx][indx] * ( vector_second[indx]- product_sum); 
 }
 // finding the reidual
 multiply_matrix(matrix_first, vector_guess, temp_b , nelem_in_array ) ;
 for( indx = 0; indx < nelem_in_array; indx++)
 {
    residual_vector[indx] = vector_second[indx] - temp_b[indx];
 }
 sum_squares=0.0;
 root_sum_squares=0.0;
 for( indx = 0; indx < nelem_in_array; indx++)
 {
    sum_squares= sum_squares +  residual_vector[indx]* residual_vector[indx];
 }
 root_sum_squares= sqrt(sum_squares);
 free(guess_old);
 free(residual_vector);
 free(temp_b);

return *vector_guess; 
} 

// This is the function which solves the Gauss- Seidal  method 

double gauss_seidal_solve(double **matrix_first, double *vector_guess,double* vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int indx,k;
 //double *vector_x_jacobi;
 //vector_x_jacobi = (double *)malloc(nelem_in_array*sizeof(double));
 //for (indx = 0; indx < nelem_in_array; indx++)
 //{
 //   vector_x_jacobi[indx]= vector_guess[indx];
 //}

 //jacobi_solve(matrix_first, vector_x_jacobi,vector_second, nelem_in_array);

 for (indx = 0; indx < nelem_in_array; indx++)
 {
    

    double product_sum= 0.;
    double product_sum2=0.;
    for (k = 0; k < nelem_in_array; k++) 
    {
        if(k<=indx-1)
        {    
            product_sum = product_sum + matrix_first[indx][k]*vector_guess[k];
        }   //product_sum = product_sum + matrix_first[indx][k]*vector_x_jacobi[k]; 
        if(k>=indx+1)
        {    
            product_sum2 = product_sum2 + matrix_first[indx][k]*vector_guess[k];
        }           
    }
    vector_guess[indx] =  1./ matrix_first[indx][indx] * ( vector_second[indx]- product_sum- product_sum2); 
 }
  
return *vector_guess; 
} 

// // This is the function which solves the Conjugate Gradient  method 


double get_steady_temperature_cg_solve(double **matrix_first, double *vector_guess,double* vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
 double *vector_r, *vector_p, *vector_cof,alpha_p, *vector_bk, *residual_vector,*temp_b, b_k_coeffic;
 double sum_squares, root_sum_squares;
 int indx,k;
 vector_r = (double *)malloc(nelem_in_array*sizeof(double));
 vector_p = (double *)malloc(nelem_in_array*sizeof(double));
 vector_cof = (double *)malloc(nelem_in_array*sizeof(double));
 vector_bk = (double *)malloc(nelem_in_array*sizeof(double));
 temp_b = (double *)malloc(nelem_in_array*sizeof(double));
 residual_vector = (double *)malloc(nelem_in_array*sizeof(double));
 

 std::ofstream output("./inputmatriks1.txt");
 output << "\n The vector guessed is \n";
 int value_sum=0;
 root_sum_squares= 1.0;
 while (root_sum_squares>= 0.00001)
 {
  value_sum= value_sum+1;
   if (value_sum==100000)
     break;
  std::cout << " I am inside the solver now hang on\n"; 
  multiply_matrix (matrix_first, vector_guess,vector_r, nelem_in_array);
 
  //#pragma omp parallel for
  for (k = 0; k < nelem_in_array; k++)
  {
    vector_r[k]= vector_second[k] - vector_r[k];
    vector_p[k]= vector_r[k];
  }  
 
  multiply_matrix(matrix_first, vector_r, vector_cof,nelem_in_array ) ;
  alpha_p = multiply_vector_transpose(vector_r,vector_r,nelem_in_array) / multiply_vector_transpose(vector_r,vector_cof, nelem_in_array) ;
  #pragma omp parallel for 
    for (k = 0; k < nelem_in_array; k++)
    {
        vector_guess[k] = vector_guess[k] + alpha_p * vector_p[k];
        multiply_matrix(matrix_first, vector_p, vector_cof,nelem_in_array );
        vector_r[k] = vector_r[k] - alpha_p * vector_cof[k]   ;
        multiply_matrix(matrix_first, vector_p, vector_bk,nelem_in_array ) ;
        b_k_coeffic = multiply_vector_transpose(vector_p,vector_r,nelem_in_array) / multiply_vector_transpose(vector_p,vector_bk, nelem_in_array) ;
    } 
 // finding the reidual
  multiply_matrix(matrix_first, vector_guess, temp_b , nelem_in_array ) ;
 // #pragma omp parallel for
  for( indx = 0; indx < nelem_in_array; indx++)
  {
    residual_vector[indx] = vector_second[indx] - temp_b[indx] ;
    output << vector_second[indx] << " " ;
    output << " \n";
    output << temp_b[indx] << " " ;
    output << "\n   ";
    output << residual_vector[indx]  << " " ;
  }
  sum_squares=0.0;
  root_sum_squares=0.0;
  for( indx = 0; indx < nelem_in_array; indx++)
  {
    sum_squares= sum_squares +  residual_vector[indx]* residual_vector[indx];
  }
  root_sum_squares= sqrt(sum_squares); 
  output << root_sum_squares << " " ;
 }
 
 free(vector_r);
 free(vector_p);
 free(vector_cof);
 free(vector_bk);
 free(temp_b);
 free(residual_vector);
return *vector_guess;
}

// // This is the function which multiplies a vector with a vector transpose and gives scalar back 

double multiply_vector_transpose(double *vector_first, double *vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int k;
 double result_vector;
 double sum = 0.;
    for (k = 0; k < nelem_in_array; k++) 
    {
        sum = sum + vector_first[k]*vector_second[k];
    }
    result_vector = sum;
 
return result_vector;
}


// // This is the function which multiplies a matrix with a vector and gives back a vector  

double multiply_matrix(double **matrix_first, double *vector_first, double *vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int indx,k;
 
 for (indx = 0; indx < nelem_in_array; indx++)
 {
    double sum = 0;
    for (k = 0; k < nelem_in_array; k++) 
    {
        sum = sum + matrix_first[indx][k]*vector_first[k];
    }
    vector_second[indx] = sum;
 }
return *vector_second;
} 



