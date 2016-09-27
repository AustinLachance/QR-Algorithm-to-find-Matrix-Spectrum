/*****************************************************************
Name: Austin Lachance                                                            
Email: austin.lachance@yale.edu                                                  
Date: 04/15/16                                                             
                                                                                 
CPSC440                                                                          
QR Algorithm with Shifts                                               
Description: This program finds the spectrum of a tridiagonal
matrix using the QR Algorithm with shifts                                                                
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void *transpose(double *matrix, int n);
void *identity(double* matrix, int n);
void matrix_print(double *matrix, int n);
int error_check(double a, double b, double error);
void qr_symmetric(double * a, int n, double * b);
void matrix_copy(double *matrix1, double* matrix2, int n);
void matrix_shift(double* matrix, double shift_val, int n, int add);
void sub_matrix(double* b, int n, int new_size, double* b_new);
void givens_matrix(double* b, int n, int j, double *q);
void row_update(double* b_new, int n, int j, double* q, int adjoint);


int main() {
	int n = 4;
	double* a = malloc(n*n*sizeof(double));
	double* b = malloc(n*n*sizeof(double));

	a[0] = 3;
	a[1] = 4;
	a[2] = 0;
	a[3] = 0;
	a[4] = 4;
	a[5] = 2;
	a[6] = 1;
	a[7] = 0;
	a[8] = 0;
	a[9] = 1;
	a[10] = 6;
	a[11] = 3;
	a[12] = 0;
	a[13] = 0;
	a[14] = 3;
	a[15] = 1;


	printf("A:\n");
	matrix_print(a, n);

	qr_symmetric(a, n, b);

	

	printf("B:\n");
	matrix_print(b, n);
	return 0;
}


void qr_symmetric(double * a, int n, double * b) {
	int i, j, row, col, new_n;
	double shift_val;

	//copy matrix a into b
	matrix_copy(a, b, n);

	//Calculate the eigenvalue at position (i, i)
	for(i = n; i > 1; i--) {
		double* b_new = malloc(i*i*sizeof(double));
		sub_matrix(b, n, i, b_new);

		while(!error_check(b_new[1], 0.000, 1e-30)) {

			double* old_b = malloc(i*i*sizeof(double));

			matrix_copy(b_new, old_b, i);
			shift_val = b_new[0];

			matrix_shift(b_new, shift_val, i, 0);
			matrix_shift(old_b, shift_val, i, 0);

			//Perform Givens Rotation on b_new to create L 
			for(j = i-2; j >= 0; j--) {
				double* q = malloc(i*i*sizeof(double));
				givens_matrix(b_new, i, j, q);
				row_update(b_new, i, j, q, 0);
				free(q);
			}
			//Apply adjoints of Givens Rotation to L
			for(j = i-2; j >= 0; j--) {

				double* q = malloc(i*i*sizeof(double));

				givens_matrix(old_b, i, j, q);
				row_update(old_b, i, j, q, 0);

				transpose(q, i);
				row_update(b_new, i, j, q, 1);

				free(q);
			}
			free(old_b);
			matrix_shift(b_new, shift_val, i, 1);
			
		}
		
		//Update values in b based on b_new
		for(row = n-i; row < n; row++) {
			for(col = n-i; col < n; col++) {
				b[row*n + col] = b_new[(row - (n-i))*i+ col - (n-i)];
			}
		}
		free(b_new);
		

	}

}

//Creates a givens rotation matrix, q
void givens_matrix(double* b, int n, int j, double *q) {
	double c, s, x, y, r;
	x = b[j*n + (j+1)];
	y = b[(j + 1)*n + (j+1)];
	r = sqrt(x*x + y*y);
	c = y/r;
	s = x/r;

	q[j*n + j] = c;
	q[j*n + (j + 1)] = -s;
	q[(j + 1)*n + j] = s;
	q[(j + 1)*n + (j+1)] = c;

}	

//Transposes matrix
void *transpose(double *matrix, int n) {
  double *transpose = (double*)malloc(n*n*sizeof(double));
  int i, j = 0;
  for(i = 0; i < n; i++) {
  	for(j = 0; j < n; j++) {
	  transpose[i*n + j] = matrix[j*n + i];
	}
  }
  for(i = 0; i < n*n; i++) {
  	matrix[i] = transpose[i];
  }
  free(transpose);
}

//Makes matrix into the identity
void *identity(double* matrix, int n) {
  int i = 0;
  for(i = 0; i < n; i++) 
    {
      matrix[i*n + i] = (double)1;
    }
}

//Prints a 1-Dimensional array as a square matrix of size n
void matrix_print(double *matrix, int n) {
  int i, j = 0;
  for(i = 0; i < n; i++) {
  	for(j = 0; j < n; j++) {
	  printf("%e\t", matrix[i*n + j]);
	}
	printf("\n");
  }
}

//Checks if the difference between a and b falls within a specified error range
int error_check(double a, double b, double error) {
  return (fabs(a-b) < error);
}


//Copies each element in matrix1 into matrix2
void matrix_copy(double *matrix1, double* matrix2, int n) {
	int i;
	for (i = 0; i < n*n; i++) {
		matrix2[i] = matrix1[i];
	}
}

//Shifts the diagonal elements of matrix by shift_val
void matrix_shift(double* matrix, double shift_val, int n, int add) {
	int row, col;
	for(row = 0; row < n; row++) {
		for(col = 0; col < n; col++) {
			if(row == col) {
				if(add) {
					matrix[row*n + col] += shift_val;
				}
				else {
					matrix[row*n + col] -= shift_val;
				}
			}
		}
	}
}

//Create b_new as submatrix of b missing first (n - new_size) rows & cols
void sub_matrix(double* b, int n, int new_size, double* b_new) {
	int i, j;
	for(i = 0; i < new_size; i++) {
		for(j = 0; j < new_size; j++) {
			b_new[i*new_size + j] = b[(n - new_size + i)*n + (n - new_size + j)];
		}
	}
}


//Updates rows in b-new by multiplying by givens rotation matrix
void row_update(double* b_new, int n, int j, double* q, int adjoint) {
	int i;
	double* temp = malloc(2*n*sizeof(double));

	if(adjoint) {
		for(i = 0; i < n; i++) {
			temp[i] = b_new[i*n + j] * q[j*n + j]
			 + b_new[i*n + (j+1)] * q[(j+1)*n + j];

			temp[i + n] = b_new[i*n + j] * q[j*n + (j+1)]
			 + b_new[i*n + (j+1)] * q[(j+1)*n + (j+1)];
		}

		for(i = 0; i < n; i++) {
			b_new[i*n + j] = temp[i];
			b_new[i*n + (j+1)] = temp[i + n];
		}
	}
	else {
		for(i = 0; i < n; i++) {
			temp[i] = q[j*n + j] * b_new[j*n + i]
			 + q[j*n + (j+1)] * b_new[(j+1)*n + i];

			temp[i + n] = q[(j+1)*n + j] * b_new[j*n + i]
			 + q[(j+1)*n + (j+1)] * b_new[(j+1)*n + i];
		}

		for(i = 0; i < n; i++) {
			b_new[j*n + i] = temp[i];
			b_new[(j+1)*n + i] = temp[i + n];
		}
	}

	free(temp);
}



