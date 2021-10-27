#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, j, k;
    int max_Aii;
    double max;
    double* temp_Aii = (double*)malloc(n * sizeof(double)); // a temp row for swap
    
    for (i = 0; i < n; i++) {
        // find the biggest |A[i,i]|
        max = fabs(A[i*n + i]);
        max_Aii = i; 
        for (j = i+1; j < n; j++) {
            if ( fabs(A[j*n + i]) > max ) { // update the biggest element in i-th col
                max = fabs(A[j*n + i]);
                max_Aii = j; // |A[j,i]| is the biggest in rest of i-th col
            }
        }
        // swap row i and max_Aii and the corrosponding vector
        if (max == 0) {
            // printf("WARNING! A is singular, or near so... \n");
            return -1;
        } else {
            if (max_Aii != i) { // bug here!!!!
                // swap i-th and max_Aii-th element of vector pivot
                // printf("Before swap max_Aii and i... \n");
                // printf("index max-Aii = %d\n", max_Aii);
                // printf("ipiv[%d] = %d\n", max_Aii, ipiv[max_Aii]);
                // printf("index i = %d\n", i);
                // printf("ipiv[%d] = %d\n", i, ipiv[i]);
                int temp = ipiv[i];
                ipiv[i] = ipiv[max_Aii];
                ipiv[max_Aii] = temp;
                // printf("After swap max_Aii and i... \n");
                // printf("ipiv[%d] = %d\n", max_Aii, ipiv[max_Aii]);
                // printf("ipiv[%d] = %d\n", i, ipiv[i]);
                // swap i-th row and max_Aii-th row
                memcpy(temp_Aii, A + i*n, n*sizeof(double));
                memcpy(A + i*n, A + max_Aii*n, n*sizeof(double));
                memcpy(A + max_Aii*n, temp_Aii, n*sizeof(double));
            }
        }

        for (j = i+1; j < n; j++) {
            A[j*n + i] = A[j*n+i] / A[i*n + i];
            for (k = i+1; k < n; k++) {
                A[j*n + k] -= A[j*n + i] * A[i*n + k];
            }
        }
    }
    // printf("finally, ipiv is like this: \n");
    // int id;    
    // for (id = 0;id < n;id++) {
    //     printf("ipiv[%d] = %d\n", id, ipiv[id]);
    // }
    free(temp_Aii);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv) { // can only get result by copy to vector B
    /* add your code here */
    int i, j;
    double* y = (double*)malloc(n * sizeof(double)); // vector y
    
    if (UPLO == 'L') { // lower triangle, calculating vector y by L*y=b
        y[0] = B[ipiv[0]]; // y0*1 + y1*0 + y2*0+...+y(n-1)*0=b0
        for (i = 1;i < n;i++) {
            double sum = 0.0;
            for (j = 0; j < i; j++) {  // only use the lower triangle of A
                sum += A[i*n + j] * y[j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    } 
    if (UPLO == 'U') { // here is a bug!!!!!
        // for (j = n-1; j >= 0; --j) {
		//     temp = y[j];
		//     for (i = j - 2; i >= 0; --i) {
		//         y[i] -= temp * A(i,j);
		// 	}
		// }
    
        y[n-1] = B[n-1]/A[n*n-1];// y[n-1] is correct!!!
        for (i = n-2;i >= 0;i--) { // the other y[i] is incorrect!!!
            double sum = 0.0;
            for (j = i+1;j < n;j++) {
                sum += A[i*n + j] * y[j];
            }
            y[i] = (B[i] - sum)/A[i*n + i];
        }
    }
    memcpy(B, y, n*sizeof(double)); // copy back to vector B
    free(y);
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b) {
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    int i1, j1, k1;

    for (i1 = i; i1 < i+b; i1 += 2) {
        for (k1 = k; k1 < k+b; k1 += 2) {
            // reg for A
            register int ra1 = i1*n + k1;// A[i1, k1]
            register int ra2 = ra1 + n;// A[(i1+1), k1]
            register double a00 = A[ra1];// A[i1, k1]
            register double a01 = A[ra1+1];// A[i1, k1+1]
            register double a10 = A[ra2];// A[(i1+1), k1]
            register double a11 = A[ra2+1];// A[(i1+1), (k1+1)]
            for (j1 = j; j1 < j+b; j1 += 2) {
                // reg for C
                register int rc1 = i1*n + j1;// C[i1, j1]
                register int rc2 = rc1 + n;// C[i1+1, j1]
                register double c00 = C[rc1];// C[i1, j1]
                register double c01 = C[rc1+1];// C[i1, j1+1]
                register double c10 = C[rc2];// C[i1+1, j1]
                register double c11 = C[rc2+1];// C[i1+1, j1+1]
                // reg for B
                register int rb1 = k1*n + j1;// B[k1, j1]
                register int rb2 = rb1 + n;// B[(k1+1), j1]
                register double b00 = B[rb1];// B[k1, j1]
                register double b01 = B[rb1+1];// B[k1, j1+1]
                register double b10 = B[rb2];// B[(k1+1), j1]
                register double b11 = B[rb2+1];// B[k1+1, j1+1]
                // calculate C
                c00 -= a00 * b00 + a01 * b10;
                c01 -= a00 * b01 + a01 * b11;
                c10 -= a10 * b00 + a11 * b10;
                c11 -= a10 * b01 + a11 * b11;
                C[rc1] = c00;
                C[rc1+1] = c01;
                C[rc2] = c10;
                C[rc2+1] = c11;
            }
        }
    }
}


/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) {
    int i, j, k;
    int max_Aii;
    int ib;
    double max;
    double* temp_Aii = (double*)malloc(n * sizeof(double)); // a temp row for swap
    
    for (ib = 0;ib < n;ib += b) {
        for (i = ib; i < ib+b; i++) {
            max = fabs(A[i*n + i]);
            max_Aii = i;
            for (j = i+1;j < n;j++) {
                if ( fabs(A[j*n + i]) > max ) {
                    max = fabs(A[j*n + i]);
                    max_Aii = j; // |A[j,i]| is the biggest in rest of i-th col
                }
            }
            // swap row i and max_Aii and the corrosponding vector
            if (max == 0) {
                return -1;
            } else {
                if (max_Aii != i) { // bug here!!!!
                    // swap i-th and max_Aii-th element of vector pivot
                    int temp = ipiv[i];
                    ipiv[i] = ipiv[max_Aii];
                    ipiv[max_Aii] = temp;
                    // swap i-th row and max_Aii-th row
                    memcpy(temp_Aii, A + i*n, n*sizeof(double));
                    memcpy(A + i*n, A + max_Aii*n, n*sizeof(double));
                    memcpy(A + max_Aii*n, temp_Aii, n*sizeof(double));
                }
            }
            for (j = i+1; j < n; j++) {
                A[j*n + i] = A[j*n + i] / A[i*n + i];
                for (k = i+1; k < ib+b; k++) {
                    A[j*n + k] -= A[j*n + i] * A[i*n + k];
                }
            }
        }
        // A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n), update next b rows of U
        for (i = ib; i < ib+b; i++) {
            for (j = ib+b;j < n;j++) { // j:[end+1=ib+b-1+1, n]
                double sum = 0.0;
                for (k = ib;k < i;k++) {
                    sum += A[i*n + k] * A[k*n + j];
                }
                A[i*n + j] -= sum;
            }
        }

        // update A(end+1=ib+b-1+1 : n, end+1:n) block by block
        for (i = ib+b;i < n;i+=b) {
            for (j = ib+b;j < n;j+=b) {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }

    return 0;
}

