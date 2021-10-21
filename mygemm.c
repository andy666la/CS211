#include "mygemm.h"

/**
 *
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 *
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i,j,k;
    for (i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                    for(k = 0; k < n; k++){
                        C[i*n+j] += A[i*n+k] * B[k*n+j];
                    }
            }
        }
}

void dgemm1(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
        for (i = 0; i < n; i++){
               for(j = 0; j < n; j++){
                double r =C[i*n+j];
                    for(k = 0; k < n; k++){
                        r += A[i*n+k] * B[k*n+j];
                        C[i*n+j] = r;
                    }
            }
        }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i+=2){
        for(j = 0; j < n; j+=2){
            for(k = 0; k < n; k+=2){
                C[i*n+j] = A[i*n+k]*B[k*n+j]
                         + A[i*n+k+1]*B[(k+1)*n+j]
                         + C[i*n+j];
                C[(i+1)*n+j] = A[(i+1)*n + k]*B[k*n+j]
                             + A[(i+1)*n + k+1]*B[(k+1)*n+j]
                             + C[(i+1)*n+j];
                C[i*n + (j+1)] = A[i*n + k]*B[k*n + (j+1)]
                               + A[i*n + k+1]*B[(k+1)*n + (j+1)]
                               + C[i*n + (j+1)];
                C[(i+1)*n + (j+1)] = A[(i+1)*n + k]*B[k*n + (j+1)]
                                   + A[(i+1)*n + k+1]*B[(k+1)*n + (j+1)]
                                   + C[(i+1)*n + (j+1)];
            }
        }
    }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i+=2){
        for(j = 0; j < n; j+=2){
            double C_ij = C[i*n+j],
                   C_i1j = C[(i+1)*n+j],
                   C_ij1 = C[i*n + (j+1)],
                   C_i1j1 =  C[(i+1)*n + (j+1)];
            for(k = 0; k < n; k+=2){
                double A_ik = A[i*n+k],
                       A_ik1 = A[i*n+k+1],
                       A_i1k = A[(i+1)*n + k],
                       A_i1k1 = A[(i+1)*n + k+1],
                       B_kj = B[k*n+j],
                       B_k1j = B[(k+1)*n+j],
                       B_kj1 = B[k*n + (j+1)],
                       B_k1j1 = B[(k+1)*n + (j+1)];
                C_ij   = A_ik*B_kj   + A_ik1*B_k1j   + C_ij;
                C_i1j  = A_i1k*B_kj  + A_i1k1*B_k1j  + C_i1j;
                C_ij1  = A_ik*B_kj1  + A_ik1*B_k1j1  + C_ij1;
                C_i1j1 = A_i1k*B_kj1 + A_i1k1*B_k1j1 + C_i1j1;
            }
            C[i*n+j] = C_ij;
            C[(i+1)*n+j] = C_i1j;
            C[i*n + (j+1)] = C_ij1;
            C[(i+1)*n + (j+1)] =  C_i1j1;
        }
    }


}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            double sum = C[i*n+j];
            for(k = 0; k < n; k++){
                sum += A[i*n+k]*B[k*n+j];
            }
            C[i*n+j] = sum;
        }
}
}

void bijk(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(i = 0; i < n; i+=b){
        for(j = 0; j < n; j+=b){
            for(k = 0; k < n; k+=b){
                for(i1 = i; i1 < i+b; i1++){
                    for(j1 = j; j1 < j+b; j1++){
                        double sum = C[i1*n+j1];
                        for(k1 = k; k1 < k+b; k1++)
                            sum += A[i1*n+k1]*B[k1*n+j1];
                        C[i1*n+j1] = sum;
                    }
                }
            }
        }
    }
}

void jik(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
for(j = 0; j < n; j++){
    for(i = 0; i < n; i++){
        double sum = C[i*n+j];
        for(k = 0; k < n; k++){
            sum += A[i*n+k]*B[k*n+j];
        }
        C[i*n+j] = sum;
    }
}
}

void bjik(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(j = 0; j < n; j+=b){
        for(i = 0; i < n; i+=b){
            for(k = 0; k < n; k+=b){
                for(j1 = j; j1 < j+b; j1++){
                    for(i1 = i; i1 < i+b; i1++){
                        double sum = C[i1*n+j1];
                        for(k1 = k; k1 < k+b; k1++)
                            sum += A[i1*n+k1]*B[k1*n+j1];
                        C[i1*n+j1] = sum;
                    }
                }
            }
        }
    }
}

void kij(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
for(k = 0; k < n; k++){
    for(i = 0; i < n; i++){
        double r = A[i*n+k];
        for(j = 0; j < n; j++){
            C[i*n+j] += r*B[k*n+j];
        }
    }
}
}

void bkij(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(k = 0; k < n; k+=b){
        for(i = 0; i < n; i+=b){
            for(j = 0; j < n; j+=b){
                for(k1 = k; k1 < k+b; k1++){
                    for(i1 = i; i1 < i+b; i1++){
                        double r = A[i1*n+k1];
                        for(j1 = j; j1 < j+b; j1++)
                            C[i1*n+j1] += r*B[k1*n+j1];
                    }
                }
            }
        }
    }
}
    

void ikj(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
for(i = 0; i < n; i++){
    for(k = 0; k < n; k++){
        double r = A[i*n+k];
        for(j = 0; j < n; j++){
            C[i*n+j] += r*B[k*n+j];
        }
    }
}
}

void bikj(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(i = 0; i < n; i+=b){
        for(k = 0; k < n; k+=b){
            for(j = 0; j < n; j+=b){
                for(i1 = i; i1 < i+b; i1++){
                    for(k1 = k; k1 < k+b; k1++){
                        double r = A[i1*n+k1];
                        for(j1 = j; j1 < j+b; j1++){
                            C[i1*n+j1] += r*B[k1*n+j1];
                        }
                    }
                }
            }
        }
    }
}
void jki(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
for(j = 0; j < n; j++){
    for(k = 0; k < n; k++){
        double r = B[k*n+j];
        for(i = 0; i < n; i++){
            C[i*n+j] += A[i*n+k]*r;
        }
    }
}
}

void bjki(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(j = 0; j < n; j+=b){
        for(k = 0; k < n; k+=b){
            for(i = 0; i < n; i+=b){
                for(j1 = j; j1 < j+b; j1++){
                    for(k1 = k; k1 < k+b; k1++){
                        double r = B[k1*n+j1];
                        for(i1 = i; i1 < i+b; i1++)
                            C[i1*n+j1] += A[i1*n+k1]*r;
                    }
                }
            }
        }
    }
}

void kji(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
for(k = 0; k < n; k++){
    for(j = 0; j < n; j++){
        double r = B[k*n+j];
        for(i = 0; i < n; i++){
            C[i*n+j] += A[i*n+k]*r;
        }
    }
}
}

void bkji(const double *A, const double *B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(k = 0; k < n; k+=b){
        for(j = 0; j < n; j+=b){
            for(i = 0; i < n; i+=b){
                for(k1 = k; k1 < k+b; k1++){
                    for(j1 = j; j1 < j+b; j1++){
                        double r = B[k1*n+j1];
                        for(i1 = i; i1 < i+b; i1++){
                            C[i1*n+j1] += A[i1*n+k1]*r;
                        }
                    }
                }
            }
        }
    }
}
//Cache Reuse part 3 End

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1, j1, k1;
    for(k = 0; k < n; k +=b)
    {
        for(i = 0; i < n; i +=b)
        {
            for(j = 0; j < n; j +=b)
            {
                /* BxB mini matrix multiplications */
                for(k1 = k; k1 < k+b; k1+=2)
                {
                    int k_0 = (k1+1);
                    int k_1 = k1;
                    for(i1 = i; i1 < i+b; i1+=2)
                    {
                        int i_0 = i1*n;
                        int i_1 = (i1+1)*n;
                        double a0 = A[i_0+k_1];
                        double a1 = A[i_0+k_0];
                        double a2 = A[i_1+k_1];
                        double a3 = A[i_1+k_0];
                        for(j1 = j; j1 < j+b; j1+=2)
                        {
                            double c0 = C[i_0 +j1];
                            double c1 = C[i_1+j1];
                            double c2 = C[i_0+(j1+1)];
                            double c3 = C[i_1+(j1+1)];
                            double b0 = B[k_1*n+j1];
                            double b1 = B[k_0*n+j1];
                            double b2 = B[k_1*n+(j1+1)];
                            double b3 = B[k_0*n+(j1+1)];
                            c0 = a0 * b0 + a1 * b1 + c0;
                            c1 = a2 * b0 + a3 * b1 + c1;
                            c2 = a0 * b2 + a1 * b3 + c2;
                            c3 = a2 * b2 + a3 * b3 + c3;
                            C[i_0 +j1] = c0;
                            C[i_1+j1] = c1;
                            C[i_0+(j1+1)] = c2;
                            C[i_1+(j1+1)] = c3;
                        }
                    }
                }
            }
        }
    }
}

