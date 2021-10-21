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
        for (i=0; i<n; i++)
                for (j=0; j<n; j++)
                        for (k=0; k<n; k++)
                                C[i*n+j] += A[i*n+k] * B[k*n+j];
}

void dgemm1(const double *A, const double *B, double *C, const int n)
{
        int i,j,k;
        for (i=0; i<n; i++)
                for (j=0; j<n; j++) {
                        register double r = C[i*n+j] ;
                         for (k=0; k<n; k++)
                                r += A[i*n+k] * B[k*n+j];
                         C[i*n+j] = r;
                        }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n)
{
int i=0;int j=0;int k=0;
         for (i=0;i<n;i+=2)
                for (j=0;j<n;j+=2)
                     for (k=0;k<n;k+=2){
                        C[i*n + j]= A[i*n +k]*B[k*n + j] + A[i*n + k+1 ]* B[(k+1)*n + j] + C[i*n + j];
                        C[(i+1)*n+j]=A[(i+1)*n+k]*B[k*n+j]+A[(i+1)*n+k+1]*B[(k+1)*n+j]+C[(i+1)*n+j];
                        C[i*n + (j+1)]= A[i*n + k]*B[k*n + (j+1)] + A[i*n + k+1]*B[(k+1)*n + (j+1)] + C[i*n + (j+1)];
                        C[(i+1)*n + (j+1)]= A[(i+1)*n + k]*B[k*n + (j+1)] + A[(i+1)*n + k+1]*B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
                        }
}
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n)
{
        int i,j,k;
for(i = 0; i < n; i += 2)
       for(j = 0; j < n; j += 2)  {
           register int t = i*n+j; register int tt = t+n;
           register double c00 = C[i*n+j]; register double c01 = C[i*n+j +1];  register double c10 = C[tt]; register double c11 = C[tt+1];

            for(k = 0; k < n; k += 2) {
                /* 2 by 2 mini matrix multiplication using registers*/
                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                register double a00 = A[ta]; register double a01 = A[ta+1]; register double a10 = A[tta]; register double a11 = A[tta+1];
                register double b00 = B[tb]; register double b01 = B[tb+1]; register double b10 = B[ttb]; register double b11 = B[ttb+1];
                c00 += a00*b00 + a01*b10;
                c01 += a00*b01 + a01*b11;
                c10 += a10*b00 + a11*b10;
                c11 += a10*b01 + a11*b11;
             }

             C[t] = c00;
             C[t+1] = c01;
             C[tt] = c10;
             C[tt+1] = c11;
        }

}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n)
{
int i,j,k;
for (i=0; i<n; i++)  {
  for (j=0; j<n; j++) {
    register double sum = C[i*n+j];
    for (k=0; k<n; k++)
      sum += A[i*n+k] * B[k*n+j];
    C[i*n+j] = sum;
  }
}

}

void bijk(const double *A, const double *B, double *C, const int n, const int b)
{
int i1,j1,k1;
int i,j,k;
for (i1 = 0; i1 < n; i1+=b)
for (j1 = 0; j1 < n; j1+=b)
for (k1 = 0; k1 < n; k1+=b)
for (i=i1; i<i1+b; i++)  {
  for (j=j1; j<j1+b; j++) {
    register double sum = C[i*n+j];
    for (k=k1; k<k1+b; k++)
      sum += A[i*n+k] * B[k*n+j];
    C[i*n+j] = sum;
  }

}
}
void jik(const double *A, const double *B, double *C, const int n)
{
int i,j,k;
for (j=0; j<n; j++) {
  for (i=0; i<n; i++) {
    register double sum = C[i*n+j];
    for (k=0; k<n; k++)
      sum += A[i*n+k] * B[k*n+j];
    C[i*n+j] = sum;
  }
}

}

void bjik(const double *A, const double *B, double *C, const int n, const int b)
{
int i1,j1,k1;
int i,j,k;
for (i1 = 0; i1 < n; i1+=b)
for (j1 = 0; j1 < n; j1+=b)
for (k1 = 0; k1 < n; k1+=b)
for (j=j1; j<j1+b; j++) {
  for (i=i1; i<i1+b; i++) {
    register double sum = C[i*n+j];
    for (k=k1; k<k1+b; k++)
      sum += A[i*n+k] * B[k*n+j];
    C[i*n+j] = sum;
  }
}

}

void kij(const double *A, const double *B, double *C, const int n)
{
int i,j,k;
for (k=0; k<n; k++) {
  for (i=0; i<n; i++) {
   register double r = A[i*n+k];
    for (j=0; j<n; j++)
      C[i*n+j] += r * B[k*n+j];
  }
}

}

void bkij(const double *A, const double *B, double *C, const int n, const int b)
{
int i1,j1,k1;
int i,j,k;
for (i1 = 0; i1 < n; i1+=b)
for (j1 = 0; j1 < n; j1+=b)
for (k1 = 0; k1 < n; k1+=b)
for (k=k1; k<k1+b; k++) {
  for (i=i1; i<i1+b; i++) {
   register double r = A[i*n+k];
    for (j=j1; j<j1+b; j++)
      C[i*n+j] += r * B[k*n+j];
  }
}

}


void ikj(const double *A, const double *B, double *C, const int n)
{
int i,j,k;
for (i=0; i<n; i++) {
  for (k=0; k<n; k++) {
    register double r = A[i*n+k];
    for (j=0; j<n; j++)
      C[i*n+j] += r * B[k*n+j];
  }
}

}

void bikj(const double *A, const double *B, double *C, const int n, const int b)
{
int i1,j1,k1;
int i,j,k;
for (i1 = 0; i1 < n; i1+=b)
for (j1 = 0; j1 < n; j1+=b)
for (k1 = 0; k1 < n; k1+=b)
for (i=i1; i<i1+b; i++) {
  for (k=k1; k<k1+b; k++) {
    register double r = A[i*n+k];
    for (j=j1; j<j1+b; j++)
      C[i*n+j] += r * B[k*n+j];
  }
}

}

void jki(const double *A, const double *B, double *C, const int n)
{
int i,j,k;
for (j=0; j<n; j++) {
  for (k=0; k<n; k++) {
   register double r = B[k*n+j];
    for (i=0; i<n; i++)
      C[i*n+j] += A[i*n+k] * r;
  }
}
}

void bjki(const double *A, const double *B, double *C, const int n, const int b)
{
int i1,j1,k1;
int i,j,k;
for (i1 = 0; i1 < n; i1+=b)
for (j1 = 0; j1 < n; j1+=b)
for (k1 = 0; k1 < n; k1+=b)
for (j=j1; j<j1+b; j++) {
  for (k=k1; k<k1+b; k++) {
   register double r = B[k*n+j];
    for (i=i1; i<i1+b; i++)
      C[i*n+j] += A[i*n+k] * r;
  }
}
}

void kji(const double *A, const double *B, double *C, const int n)
{
int i,j,k;
for (k=0; k<n; k++) {
  for (j=0; j<n; j++) {
    register double r = B[k*n+j];
    for (i=0; i<n; i++)
      C[i*n+j] += A[i*n+k] * r;
  }
}
}

void bkji(const double *A, const double *B, double *C, const int n, const int b)
{
int i1,j1,k1;
int i,j,k;
for (i1 = 0; i1 < n; i1+=b)
for (j1 = 0; j1 < n; j1+=b)
for (k1 = 0; k1 < n; k1+=b)
for (k=k1; k<k1+b; k++) {
  for (j=j1; j<j1+b; j++) {
    register double r = B[k*n+j];
    for (i=i1; i<i1+b; i++)
      C[i*n+j] += A[i*n+k] * r;
  }
}
}

//Cache Reuse part 3 End

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i, j, k;
    int i1,j1,k1;
    for (i = 0; i < n; i+=b)
        for (j = 0; j < n; j+=b)
             for (k = 0; k < n; k+=b)
                 /* B x B */
                for (i1=i; i1<i+b; i1++) {
                        for (k1=k; k1<k+b; k1++) {
                                register double r = A[i1*n+k1];
                                for (j1=j; j1<j+b; j1++)
                                        C[i1*n+j1] += r * B[k1*n+j1];
                        }
                }
}
//strassen
/*
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i, j, k;
    double m1,m2,m3,m4,m5,m6,m7;
    for (i = 0; i < n; i+=2)
        for (j = 0; j < n; j+=2)
             for (k = 0; k < n; k+=2)
                 {
                register a00=A[i*n+k];          register b00=B[k*n+j];
                register a01=A[i*n+k+1];        register b01=B[k*n+j+1];
                register a10=A[i*n+k+n];        register b10=B[k*n+j+n];
                register a11=A[i*n+k+n+1];      register b11=B[k*n+j+n+1];

                m1 = (a00 + a11) * (b00 + b11);
                m2 = (a10 + a11) * b00;
                m3 =  a00 * (b01 - b11);
                m4 =  a11 * (b10 - b00);
                m5 = (a00 + a01) * b11;
                m6 = (a10 - a00) * (b00 + b01);
                m7 = (a01 - a11) * (b10 + b11);

               register C[i*n+j]        += m1 + m4 - m5 + m7;
               register C[i*n+j+1]      += m3 + m5 ;
               register C[i*n+j+n]      += m2 + m4 ;
               register C[i*n+j+n+1]    += m1 + m3 - m2 + m6 ;
                }
}
*/
