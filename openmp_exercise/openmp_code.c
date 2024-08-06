#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>


void matrix_multi(int m, int n, int k, double *a, double *b, double *c){
     for(int i = 0; i < m; ++i){
       for(int j = 0; j < n; ++j){
          for(int l = 0; l < k; ++l){
              c[i * n + j] += a[i * k + l] * b[l * n + j];
         }
       }
     }
  }

void test_parallel_code(int m, int n, int k, double *a, double *b, double *c, int threads){
     for(int i = 0; i < m; ++i){
       #pragma omp parallel for num_threads(threads)	     
       for(int j = 0; j < n; ++j){
          for(int l = 0; l < k; ++l){
              c[i * n + j]+= a[i * k + l] * b[l * n + j];
         }
       }
     }
  }
int main(int argc, char **argv){
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int k = atoi(argv[3]);
    int t = atoi(argv[4]);
    int w = atoi(argv[5]);
    int h = atoi(argv[6]);

    double *a, *b, *c0, *c1;
    a = (double *)malloc(m * k * sizeof(double));
    b = (double *)malloc(k * n * sizeof(double));
    c0 = (double *)malloc(m * n * sizeof(double));
    c1 = (double *)malloc(m * n * sizeof(double));

    for (int i = 0; i < m; ++i) {
        for (int l = 0; l < k; ++l) {
            a[i * k + l] = (i + 1);
        }
    }

    for (int l = 0; l < k; ++l) {
        for (int j = 0; j < n; ++j) {
            b[l * n + j] = (l + 2);
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            c0[i * n + j] = 0;
            c1[i * n + j] = 0;
        }
    }

    long long sum0 = 0, sum1=0;
    struct timeval t0, t1;

    for(int r = 0; r < (w + h); ++r) {	    
       gettimeofday(&t0, NULL);
       matrix_multi(m, n, k, a, b, c0);
       gettimeofday(&t1, NULL);

    if (r >= w)
      sum0 += (t1.tv_sec - t0.tv_sec) * 1000000 + (t1.tv_usec - t0.tv_usec);
    }

    for(int r = 0; r < (w + h); ++r) {	    
       gettimeofday(&t0, NULL);
       test_parallel_code(m, n, k, a, b, c1, t);
       gettimeofday(&t1, NULL);

    if (r >= w)
      sum1 += (t1.tv_sec - t0.tv_sec) * 1000000 + (t1.tv_usec - t0.tv_usec);
    }    

    int correct = 1;
    for(int i = 0; i < m; ++i){
       for(int j = 0; j< n; ++j){      
         if(fabs(c0[i * n + j] - c1[i * n + j]) > 1e-6) correct = 0;
       }
     }
    printf("%d\t%lf\t%lf\n", correct, sum0/((double) (h * 1.0)), sum1/((double) (h * 1.0)));

    free(a);
    free(b);
    free(c0);
    free(c1);

   return 0;
}

