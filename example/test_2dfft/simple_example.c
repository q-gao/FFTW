/* cc -I$HOME/usr/include simple_example.c -L$HOME/usr/lib -lfftw3 -o simple_example */

#include <fftw3.h>
#include <math.h>
#include <stdint.h>
#include <dlfcn.h>
#include <time.h>
uint64_t get_time_ns() {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return (uint64_t)(now.tv_sec * 1000000000ULL + now.tv_nsec);
}
     
int main(int argc, char **argv){
  const ptrdiff_t N0 = 16, N1 = 16;
  fftw_plan plan;
  fftw_complex *data;

  printf("fftw_plan = %lu bytes, fftw_complex = %lu bytes\n", sizeof(fftw_plan), sizeof(fftw_complex));
  uint64_t start_ns = get_time_ns();
  data = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 * N1);

  /* create plan for forward DFT */
  plan = fftw_plan_dft_2d(N0, N1, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
  
  printf("FFTW init time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns)));
  /* initialize data to some function my_function(x,y) */
  int i, j;
  double pdata=0;
  for (i = 0; i < N0; ++i){
    for (j = 0; j < N1; ++j){
      data[i*N1 + j][0]=i; 
      data[i*N1 + j][1]=0;
      //pdata+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
    }
  }
  //printf("power of original data is %f\n", pdata);

  const int num_run = 2000000;
  start_ns = get_time_ns();
  for(int i = 0; i < num_run; i++)
    /* compute transforms, in-place, as many times as desired */
    fftw_execute(plan);

  printf("Total 2D FFT Time %llu us\n", (unsigned long long)((get_time_ns() - start_ns)/1000));
  printf("Avg 2D FFT Time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns)/num_run));
  // double normalization=sqrt((double)N0*N1);
  // double ptransform = 0;
  // for (i = 0; i < N0; ++i){
  //   for (j = 0; j < N1; ++j){
  //     data[i*N1+j][0]/=normalization;
  //     data[i*N1+j][1]/=normalization;
  //     ptransform+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
  //   }
  // }

  //printf("power of transform is %f\n", ptransform);
 
  fftw_destroy_plan(plan);

  return 0;
}

