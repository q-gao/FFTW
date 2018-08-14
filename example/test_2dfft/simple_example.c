/* cc -I$HOME/usr/include simple_example.c -L$HOME/usr/lib -lfftw3 -o simple_example */

#include <fftw3.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#ifdef _MSC_VER
#include <windows.h>
#define uint64_t UINT64
uint64_t get_time_ns() {
  return GetTickCount64() * 1000000UL; // ms -> ns
}

uint64_t get_time_ms() {
  return GetTickCount64();
}
#else // work in ANDROID and Linux(?)
#include <dlfcn.h>
uint64_t get_time_ns() {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return static_cast<uint64_t>(now.tv_sec) * 1000000000ULL + now.tv_nsec;
}

uint64_t get_time_ms() {
  struct timespec now;
  clock_gettime(CLOCK_MONOTONIC, &now);
  return static_cast<uint64_t>(now.tv_sec) * 1000ULL + now.tv_nsec/1000000ULL;
}
#endif

typedef double Float_T;
typedef struct {
	Float_T real, imag;
} Complex_T;

#define MAX(a,b) ((a)>(b)? (a) : (b))
#define ABS(v) ( (v) >=0? (v) : -(v))

int FFT(int dir, int m, Float_T *x, Float_T *y)
{
	long nn, i, i1, j, k, i2, l, l1, l2;
	Float_T c1, c2, tx, ty, t1, t2, u1, u2, z;

	/* Calculate the number of points */
	nn = 1;
	for (i = 0; i<m; i++)
		nn *= 2;

	/* Do the bit reversal */
	i2 = nn >> 1;
	j = 0;
	for (i = 0; i<nn - 1; i++) {
		if (i < j) {
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l = 0; l<m; l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j = 0; j<l1; j++) {
			for (i = j; i<nn; i += l2) {
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = x[i] - t1;
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
			}
			z = u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = (Float_T)sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = (Float_T)sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for reverse transform */
	if (dir == -1) {
		for (i = 0; i<nn; i++) {
			x[i] /= (Float_T)nn;
			y[i] /= (Float_T)nn;
		}
	}

	return 1;
}

/*-------------------------------------------------------------------------
Calculate the closest but lower power of two of a number
twopm = 2**m <= n
Return TRUE if 2**m == n
*/
inline int Powerof2(int n, int *m, int *twopm)
{
	if (n <= 1) {
		*m = 0;
		*twopm = 1;
		return 0;
	}

	*m = 1;
	*twopm = 2;
	do {
		(*m)++;
		(*twopm) *= 2;
	} while (2 * (*twopm) <= n);

	if (*twopm != n)
		return 0;
	else
		return 1;
}

/*-------------------------------------------------------------------------
Perform a 2D FFT inplace given a complex 2D array
The direction dir, 1 for forward, -1 for reverse
The size of the array (nx,ny)
Return false if there are memory problems or
the dimensions are not powers of 2
*/
int FFT2D(Complex_T **c, int nx, int ny, int dir, Float_T *real, Float_T * imag)
/*
ARGUMENTS:
- c: input data in Complex
- real, imag: work buff of one row/column length for 1D FFT
*/
{
	int i, j;
	int m, twopm;

	/* Transform the rows */
	Powerof2(nx, &m, &twopm);
	for (j = 0; j<ny; j++) {
		for (i = 0; i<nx; i++) {
			real[i] = c[i][j].real;
			imag[i] = c[i][j].imag;
		}
		FFT(dir, m, real, imag);
		for (i = 0; i<nx; i++) {
			c[i][j].real = real[i];
			c[i][j].imag = imag[i];
		}
	}

	for (i = 0; i<nx; i++) {
		for (j = 0; j<ny; j++) {
			real[j] = c[i][j].real;
			imag[j] = c[i][j].imag;
		}
		FFT(dir, m, real, imag);
		for (j = 0; j<ny; j++) {
			c[i][j].real = real[j];
			c[i][j].imag = imag[j];
		}
	}

	return 1;
}

int FFT2D_1DInput(Float_T * data_real, Float_T* data_imag, int nx, int ny, int dir, Float_T *real, Float_T * imag)
/*
ARGUMENTS:
- c: input data in Complex
- real, imag: work buff of one row/column length for 1D FFT
*/
{
	int i, j;
	int m, twopm;

	/* Transform the rows */
	Powerof2(nx, &m, &twopm);

	Float_T * pr, *pi;
	pr = data_real; pi = data_imag;
	for (j = 0; j < ny; j++) {
		FFT(dir, m, pr, pi);
		pr += nx;
		pi += nx;
	}

	for (i = 0; i<nx; i++) {
		pr = data_real + i; pi = data_imag + i;
		for (j = 0; j<ny; j++) {
			real[j] = *pr;
			imag[j] = *pi;
			pr += nx;
			pi += nx;
		}
		FFT(dir, m, real, imag);
		pr = data_real + i; pi = data_imag + i;
		for (j = 0; j<ny; j++) {
			*pr = real[j];
			*pi = imag[j];
			pr += nx;
			pi += nx;
		}
	}

	return 1;
}

// TODO: assume N0 == N1
#define N0  16
#define N1  16
#define VERIFY_RESULT
int main(int argc, char **argv){
  //const ptrdiff_t N0 = 16, N1 = 16;  
  fftw_plan plan, plan_realin;
  fftw_complex *fftwCplxData=(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N0 * N1);
  double       *fftwRealData = (double*)fftw_malloc(sizeof(double) * N0 * N1);
#ifdef VERIFY_RESULT
  Float_T real[MAX(N0, N1)], imag[MAX(N0, N1)];
  Float_T data_real[N0 * N1], data_imag[N0 * N1];
  Float_T * pr, *pi;
#endif

  //printf("fftw_plan = %lu bytes, fftw_complex = %lu bytes\n", sizeof(fftw_plan), sizeof(fftw_complex));
  uint64_t start_ns = get_time_ns();
  /* create plan for forward DFT */
#ifndef VERIFY_RESULT
  plan        = fftw_plan_dft_2d    (N0, N1, fftwCplxData, fftwCplxData,     FFTW_FORWARD, FFTW_ESTIMATE);
#endif
  plan_realin = fftw_plan_dft_r2c_2d(N0, N1, fftwRealData, fftwCplxData,FFTW_ESTIMATE);
  printf("FFTW init time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns)));

  /* initialize fftwCplxData to some function my_function(x,y) */
  int i, j;
  double pdata=0;
#ifdef VERIFY_RESULT
  pr = data_real;
  pi = data_imag;
#endif
  for (i = 0; i < N0; ++i){
    for (j = 0; j < N1; ++j){
      fftwCplxData[i*N1 + j][0]=i; 
      fftwCplxData[i*N1 + j][1]=0;
      fftwRealData[i*N1 + j] = i;
      //pdata+=fftwCplxData[i*N1 + j][0]*fftwCplxData[i*N1 + j][0]+fftwCplxData[i*N1 + j][1]*fftwCplxData[i*N1 + j][1];
#ifdef VERIFY_RESULT
	  *pr++ = i;
	  *pi++ = 0;
#endif
    }
  }
  //printf("power of original fftwCplxData is %f\n", pdata);

#ifdef VERIFY_RESULT
  fftw_execute(plan_realin);
  FFT2D_1DInput(data_real, data_imag, N0, N1, 1, real, imag);

  Float_T d, max_diff = (Float_T)0.0;
  int max_diff_i=-1, max_diff_j=-1;
  pr = data_real;
  pi = data_imag;
  for (i = 0; i < N0; ++i) {
	  for (j = 0; j < N1; ++j) {
		  printf("(%d,%d): (%lf, %lf) vs. (%lf, %lf) \n", i, j, (double)*pr, (double)*pi,
				(double)fftwCplxData[i*N1 + j][0],
				(double)fftwCplxData[i*N1 + j][1]
		  );
		  d = ABS(fftwCplxData[i*N1 + j][0] - *pr);
		  if (max_diff < d) {
			  max_diff = d;
			  max_diff_i = i; max_diff_j = j;
		  }
		  d = ABS(fftwCplxData[i*N1 + j][1] - *pi);
		  if (max_diff < d) {
			  max_diff = d;
			  max_diff_i = i; max_diff_j = j;
		  }
		  pr++; pi++;
		  //pdata+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
	  }
  }
  printf("Max Diff %f at (%d, %d)\n", max_diff, max_diff_i, max_diff_j);
#else
  const int num_run = 2000000;

  printf("Complext 2D FFT:\n");
  start_ns = get_time_ns();
  for(int i = 0; i < num_run; i++)
    /* compute transforms, in-place, as many times as desired */
    fftw_execute(plan);
  printf("  Total 2D FFT Time %llu us\n", (unsigned long long)((get_time_ns() - start_ns)/1000));
  printf("  Avg 2D FFT Time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns)/num_run));

  printf("Real 2D FFT:\n");
  start_ns = get_time_ns();
  for(int i = 0; i < num_run; i++)
    /* compute transforms, in-place, as many times as desired */
    fftw_execute(plan_realin);
  printf("  Total 2D FFT Time %llu us\n", (unsigned long long)((get_time_ns() - start_ns)/1000));
  printf("  Avg 2D FFT Time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns)/num_run));
#endif
  // double normalization=sqrt((double)N0*N1);
  // double ptransform = 0;
  // for (i = 0; i < N0; ++i){
  //   for (j = 0; j < N1; ++j){
  //     fftwCplxData[i*N1+j][0]/=normalization;
  //     fftwCplxData[i*N1+j][1]/=normalization;
  //     ptransform+=fftwCplxData[i*N1 + j][0]*fftwCplxData[i*N1 + j][0]+fftwCplxData[i*N1 + j][1]*fftwCplxData[i*N1 + j][1];
  //   }
  // }

  //printf("power of transform is %f\n", ptransform);
 
#ifndef VERIFY_RESULT
  fftw_destroy_plan(plan);
#endif
  fftw_destroy_plan(plan_realin);

  return 0;
}

