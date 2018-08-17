#include <stdio.h>
#include <stdint.h>
#include <math.h>

typedef double FLOAT_T;
typedef struct {
	FLOAT_T real, imag;
} Complex_T;

#define ABS(v) ( (v) >=0? (v) : -(v))

int FFT(int dir, int m, FLOAT_T *x, FLOAT_T *y)
{
	long nn, i, i1, j, k, i2, l, l1, l2;
	FLOAT_T c1, c2, tx, ty, t1, t2, u1, u2, z;

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
		c2 = (FLOAT_T)sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = (FLOAT_T) sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for forward transform */
	if (dir == -1) {
		for (i = 0; i<nn; i++) {
			x[i] /= (FLOAT_T)nn;
			y[i] /= (FLOAT_T)nn;
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
int FFT2D(Complex_T **c, int nx, int ny, int dir, FLOAT_T *real, FLOAT_T * imag)
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
#if 0
int FFT2D_1DInput(FLOAT_T * data_real, FLOAT_T* data_imag, int nx, int ny, int dir, FLOAT_T *real, FLOAT_T * imag)
/*
ARGUMENTS:
- data_real, data_imag: input data
- real, imag: work buff of one row/column length for 1D FFT
*/
{
	int i, j;
	int m, twopm;

	/* Transform the rows */
	Powerof2(nx, &m, &twopm);

	FLOAT_T * pr, *pi;
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
#endif
int FFT2D_1DInput(Complex_T * data, int nx, int ny, int dir, FLOAT_T *real, FLOAT_T * imag)
/*
ARGUMENTS:
- c: input data in Complex
- real, imag: work buff of one row/column length for 1D FFT
*/
{
	int i, j;
	int m, twopm;
	Complex_T * prow, *pcol;

	/* Transform the rows */
	Powerof2(nx, &m, &twopm);

	prow = data;
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			real[j] = prow[j].real;
			imag[j] = prow[j].imag;
		}
		FFT(dir, m, real, imag);

		for (j = 0; j < ny; j++) {
			prow[j].real = real[j];
			prow[j].imag = imag[j];
		}
		prow += nx;
	}

	pcol = data;
	for (i = 0; i < nx; i++) {
		prow = pcol;
		for (j = 0; j < ny; j++) {
			real[j] = prow->real;
			imag[j] = prow->imag;
			prow += nx;
		}
		FFT(dir, m, real, imag);

		prow = pcol;
		for (j = 0; j < ny; j++) {
			prow->real = real[j];
			prow->imag = imag[j];
			prow += nx;
		}

		pcol++;
	}

	return 1;
}

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
#include <time.h>
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

     
int main(int argc, char **argv)
{
	const int N = 16;
	uint64_t start_ns;
	int i, j;
	FLOAT_T real[N], imag[N];
	FLOAT_T data_real[N*N], data_imag[N*N];
	Complex_T ** data;
	data = new Complex_T*[N];
	for (int i = 0; i < N; i++)
		data[i] = new Complex_T[N];
	Complex_T data_c[N*N];

	FLOAT_T pdata=0;
	FLOAT_T * pr, *pi;
	pr = data_real;
	pi = data_imag;
	for (i = 0; i < N; ++i){
		for (j = 0; j < N; ++j){
		  data[i][j].real = rand(); 
		  data[i][j].imag = rand();

		  *pr++ = i;
		  *pi++ = 0;

		  data_c[i * N + j].real = data[i][j].real;
		  data_c[i * N + j].imag = data[i][j].imag;
		  //pdata+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
		}
	}
	//printf("power of original data is %f\n", pdata);

	FFT2D(data, N, N, 1, real, imag);
	FFT2D_1DInput(data_c, N, N, 1, real, imag);
	double max_diff = -10.0, d;
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			d = ABS( data[i][j].real - data_c[i*N+j].real);
			if (max_diff < d)
				max_diff = d;

			d = ABS(data[i][j].imag - data_c[i*N + j].imag);
			if (max_diff < d)
				max_diff = d;
		}
	}

	printf("Max Diff %lf\n", max_diff);
	//const int num_run = 2000000;
	//start_ns = get_time_ns();
	//for(int i = 0; i < num_run; i++)
	//	FFT2D(data, N, N, 1, real, imag);
	//printf("Total 2D FFT Time %llu us\n", (unsigned long long)((get_time_ns() - start_ns)/1000));
	//printf("Avg 2D FFT Time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns)/num_run));

	//printf("\n");
	//start_ns = get_time_ns();
	//for (int i = 0; i < num_run; i++)
	//	FFT2D_1DInput(data_c, N, N, 1, real, imag);
	//printf("Total 2D FFT Time %llu us\n", (unsigned long long)((get_time_ns() - start_ns) / 1000));
	//printf("Avg 2D FFT Time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns) / num_run));

	//for (int i = 0; i < num_run; i++)
	//	FFT2D_1DInput(data_real, data_imag, N, N, 1, real, imag);
	//printf("Total 2D FFT Time %llu us\n", (unsigned long long)((get_time_ns() - start_ns) / 1000));
	//printf("Avg 2D FFT Time %llu ns\n", (unsigned long long)((get_time_ns() - start_ns) / num_run));


	// FLOAT_T normalization=sqrt((FLOAT_T)N0*N1);
	// FLOAT_T ptransform = 0;
	// for (i = 0; i < N0; ++i){
	//   for (j = 0; j < N1; ++j){
	//     data[i*N1+j][0]/=normalization;
	//     data[i*N1+j][1]/=normalization;
	//     ptransform+=data[i*N1 + j][0]*data[i*N1 + j][0]+data[i*N1 + j][1]*data[i*N1 + j][1];
	//   }
	// }

	//printf("power of transform is %f\n", ptransform);
#if 0
	FFT2D(data, N, N, 1, real, imag);
	FFT2D_1DInput(data_real, data_imag, N, N, 1, real, imag);


#endif
	// int v = -10;
	// int bs = 1;
	// printf("Arithmetic shift test: %d << %d = %d\n", v, bs, v << bs);
	// printf("Arithmetic shift test: %d >> %d = %d\n", v, bs, v >> bs);

	return 0;
}