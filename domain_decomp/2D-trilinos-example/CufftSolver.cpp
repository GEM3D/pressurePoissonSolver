#include "CufftSolver.h"
#include <fftw3.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
using namespace std;
// intialize static members
map<CDomainKey, valarray<double>> CufftSolver::denoms;

__global__ void scale(cufftDoubleReal *in, double scale, int n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		in[i] *= scale;
	}
}
__global__ void dct2_pre(cufftDoubleReal *in, cufftDoubleComplex *out, cufftDoubleComplex *w, int n,
                         int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (row == 0) {
			double x = col * M_PI / (2 * n);
			w[col].x = 2 * cos(x);
			w[col].y = 2 * sin(x);
		}
		if (dim == 0) {
			// rows
			out[2 * n * row + col].x               = in[i];
			out[2 * n * row + (2 * n - 1 - col)].x = in[i];
		} else {
			// cols
			out[n * row + col].x               = in[i];
			out[n * (2 * n - 1 - row) + col].x = in[i];
		}
	}
}

__global__ void dct2_post(cufftDoubleComplex *in, cufftDoubleReal *out, cufftDoubleComplex *w,
                          int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			// rows
			double x      = in[2 * n * row + col].x;
			double y      = in[2 * n * row + col].y;
			double w_x    = w[col].x;
			double w_y    = w[col].y;
			double result = x * w_x - x * w_y - y * w_x - y * w_y;
			result *= n / 2.0;
			out[n * row + col] = result;
		} else {
			// cols
			double x      = in[n * row + col].x;
			double y      = in[n * row + col].y;
			double w_x    = w[row].x;
			double w_y    = w[row].y;
			double result = x * w_x - x * w_y - y * w_x - y * w_y;
			result *= n / 2.0;
			out[n * row + col] = result;
		}
	}
}
cufftDoubleReal *dct2(cufftDoubleReal *x, int n, int dim)
{
	// allocate fft input and memset
	cufftDoubleComplex *in = nullptr;
	cudaMalloc((void **) &in, sizeof(cufftDoubleComplex) * n * n * 2);
	cudaMemset((void *) in, 0, sizeof(cufftDoubleComplex) * n * n * 2);

	// allocate weight vector
	cufftDoubleComplex *w = nullptr;
	cudaMalloc((void **) &w, sizeof(cufftDoubleComplex) * n);

	// prepare fft input
	dct2_pre<<<n, n>>>(x, in, w, n, dim);

	// create plan
	cufftHandle plan;
	if (dim == 0) {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = 1;
		int idist          = 2 * n;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	} else {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = n;
		int idist          = 1;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	}

	// allocate fft output
	cufftDoubleComplex *out = nullptr;
	cudaMalloc((void **) &out, sizeof(cufftDoubleComplex) * n * n * 2);

	// perform transform
	int i = cufftExecZ2Z(plan, in, out, CUFFT_INVERSE);
	cufftDestroy(plan);

	// free input
	cudaFree(in);

	// allocate final output
	cufftDoubleReal *retval = nullptr;
	cudaMalloc((void **) &retval, sizeof(cufftDoubleReal) * n * n);

	// get final output
	dct2_post<<<n, n>>>(out, retval, w, n, dim);

	// free weights and fft output
	cudaFree(w);
	cudaFree(out);

	return retval;
}
__global__ void idct2_pre(cufftDoubleReal *in, cufftDoubleComplex *out, cufftDoubleComplex *w,

                          int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (row == 0) {
			double x = col * M_PI / (2 * n);
			w[col].x = cos(-x) / 2;
			w[col].y = sin(-x) / 2;
			if (col == 0) {
				w[col].x *= 2;
				w[col].y *= 2;
			}
			x            = (n + col) * M_PI / (2 * n);
			w[n + col].x = -cos(-x) / 2;
			w[n + col].y = -sin(-x) / 2;
			if (col == 0) {
				w[n].x = 0;
				w[n].y = 0;
			}
		}
		if (dim == 0) {
			// rows
			if (col == 0) {
				out[2 * n * row].x     = in[i] / 2.0;
				out[2 * n * row + n].x = 1;
			} else {
				out[2 * n * row + col].x           = in[i];
				out[2 * n * row + (2 * n - col)].x = in[i];
			}
		} else {
			// cols
			if (row == 0) {
				out[col].x           = in[i] / 2.0;
				out[n * (n) + col].x = 1;
			} else {
				out[n * row + col].x           = in[i];
				out[n * (2 * n - row) + col].x = in[i];
			}
		}
	}
}

__global__ void idct2_weight(cufftDoubleComplex *in, cufftDoubleComplex *w, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			double x, y, w_x, w_y;
			x                           = in[2 * n * row + col].x;
			y                           = in[2 * n * row + col].y;
			w_x                         = w[col].x;
			w_y                         = w[col].y;
			in[2 * n * row + col].x     = x * w_x - y * w_y;
			in[2 * n * row + col].y     = x * w_y + y * w_x;
			x                           = in[2 * n * row + n + col].x;
			y                           = in[2 * n * row + n + col].y;
			w_x                         = w[col + n].x;
			w_y                         = w[col + n].y;
			in[2 * n * row + col + n].x = x * w_x - y * w_y;
			in[2 * n * row + col + n].y = x * w_y + y * w_x;
		} else {
			double x, y, w_x, w_y;
			x                         = in[n * row + col].x;
			y                         = in[n * row + col].y;
			w_x                       = w[row].x;
			w_y                       = w[row].y;
			in[n * row + col].x       = x * w_x - y * w_y;
			in[n * row + col].y       = x * w_y + y * w_x;
			x                         = in[n * (n + row) + col].x;
			y                         = in[n * (n + row) + col].y;
			w_x                       = w[n + row].x;
			w_y                       = w[n + row].y;
			in[n * (n + row) + col].x = x * w_x - y * w_y;
			in[n * (n + row) + col].y = x * w_y + y * w_x;
		}
	}
}
__global__ void idct2_post(cufftDoubleComplex *in, cufftDoubleReal *out, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			// rows
			double x      = in[2 * n * row + col].x;
			double y      = in[2 * n * row + col].y;
			double result = x - y;
			result *= 2.0 / n;
			out[n * row + col] = result;
		} else {
			// cols
			double x      = in[n * row + col].x;
			double y      = in[n * row + col].y;
			double result = x - y;
			result *= 2.0 / n;
			out[n * row + col] = result;
		}
	}
}

cufftDoubleReal *idct2(cufftDoubleReal *x, int n, int dim)
{
	// allocate fft input and memset
	cufftDoubleComplex *in = nullptr;
	cudaMalloc((void **) &in, sizeof(cufftDoubleComplex) * n * n * 2);
	cudaMemset((void *) in, 0, sizeof(cufftDoubleComplex) * n * n * 2);

	// allocate weight vector
	cufftDoubleComplex *w = nullptr;
	cudaMalloc((void **) &w, 2 * sizeof(cufftDoubleComplex) * n);

	// prepare fft input
	idct2_pre<<<n, n>>>(x, in, w, n, dim);
	idct2_weight<<<n, n>>>(in, w, n, dim);

	// create plan
	cufftHandle plan;
	if (dim == 0) {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = 1;
		int idist          = 2 * n;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	} else {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = n;
		int idist          = 1;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	}

	// allocate fft output
	cufftDoubleComplex *out = nullptr;
	cudaMalloc((void **) &out, sizeof(cufftDoubleComplex) * n * n * 2);

	// perform transform
	int i = cufftExecZ2Z(plan, in, out, CUFFT_FORWARD);
	cufftDestroy(plan);

	// free input
	cudaFree(in);

	// allocate final output
	cufftDoubleReal *retval = nullptr;
	cudaMalloc((void **) &retval, sizeof(cufftDoubleReal) * n * n);

	// get final output
	idct2_post<<<n, n>>>(out, retval, n, dim);
	scale<<<n, n>>>(retval, 1.0/(2.0*n), n);

	// free weights and fft output
	cudaFree(w);
	cudaFree(out);

	return retval;
}
__global__ void dst2_pre(cufftDoubleReal *in, cufftDoubleComplex *out, cufftDoubleComplex *w, int n,
                         int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (row == 0) {
			double x = col * M_PI / (2 * n);
			w[col].x = 2 * cos(x);
			w[col].y = 2 * sin(x);
		}
		if (dim == 0) {
			// rows
			if (col % 2 == 0) {
				out[2 * n * row + (n - col - 1)].x = -in[i];
				out[2 * n * row + (n + col)].x     = -in[i];
			} else {
				out[2 * n * row + (n - col - 1)].x = in[i];
				out[2 * n * row + (n + col)].x     = in[i];
			}
		} else {
			// cols
			if (row % 2 == 0) {
				out[n * (n - row - 1) + col].x = -in[i];
				out[n * (n + row) + col].x     = -in[i];
			} else {
				out[n * (n - row - 1) + col].x = in[i];
				out[n * (n + row) + col].x     = in[i];
			}
		}
	}
}

__global__ void dst2_post(cufftDoubleComplex *in, cufftDoubleReal *out, cufftDoubleComplex *w,
                          int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			// rows
			double x      = in[2 * n * row + col].x;
			double y      = in[2 * n * row + col].y;
			double w_x    = w[col].x;
			double w_y    = w[col].y;
			double result = x * w_x - x * w_y - y * w_x - y * w_y;
			result *= n / 2.0;
			if (col % 2 == 0) {
				out[n * row + (n - col - 1)] = -result;
			} else {
				out[n * row + (n - col - 1)] = result;
			}
		} else {
			// cols
			double x      = in[n * row + col].x;
			double y      = in[n * row + col].y;
			double w_x    = w[row].x;
			double w_y    = w[row].y;
			double result = x * w_x - x * w_y - y * w_x - y * w_y;
			result *= n / 2.0;
			if (row % 2 == 0) {
				out[n * (n - row - 1) + col] = -result;
			} else {
				out[n * (n - row - 1) + col] = result;
			}
		}
	}
}
cufftDoubleReal *dst2(cufftDoubleReal *x, int n, int dim)
{
	// allocate fft input and memset
	cufftDoubleComplex *in = nullptr;
	cudaMalloc((void **) &in, sizeof(cufftDoubleComplex) * n * n * 2);
	cudaMemset((void *) in, 0, sizeof(cufftDoubleComplex) * n * n * 2);

	// allocate weight vector
	cufftDoubleComplex *w = nullptr;
	cudaMalloc((void **) &w, sizeof(cufftDoubleComplex) * n);

	// prepare fft input
	dst2_pre<<<n, n>>>(x, in, w, n, dim);

	// create plan
	cufftHandle plan;
	if (dim == 0) {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = 1;
		int idist          = 2 * n;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	} else {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = n;
		int idist          = 1;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	}

	// allocate fft output
	cufftDoubleComplex *out = nullptr;
	cudaMalloc((void **) &out, sizeof(cufftDoubleComplex) * n * n * 2);

	// perform transform
	int i = cufftExecZ2Z(plan, in, out, CUFFT_INVERSE);
	cufftDestroy(plan);

	// free input
	cudaFree(in);

	// allocate final output
	cufftDoubleReal *retval = nullptr;
	cudaMalloc((void **) &retval, sizeof(cufftDoubleReal) * n * n);

	// get final output
	dst2_post<<<n, n>>>(out, retval, w, n, dim);

	// free weights and fft output
	cudaFree(w);
	cudaFree(out);

	return retval;
}
__global__ void idst2_pre(cufftDoubleReal *in, cufftDoubleComplex *out, cufftDoubleComplex *w,

                          int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (row == 0) {
			double x = col * M_PI / (2 * n);
			w[col].x = cos(-x) / 2;
			w[col].y = sin(-x) / 2;
			if (col == 0) {
				w[col].x *= 2;
				w[col].y *= 2;
			}
			x            = (n + col) * M_PI / (2 * n);
			w[n + col].x = -cos(-x) / 2;
			w[n + col].y = -sin(-x) / 2;
			if (col == 0) {
				w[n].x = 0;
				w[n].y = 0;
			}
		}
		if (dim == 0) {
			// rows
			if (n - col - 1 == 0) {
				if (col % 2 == 0) {
					out[2 * n * row].x = -in[i] / 2.0;
				} else {
					out[2 * n * row].x = in[i] / 2.0;
				}
				out[2 * n * row + n].x = 1;
			} else {
				if (col % 2 == 0) {
					out[2 * n * row + (n - col - 1)].x = -in[i];
					out[2 * n * row + (n + col + 1)].x = -in[i];
				} else {
					out[2 * n * row + (n - col - 1)].x = in[i];
					out[2 * n * row + (n + col + 1)].x = in[i];
				}
			}
		} else {
			// cols
			if (n - row - 1 == 0) {
				if (row % 2 == 0) {
					out[col].x = -in[i] / 2.0;
				} else {
					out[col].x = in[i] / 2.0;
				}
				out[n * (n) + col].x = 1;
			} else {
				if (row % 2 == 0) {
					out[n * (n - row - 1) + col].x = -in[i];
					out[n * (n + row + 1) + col].x = -in[i];
				} else {
					out[n * (n - row - 1) + col].x = in[i];
					out[n * (n + row + 1) + col].x = in[i];
				}
			}
		}
	}
}

__global__ void idst2_weight(cufftDoubleComplex *in, cufftDoubleComplex *w, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			double x, y, w_x, w_y;
			x                           = in[2 * n * row + col].x;
			y                           = in[2 * n * row + col].y;
			w_x                         = w[col].x;
			w_y                         = w[col].y;
			in[2 * n * row + col].x     = x * w_x - y * w_y;
			in[2 * n * row + col].y     = x * w_y + y * w_x;
			x                           = in[2 * n * row + n + col].x;
			y                           = in[2 * n * row + n + col].y;
			w_x                         = w[col + n].x;
			w_y                         = w[col + n].y;
			in[2 * n * row + col + n].x = x * w_x - y * w_y;
			in[2 * n * row + col + n].y = x * w_y + y * w_x;
		} else {
			double x, y, w_x, w_y;
			x                         = in[n * row + col].x;
			y                         = in[n * row + col].y;
			w_x                       = w[row].x;
			w_y                       = w[row].y;
			in[n * row + col].x       = x * w_x - y * w_y;
			in[n * row + col].y       = x * w_y + y * w_x;
			x                         = in[n * (n + row) + col].x;
			y                         = in[n * (n + row) + col].y;
			w_x                       = w[n + row].x;
			w_y                       = w[n + row].y;
			in[n * (n + row) + col].x = x * w_x - y * w_y;
			in[n * (n + row) + col].y = x * w_y + y * w_x;
		}
	}
}
__global__ void idst2_post(cufftDoubleComplex *in, cufftDoubleReal *out, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			// rows
			double x      = in[2 * n * row + col].x;
			double y      = in[2 * n * row + col].y;
			double result = x - y;
			result *= 2.0 / n;
			if (col % 2 == 0) {
				out[n * row + (n - col - 1)] = -result;
			} else {
				out[n * row + (n - col - 1)] = result;
			}
		} else {
			// cols
			double x      = in[n * row + col].x;
			double y      = in[n * row + col].y;
			double result = x - y;
			result *= 2.0 / n;
			if (row % 2 == 0) {
				out[n * (n - row - 1) + col] = -result;
			} else {
				out[n * (n - row - 1) + col] = result;
			}
		}
	}
}

cufftDoubleReal *idst2(cufftDoubleReal *x, int n, int dim)
{
	// allocate fft input and memset
	cufftDoubleComplex *in = nullptr;
	cudaMalloc((void **) &in, sizeof(cufftDoubleComplex) * n * n * 2);
	cudaMemset((void *) in, 0, sizeof(cufftDoubleComplex) * n * n * 2);

	// allocate weight vector
	cufftDoubleComplex *w = nullptr;
	cudaMalloc((void **) &w, 2 * sizeof(cufftDoubleComplex) * n);

	// prepare fft input
	idst2_pre<<<n, n>>>(x, in, w, n, dim);
	idst2_weight<<<n, n>>>(in, w, n, dim);

	// create plan
	cufftHandle plan;
	if (dim == 0) {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = 1;
		int idist          = 2 * n;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	} else {
		int transform_size = 2 * n;
		int inembed        = 2 * n;
		int istride        = n;
		int idist          = 1;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	}

	// allocate fft output
	cufftDoubleComplex *out = nullptr;
	cudaMalloc((void **) &out, sizeof(cufftDoubleComplex) * n * n * 2);

	// perform transform
	int i = cufftExecZ2Z(plan, in, out, CUFFT_FORWARD);
	cufftDestroy(plan);

	// free input
	cudaFree(in);

	// allocate final output
	cufftDoubleReal *retval = nullptr;
	cudaMalloc((void **) &retval, sizeof(cufftDoubleReal) * n * n);

	// get final output
	idst2_post<<<n, n>>>(out, retval, n, dim);
	scale<<<n, n>>>(retval, 1.0/(2.0*n), n);

	// free weights and fft output
	cudaFree(w);
	cudaFree(out);

	return retval;
}
__global__ void dct4_pre(cufftDoubleReal *in, cufftDoubleComplex *out, cufftDoubleComplex *w,

                         int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (row == 0) {
			double x = col * M_PI / (4 * n);
			w[col].x = cos(-x) / 2;
			w[col].y = sin(-x) / 2;
			if (col == 0) {
				w[0].x *= 2;
				w[0].y *= 2;
			}
			x                = (n + col) * M_PI / (4 * n);
			w[n + col].x     = cos(-x) / 2;
			w[n + col].y     = sin(-x) / 2;
			x                = (2 * n + col) * M_PI / (4 * n);
			w[2 * n + col].x = -cos(-x) / 2;
			w[2 * n + col].y = -sin(-x) / 2;
			if (col == 0) {
				w[2 * n].x = 0;
				w[2 * n].y = 0;
			}
			x                = (3 * n + col) * M_PI / (4 * n);
			w[3 * n + col].x = -cos(-x) / 2;
			w[3 * n + col].y = -sin(-x) / 2;
		}
		if (dim == 0) {
			// rows
			out[4 * n * row + 2 * col + 1].x         = in[i];
			out[4 * n * row + 4 * n - 2 * col - 1].x = in[i];
		} else {
			// cols
			out[n * (2 * row + 1) + col].x         = in[i];
			out[n * (4 * n - 2 * row - 1) + col].x = in[i];
		}
	}
}

__global__ void dct4_weight(cufftDoubleComplex *in, cufftDoubleComplex *w, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			double x, y, w_x, w_y;
			x                               = in[4 * n * row + col].x;
			y                               = in[4 * n * row + col].y;
			w_x                             = w[col].x;
			w_y                             = w[col].y;
			in[4 * n * row + col].x         = x * w_x - y * w_y;
			in[4 * n * row + col].y         = x * w_y + y * w_x;
			x                               = in[4 * n * row + n + col].x;
			y                               = in[4 * n * row + n + col].y;
			w_x                             = w[col + n].x;
			w_y                             = w[col + n].y;
			in[4 * n * row + col + n].x     = x * w_x - y * w_y;
			in[4 * n * row + col + n].y     = x * w_y + y * w_x;
			x                               = in[4 * n * row + 2 * n + col].x;
			y                               = in[4 * n * row + 2 * n + col].y;
			w_x                             = w[col + 2 * n].x;
			w_y                             = w[col + 2 * n].y;
			in[4 * n * row + col + 2 * n].x = x * w_x - y * w_y;
			in[4 * n * row + col + 2 * n].y = x * w_y + y * w_x;
			x                               = in[4 * n * row + 3 * n + col].x;
			y                               = in[4 * n * row + 3 * n + col].y;
			w_x                             = w[col + 3 * n].x;
			w_y                             = w[col + 3 * n].y;
			in[4 * n * row + col + 3 * n].x = x * w_x - y * w_y;
			in[4 * n * row + col + 3 * n].y = x * w_y + y * w_x;
		} else {
			double x, y, w_x, w_y;
			x                             = in[n * row + col].x;
			y                             = in[n * row + col].y;
			w_x                           = w[row].x;
			w_y                           = w[row].y;
			in[n * row + col].x           = x * w_x - y * w_y;
			in[n * row + col].y           = x * w_y + y * w_x;
			x                             = in[n * (n + row) + col].x;
			y                             = in[n * (n + row) + col].y;
			w_x                           = w[n + row].x;
			w_y                           = w[n + row].y;
			in[n * (n + row) + col].x     = x * w_x - y * w_y;
			in[n * (n + row) + col].y     = x * w_y + y * w_x;
			x                             = in[n * (2 * n + row) + col].x;
			y                             = in[n * (2 * n + row) + col].y;
			w_x                           = w[2 * n + row].x;
			w_y                           = w[2 * n + row].y;
			in[n * (2 * n + row) + col].x = x * w_x - y * w_y;
			in[n * (2 * n + row) + col].y = x * w_y + y * w_x;
			x                             = in[n * (3 * n + row) + col].x;
			y                             = in[n * (3 * n + row) + col].y;
			w_x                           = w[3 * n + row].x;
			w_y                           = w[3 * n + row].y;
			in[n * (3 * n + row) + col].x = x * w_x - y * w_y;
			in[n * (3 * n + row) + col].y = x * w_y + y * w_x;
		}
	}
}
__global__ void dct4_post(cufftDoubleComplex *in, cufftDoubleReal *out, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			// rows
			double x           = in[4 * n * row + col].x;
			double y           = in[4 * n * row + col].y;
			double result      = x - y;
			out[n * row + col] = result;
		} else {
			// cols
			double x           = in[n * row + col].x;
			double y           = in[n * row + col].y;
			double result      = x - y;
			out[n * row + col] = result;
		}
	}
}

cufftDoubleReal *dct4(cufftDoubleReal *x, int n, int dim)
{
	// allocate fft input and memset
	cufftDoubleComplex *in = nullptr;
	cudaMalloc((void **) &in, sizeof(cufftDoubleComplex) * n * n * 4);
	cudaMemset((void *) in, 0, sizeof(cufftDoubleComplex) * n * n * 4);

	// allocate weight vector
	cufftDoubleComplex *w = nullptr;
	cudaMalloc((void **) &w, 4 * sizeof(cufftDoubleComplex) * n);

	// prepare fft input
	dct4_pre<<<n, n>>>(x, in, w, n, dim);
	dct4_weight<<<n, n>>>(in, w, n, dim);

	// create plan
	cufftHandle plan;
	if (dim == 0) {
		int transform_size = 4 * n;
		int inembed        = 4 * n;
		int istride        = 1;
		int idist          = 4 * n;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	} else {
		int transform_size = 4 * n;
		int inembed        = 4 * n;
		int istride        = n;
		int idist          = 1;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	}

	// allocate fft output
	cufftDoubleComplex *out = nullptr;
	cudaMalloc((void **) &out, sizeof(cufftDoubleComplex) * n * n * 4);

	// perform transform
	int i = cufftExecZ2Z(plan, in, out, CUFFT_FORWARD);
	cufftDestroy(plan);

	// free input
	cudaFree(in);

	// allocate final output
	cufftDoubleReal *retval = nullptr;
	cudaMalloc((void **) &retval, sizeof(cufftDoubleReal) * n * n);

	// get final output
	dct4_post<<<n, n>>>(out, retval, n, dim);

	// free weights and fft output
	cudaFree(w);
	cudaFree(out);

	return retval;
}
__global__ void dst4_pre(cufftDoubleReal *in, cufftDoubleComplex *out, cufftDoubleComplex *w,

                         int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (row == 0) {
			double x = col * M_PI / (4 * n);
			w[col].x = cos(-x) / 2;
			w[col].y = sin(-x) / 2;
			if (col == 0) {
				w[0].x *= 2;
				w[0].y *= 2;
			}
			x                = (n + col) * M_PI / (4 * n);
			w[n + col].x     = cos(-x) / 2;
			w[n + col].y     = sin(-x) / 2;
			x                = (2 * n + col) * M_PI / (4 * n);
			w[2 * n + col].x = -cos(-x) / 2;
			w[2 * n + col].y = -sin(-x) / 2;
			if (col == 0) {
				w[2 * n].x = 0;
				w[2 * n].y = 0;
			}
			x                = (3 * n + col) * M_PI / (4 * n);
			w[3 * n + col].x = -cos(-x) / 2;
			w[3 * n + col].y = -sin(-x) / 2;
		}
		if (dim == 0) {
			// rows
			out[4 * n * row + 2 *(n- col-1) + 1].x         = in[i];
			out[4 * n * row + 4 * n - 2 * (n-col-1) - 1].x = in[i];
		} else {
			// cols
			out[n * (2 * (n-row-1) + 1) + col].x         = in[i];
			out[n * (4 * n - 2 * (n-row-1) - 1) + col].x = in[i];
		}
	}
}

__global__ void dst4_weight(cufftDoubleComplex *in, cufftDoubleComplex *w, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			double x, y, w_x, w_y;
			x                               = in[4 * n * row + col].x;
			y                               = in[4 * n * row + col].y;
			w_x                             = w[col].x;
			w_y                             = w[col].y;
			in[4 * n * row + col].x         = x * w_x - y * w_y;
			in[4 * n * row + col].y         = x * w_y + y * w_x;
			x                               = in[4 * n * row + n + col].x;
			y                               = in[4 * n * row + n + col].y;
			w_x                             = w[col + n].x;
			w_y                             = w[col + n].y;
			in[4 * n * row + col + n].x     = x * w_x - y * w_y;
			in[4 * n * row + col + n].y     = x * w_y + y * w_x;
			x                               = in[4 * n * row + 2 * n + col].x;
			y                               = in[4 * n * row + 2 * n + col].y;
			w_x                             = w[col + 2 * n].x;
			w_y                             = w[col + 2 * n].y;
			in[4 * n * row + col + 2 * n].x = x * w_x - y * w_y;
			in[4 * n * row + col + 2 * n].y = x * w_y + y * w_x;
			x                               = in[4 * n * row + 3 * n + col].x;
			y                               = in[4 * n * row + 3 * n + col].y;
			w_x                             = w[col + 3 * n].x;
			w_y                             = w[col + 3 * n].y;
			in[4 * n * row + col + 3 * n].x = x * w_x - y * w_y;
			in[4 * n * row + col + 3 * n].y = x * w_y + y * w_x;
		} else {
			double x, y, w_x, w_y;
			x                             = in[n * row + col].x;
			y                             = in[n * row + col].y;
			w_x                           = w[row].x;
			w_y                           = w[row].y;
			in[n * row + col].x           = x * w_x - y * w_y;
			in[n * row + col].y           = x * w_y + y * w_x;
			x                             = in[n * (n + row) + col].x;
			y                             = in[n * (n + row) + col].y;
			w_x                           = w[n + row].x;
			w_y                           = w[n + row].y;
			in[n * (n + row) + col].x     = x * w_x - y * w_y;
			in[n * (n + row) + col].y     = x * w_y + y * w_x;
			x                             = in[n * (2 * n + row) + col].x;
			y                             = in[n * (2 * n + row) + col].y;
			w_x                           = w[2 * n + row].x;
			w_y                           = w[2 * n + row].y;
			in[n * (2 * n + row) + col].x = x * w_x - y * w_y;
			in[n * (2 * n + row) + col].y = x * w_y + y * w_x;
			x                             = in[n * (3 * n + row) + col].x;
			y                             = in[n * (3 * n + row) + col].y;
			w_x                           = w[3 * n + row].x;
			w_y                           = w[3 * n + row].y;
			in[n * (3 * n + row) + col].x = x * w_x - y * w_y;
			in[n * (3 * n + row) + col].y = x * w_y + y * w_x;
		}
	}
}
__global__ void dst4_post(cufftDoubleComplex *in, cufftDoubleReal *out, int n, int dim)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		int row = i / n;
		int col = i % n;
		if (dim == 0) {
			// rows
			double x           = in[4 * n * row + col].x;
			double y           = in[4 * n * row + col].y;
			double result      = x - y;
            if(col%2==0){
			out[n * row + col] = result;
            }else{
			out[n * row + col] = -result;
            }
		} else {
			// cols
			double x           = in[n * row + col].x;
			double y           = in[n * row + col].y;
			double result      = x - y;
            if(row%2==0){
			out[n * row + col] = result;
            }else{
			out[n * row + col] = -result;
            }
		}
	}
}

cufftDoubleReal *dst4(cufftDoubleReal *x, int n, int dim)
{
	// allocate fft input and memset
	cufftDoubleComplex *in = nullptr;
	cudaMalloc((void **) &in, sizeof(cufftDoubleComplex) * n * n * 4);
	cudaMemset((void *) in, 0, sizeof(cufftDoubleComplex) * n * n * 4);

	// allocate weight vector
	cufftDoubleComplex *w = nullptr;
	cudaMalloc((void **) &w, 4 * sizeof(cufftDoubleComplex) * n);

	// prepare fft input
	dst4_pre<<<n, n>>>(x, in, w, n, dim);
	dst4_weight<<<n, n>>>(in, w, n, dim);

	// create plan
	cufftHandle plan;
	if (dim == 0) {
		int transform_size = 4 * n;
		int inembed        = 4 * n;
		int istride        = 1;
		int idist          = 4 * n;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	} else {
		int transform_size = 4 * n;
		int inembed        = 4 * n;
		int istride        = n;
		int idist          = 1;
		cufftPlanMany(&plan, 1, &transform_size, &inembed, istride, idist, &inembed, istride, idist,
		              CUFFT_Z2Z, n);
	}

	// allocate fft output
	cufftDoubleComplex *out = nullptr;
	cudaMalloc((void **) &out, sizeof(cufftDoubleComplex) * n * n * 4);

	// perform transform
	int i = cufftExecZ2Z(plan, in, out, CUFFT_FORWARD);
	cufftDestroy(plan);

	// free input
	cudaFree(in);

	// allocate final output
	cufftDoubleReal *retval = nullptr;
	cudaMalloc((void **) &retval, sizeof(cufftDoubleReal) * n * n);

	// get final output
	dst4_post<<<n, n>>>(out, retval, n, dim);

	// free weights and fft output
	cudaFree(w);
	cudaFree(out);

	return retval;
}


__global__ void gpu_div(cufftDoubleReal *f_hat, cufftDoubleReal *denom, int n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n * n) {
		f_hat[i] /= denom[i];
	}
}
CufftSolver::CufftSolver(Domain *dom)
{
	this->d   = dom;
	x_forward = &dst2;
	x_inverse = &idst2;
	y_forward = &dst2;
	y_inverse = &idst2;
	if (d->neumann) {
		if (!d->hasNbr(Side::east) && !d->hasNbr(Side::west)) {
			x_forward = &dct2;
			x_inverse = &idct2;
		} else if (!d->hasNbr(Side::west)) {
			x_forward = &dct4;
			x_inverse = &dct4;
		} else if (!d->hasNbr(Side::east)) {
			x_forward = &dst4;
			x_inverse = &dst4;
		}
		if (!d->hasNbr(Side::north) && !d->hasNbr(Side::south)) {
			y_forward = &dct2;
			y_inverse = &idct2;
		} else if (!d->hasNbr(Side::south)) {
			y_forward = &dct4;
			y_inverse = &dct4;
		} else if (!d->hasNbr(Side::north)) {
			y_forward = &dst4;
			y_inverse = &dst4;
		}
	}
	if (denoms.count(d)) {
		denom_ptr = &denoms[d];
	} else {
		denom_ptr               = &denoms[d];
		valarray<double> &denom = *denom_ptr;
		denom.resize(d->n * d->n);
		// create denom vector
		if (!d->neumann) {
			for (int yi = 0; yi < d->n; yi++) {
				denom[slice(yi * d->n, d->n, 1)]
				= -4 / (d->h_x * d->h_x) * pow(sin((yi + 1) * M_PI / (2 * d->n)), 2);
			}
		} else {
			if (!d->hasNbr(Side::north) && !d->hasNbr(Side::south)) {
				for (int yi = 0; yi < d->n; yi++) {
					denom[slice(yi * d->n, d->n, 1)]
					= -4 / (d->h_x * d->h_x) * pow(sin(yi * M_PI / (2 * d->n)), 2);
				}
			} else if (!d->hasNbr(Side::south) || !d->hasNbr(Side::north)) {
				for (int yi = 0; yi < d->n; yi++) {
					denom[slice(yi * d->n, d->n, 1)]
					= -4 / (d->h_x * d->h_x) * pow(sin((yi + 0.5) * M_PI / (2 * d->n)), 2);
				}
			} else {
				for (int yi = 0; yi < d->n; yi++) {
					denom[slice(yi * d->n, d->n, 1)]
					= -4 / (d->h_x * d->h_x) * pow(sin((yi + 1) * M_PI / (2 * d->n)), 2);
				}
			}
		}

		valarray<double> ones(d->n);
		ones = 1;

		if (!d->neumann) {
			for (int xi = 0; xi < d->n; xi++) {
				denom[slice(xi, d->n, d->n)]
				-= 4 / (d->h_y * d->h_y) * pow(sin((xi + 1) * M_PI / (2 * d->n)), 2) * ones;
			}
		} else {
			if (!d->hasNbr(Side::east) && !d->hasNbr(Side::west)) {
				for (int xi = 0; xi < d->n; xi++) {
					denom[slice(xi, d->n, d->n)]
					-= 4 / (d->h_y * d->h_y) * pow(sin(xi * M_PI / (2 * d->n)), 2) * ones;
				}
			} else if (!d->hasNbr(Side::west) || !d->hasNbr(Side::east)) {
				for (int xi = 0; xi < d->n; xi++) {
					denom[slice(xi, d->n, d->n)]
					-= 4 / (d->h_y * d->h_y) * pow(sin((xi + 0.5) * M_PI / (2 * d->n)), 2) * ones;
				}
			} else {
				for (int xi = 0; xi < d->n; xi++) {
					denom[slice(xi, d->n, d->n)]
					-= 4 / (d->h_y * d->h_y) * pow(sin((xi + 1) * M_PI / (2 * d->n)), 2) * ones;
				}
			}
		}
	}
}
CufftSolver::~CufftSolver()
{
	/*
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
	*/
}
void CufftSolver::solve()
{
	valarray<double> f_copy(d->n * d->n);
	f_copy[slice(0, d->n * d->n, 1)] = d->f;
	if (!d->hasNbr(Side::north) && d->neumann) {
		f_copy[slice(d->n * (d->n - 1), d->n, 1)] -= 1 / d->h_y * d->boundary_north;
	} else {
		f_copy[slice(d->n * (d->n - 1), d->n, 1)] -= 2 / (d->h_y * d->h_y) * d->boundary_north;
	}
	if (!d->hasNbr(Side::east) && d->neumann) {
		f_copy[slice((d->n - 1), d->n, d->n)] -= 1 / d->h_x * d->boundary_east;
	} else {
		f_copy[slice((d->n - 1), d->n, d->n)] -= 2 / (d->h_x * d->h_x) * d->boundary_east;
	}
	if (!d->hasNbr(Side::south) && d->neumann) {
		f_copy[slice(0, d->n, 1)] += 1 / d->h_y * d->boundary_south;
	} else {
		f_copy[slice(0, d->n, 1)] -= 2 / (d->h_y * d->h_y) * d->boundary_south;
	}
	if (!d->hasNbr(Side::west) && d->neumann) {
		f_copy[slice(0, d->n, d->n)] += 1 / d->h_x * d->boundary_west;
	} else {
		f_copy[slice(0, d->n, d->n)] -= 2 / (d->h_x * d->h_x) * d->boundary_west;
	}

	// copy f to gpu
	cufftDoubleReal *dev_f = nullptr;
	cudaMalloc((void **) &dev_f, sizeof(cufftDoubleReal) * d->n * d->n);
	cudaMemcpy(dev_f, &f_copy[0], sizeof(cufftDoubleReal) * d->n * d->n, cudaMemcpyHostToDevice);

	// copy denom to gpu
	cufftDoubleReal *dev_denom = nullptr;
	cudaMalloc((void **) &dev_denom, sizeof(cufftDoubleReal) * d->n * d->n);
	cudaMemcpy(dev_denom, &(*denom_ptr)[0], sizeof(cufftDoubleReal) * d->n * d->n,
	           cudaMemcpyHostToDevice);

	// transform
	cufftDoubleReal *row_trans = x_forward(dev_f, d->n, 0);
	cufftDoubleReal *f_hat     = y_forward(row_trans, d->n, 1);
	cudaFree(row_trans);
	cudaFree(dev_f);

	// divide
	gpu_div<<<d->n, d->n>>>(f_hat, dev_denom, d->n);
	cudaFree(dev_denom);

	if (d->neumann
	    && !(d->hasNbr(Side::north) || d->hasNbr(Side::east) || d->hasNbr(Side::south)
	         || d->hasNbr(Side::west))) {
		cudaMemset((void *) f_hat, 0, sizeof(cufftDoubleReal));
	}

	// inverse transform
	row_trans              = x_inverse(f_hat, d->n, 0);
	cufftDoubleReal *dev_u = y_inverse(row_trans, d->n, 1);
	cudaFree(row_trans);
	cudaFree(f_hat);

	// copy solution to host
	cudaMemcpy(&d->u[0], dev_u, sizeof(cufftDoubleReal) * d->n * d->n, cudaMemcpyDeviceToHost);
	cudaFree(dev_u);

	// d->u /= (4.0 * d->n * d->n);

	/*if (d->ds.zero_patch) {
	    d->u -= d->u.sum() / d->u.size();
	}*/
}
