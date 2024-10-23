#include <cuda_runtime.h>
#include <cublas.h>
#include <cuda_fp16.h>
#include <stdio.h>

#define CHUNKSIZE 32

__global__ void float2double_GPU_kernel( int nrows, int ncols,
                const float *F, int ldf,
                double *D, int ldh ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

        D[idy*ldh+idx]= (double)F[idy*ldf+idx];
	//printf("D %d %d : %d %d : %lf\n", idx, idy, nrows, ncols, D[idy*ldh+idx]);
	//printf("F %d %d : %d %d : %f\n", idx, idy, nrows, ncols, F[idy*ldh+idx]);
}

extern "C"
void float2double_GPU( int nrows, int ncols,
                const float *F, int ldf,
                double *D, int ldh,
                cudaStream_t stream ) {
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);
        float2double_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, F, ldf, D, ldh);
}

/****************************************************************************************************/

__global__ void double2float_GPU_kernel( int nrows, int ncols,
                const double *D, int ldh,
                float *F, int ldf ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

        F[idy*ldf+idx]=__double2float_rn( D[idy*ldh+idx] ); 
	//printf("D %d %d : %d %d : %lf\n", idx, idy, nrows, ncols, D[idy*ldh+idx]);
	//printf("F %d %d : %d %d : %f\n", idx, idy, nrows, ncols, F[idy*ldh+idx]);
}

extern "C"
void double2float_GPU( int nrows, int ncols,
                const double *D, int ldh,
                float *F, int ldf,
                cudaStream_t stream ) {

        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);
        double2float_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, D, ldh, F, ldf);

}

/****************************************************************************************************/

__global__ void half2float_GPU_kernel( int nrows, int ncols,
                const __half *H, int ldh,
                float *F, int ldf ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

        F[idy*ldf+idx]= __half2float( H[idy*ldh+idx] );
}

extern "C"
void half2float_GPU( int nrows, int ncols,
                const void *_H, int ldh,
                float *F, int ldf,
                cudaStream_t stream ) {
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);
        __half *H = (__half *)_H;
        half2float_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, H, ldh, F, ldf);
}

/****************************************************************************************************/

__global__ void float2half_GPU_kernel( int nrows, int ncols,
                const float *F, int ldf,
                __half *H, int ldh ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

        H[idy*ldh+idx]= __float2half_rn( F[idy*ldf+idx] );
}

extern "C"
void float2half_GPU( int nrows, int ncols,
                const float *F, int ldf,
                void *_H, int ldh,
                cudaStream_t stream ) {
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);
        __half *H = (__half *)_H;
        float2half_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, F, ldf, H, ldh);
}


/****************************************************************************************************/

__global__ void memcpy_half_GPU_kernel( int nrows, int ncols,
                __half *src, __half *dest ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

        dest[idy*nrows+idx] = src[idy*nrows+idx];
}

extern "C"
void memcpy_half_GPU( int nrows, int ncols, void *_src, void *_dest, cudaStream_t stream ) { 
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);
        __half *src = (__half *)_src;
        __half *dest = (__half *)_dest;
        memcpy_half_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, src, dest);
}

/****************************************************************************************************/

__global__ void memcpy_float_GPU_kernel( int nrows, int ncols,
               float *src, float *dest ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;
                           
        if( idx >= nrows || idy >= ncols ) { return; }
                           
        dest[idy*nrows+idx] = src[idy*nrows+idx];
}                          
                           
extern "C"                 
void memcpy_float_GPU( int nrows, int ncols, void *_src, void *_dest, cudaStream_t stream ) {                         
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);
        float *src = (float *)_src;
        float *dest = (float *)_dest;
        memcpy_float_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, src, dest);
}

/****************************************************************************************************/

__global__ void matrix_print_float_GPU_kernel( int nrows, int ncols, float *A ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

	printf("SINGLE_PRINT %d %d : %g\n", idx, idy, A[idy*nrows+idx]);
}

extern "C"
void matrix_print_float_GPU( int nrows, int ncols, float *A, cudaStream_t stream ) {
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);

	matrix_print_float_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, A);
}

/****************************************************************************************************/

__global__ void matrix_print_double_GPU_kernel( int nrows, int ncols, double *A ) {
        const int tx=threadIdx.x;
        const int ty=threadIdx.y;
        const int idx= blockIdx.x * blockDim.x + tx;
        const int idy= blockIdx.y * blockDim.y + ty;

        if( idx >= nrows || idy >= ncols ) { return; }

        printf("DOUBLE_PRINT %d %d : %g\n", idx, idy, A[idy*nrows+idx]);
}

extern "C"
void matrix_print_double_GPU( int nrows, int ncols, double *A, cudaStream_t stream ) {
        int nBlockx= (nrows+CHUNKSIZE-1)/CHUNKSIZE;
        int nBlocky= (ncols+CHUNKSIZE-1)/CHUNKSIZE;
        dim3 dimBlock(CHUNKSIZE, CHUNKSIZE);
        dim3 dimGrid(nBlockx, nBlocky);

        matrix_print_double_GPU_kernel<<<dimGrid, dimBlock, 0, stream>>>(nrows, ncols, A);
}

/****************************************************************************************************/
extern "C"
cublasStatus_t my_cublasGemmEx(cublasHandle_t handle,
                           cublasOperation_t transa,
                           cublasOperation_t transb,
                           int m,
                           int n,
                           int k,
                           const void    *alpha,
                           const void     *A, 
                           cudaDataType_t Atype,
                           int lda,
                           const void     *B,
                           cudaDataType_t Btype,  
                           int ldb,
                           const void    *beta, 
                           void           *C,
                           cudaDataType_t Ctype,
                           int ldc,
                           cudaDataType_t computeType,
                           cublasGemmAlgo_t algo)
{
    return cublasGemmEx(handle, transa, transb, m, n, k,  
                        (__half *)alpha, A, Atype, lda,
                                         B, Btype, ldb,
                        (__half *)beta,  C, Ctype, ldc,
                        computeType, algo);
}

/****************************************************************************************************/

__global__ void dcmg_array_GPU_kernel(double *A, int m, int n, int m0,
        int n0, double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
        double localtheta0, double localtheta1, double localtheta2, int distance_metric, int lda)
{
	const int tx  = threadIdx.x;
	const int ty  = threadIdx.y;
	const int idx = blockIdx.x * blockDim.x + tx;
	const int idy = blockIdx.y * blockDim.y + ty;
	if(idx>=m || idy >=n){return;}

	//double x0, y0;
	double expr  = 0.0;
	double expr1 = 0.0;

	double sigma_square = localtheta0;// * localtheta[0];

	expr = sqrt(pow((l2_x_cuda[idx] - l1_x_cuda[idy]), 2) +
			pow((l2_y_cuda[idx] - l1_y_cuda[idy]), 2));

	expr1 = pow(expr, localtheta2);
	if(expr == 0)
		A[idx + idy * lda] = sigma_square /*+ 1e-4*/;
	else
		A[idx + idy * lda] = sigma_square *  exp(-(expr1/localtheta1)); // power-exp kernel



}

extern "C"
void dcmg_array_GPU( double *A, int m, int n, int m0,
		int n0, double* l1_x_cuda, double* l1_y_cuda, double* l2_x_cuda, double* l2_y_cuda,
		double *localtheta, int distance_metric, int lda, cudaStream_t stream){

	int nBlockx= (m+CHUNKSIZE-1)/CHUNKSIZE;
	int nBlocky= (n+CHUNKSIZE-1)/CHUNKSIZE;
	dim3 dimBlock(CHUNKSIZE,CHUNKSIZE);
	dim3 dimGrid(nBlockx,nBlocky);

	dcmg_array_GPU_kernel<<<dimGrid,dimBlock,0,stream>>>(A, m, n, m0, n0, l1_x_cuda, l1_y_cuda, l2_x_cuda, l2_y_cuda, localtheta[0],localtheta[1],localtheta[2], distance_metric, lda);
}
