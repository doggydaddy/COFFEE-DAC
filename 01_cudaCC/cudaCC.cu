// standard includes
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

// nifticlib includes
#include <nifti1.h>
#include <fslio.h>
#include <nifti1_io.h>
#include <math.h>

// for CUDA computing, including the following ...
#include <cuda.h>
#include <cuda_runtime.h> 

/* ----------------- */
/* support functions */
/* ----------------- */

/* return the total system memory (*nix systems )*/
unsigned long long getTotalSystemMemory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages*page_size;
}

/* return "Good-Values", or number of voxels inside the mask */
size_t getGoodValues(double ***mask, int dimx, int dimy, int dimz)
{
    int i, j, k;
    size_t gV = 0;
    for ( i=0; i<dimx; ++i) {
        for ( j=0; j<dimy; ++j ) {
            for (k=0; k<dimz; ++k) {
                if ( mask[i][j][k] != 0 )
                    gV++;
            }
        }
    }
    return gV;
}

/* allocate matrix of r rows and c columns */
float** createMatrix(int r, int c) {
    float **M;
    int i;
    M = (float**)malloc(r*sizeof(float*));
    if ( M==NULL) {
        printf("Error! Failed to allocate memory on host!\n");
        exit(-1);
    }
    for ( i=0; i<r; ++i ) {
        M[i] = (float*)malloc(c*sizeof(float));
        if ( M[i]==NULL ) {
            printf("Error! Failed to allocate memory on host!");
            exit(-1);
        }
    }
    return M;
}

/* reformat 4D array as obtained from fsl nifti-io library to 1d array
 * 
 * output 1d array is masked with what is only inside mask
 * 
 * indexing occurs from last to second index, i.e. i first, then j, finally k,
 * where input 4D data is input_data[time][k][j][i]
 */
void inputData4DtoArray(double ****input4D, float *h_arr, double ***mask, int dimx, int dimy, int dimz, int dimt)
{
    int i, j, k, t;
    int c = 0;
    for (k=0; k<dimz; ++k) {
        for (j=0; j<dimy; ++j) {
            for (i=0; i<dimx; ++i) {
                if (mask[k][j][i] != 0) {
                    for (t=0; t<dimt; ++t) {
                        h_arr[c] = (float)input4D[t][k][j][i];
                        c++;
                    }
                }
            }
        }
    }
}

/* save (2D full) cc-matrix to text file */
void saveToText_FULL(float **outputData, size_t gV, char *fileName)
{
    int i, j;
    FILE *output = fopen(fileName, "w");
    for( i=0; i<gV; ++i ) {
        for( j=0; j<gV; ++j ) {
            fprintf(output, "%lf ", outputData[i][j]);
        }
        fprintf(output, "\n");
    }
    fclose(output);
}

/* save (upper triangular without diagnoal) cc-matrix to text file */
void saveToText_TRIA(float *outputData, size_t gV, char *fileName)
{
    int i, j;
    FILE *output = fopen(fileName, "w");
	size_t c = 0;
    for( i=0; i<gV; ++i ) {
        for( j=i+1; j<gV; ++j ) {
            fprintf(output, "%lf ", outputData[c]);
			c++;
        }
        fprintf(output, "\n");
    }
    fclose(output);
}

/* converts input upper-triangular matrix without diagnoal to full 2D matrix
 * (for output to text file).
 * 
 * forces diagonal to be 0 (instead of 1 for pure cc-matrix)
 */
void upperTriangularToFull(float *input, float **output, size_t gV)
{
    int i, j;
    size_t c = 0;
    for( i=0; i<gV; ++i ) {
        for( j=i; j<gV; ++j ) {
            if ( i == j ) {
                output[i][j] = 0;
            } else {
                output[i][j] = input[c];
                output[j][i] = input[c];
                c++;
            }
        }
    }
}

/* 2D CUDA Kernel for calculating the cross-correlation matrix 
 * of data with size (gV x nt). 
 */
__global__ void CUDA_CC(float *data, float *answer, size_t gV, int nt) 
{
	float mean_t1 = 0;
	float mean_t2 = 0;
	float nom = 0;
	float de1 = 0;
	float de2 = 0;
	size_t k = 0;

	int n = blockIdx.x * blockDim.x + threadIdx.x;
	int m = blockIdx.y * blockDim.y + threadIdx.y;

	if ( n >= gV || m >= gV ) return;
	if ( m <= n ) return; 

	for (int i=0; i<nt; ++i) 
	{
		mean_t1 += data[nt*n+i];	
		mean_t2 += data[nt*m+i];
	}
	mean_t1 /= nt;
	mean_t2 /= nt;

	for (int i=0; i<nt; ++i) 
	{
		nom += (data[nt*n+i]-mean_t1)*(data[nt*m+i]-mean_t2);
		de1 += (data[nt*n+i]-mean_t1)*(data[nt*n+i]-mean_t1);
		de2 += (data[nt*m+i]-mean_t2)*(data[nt*m+i]-mean_t2);
	}

	/* 2D array index to upper triangular without diagonal */
	k = (gV*(gV-1)/2)-(gV-n)*((gV-n)-1)/2+m-n-1;
	//k = gV*n+m;
	answer[k] = 1 - ( nom / sqrt(de1*de2) ); 
}

/* main function */
int main ( int argc , char * argv[] ) 
{
	/* 
	 * -------------------------------------------------------
     * Initialize, parse inputs, and calculate memory required
	 * -------------------------------------------------------
	 */

    /* parsing input parametres */
    char *inputFile, *maskFile, *outputFile;
    inputFile = argv[1];
    maskFile = argv[2];
    outputFile = argv[3];

    /* loading datasets and mask, fsl_io library is used */ 
    FSLIO *fslio = FslInit();
    void *buffer = FslReadAllVolumes(fslio, maskFile);
    double ***mask = FslGetVolumeAsScaledDouble(fslio, 0);
    int host_err = FslClose(fslio);
    /* loading the input dataset */
    fslio = FslInit();
    buffer = FslReadAllVolumes(fslio, inputFile);
    double ****data = FslGetBufferAsScaledDouble(fslio);
    int fsl_err = FslClose(fslio);
    int nx_data = fslio->niftiptr->nx;
    int ny_data = fslio->niftiptr->ny;
    int nz_data = fslio->niftiptr->nz;
    int nt_data = fslio->niftiptr->nt;

	/* Error code to check my CUDA code */
	cudaError_t err = cudaSuccess;		

	/* Calculate sizes */
	/* calculate number of voxels inside the mask */
    size_t gV = getGoodValues(mask, nx_data, ny_data, nz_data);
	printf("number of voxels within the mask: %zu\n", gV);

	unsigned long num_elem_input = gV*nt_data;
	unsigned long num_elem_output = (unsigned long)((gV*(gV-1))/2);
	size_t input_size = num_elem_input*sizeof(float);
	size_t output_size = num_elem_output*sizeof(float);
	printf("Number of elements in input: %lu, or %zu bytes of memory in floats.\n", num_elem_input, input_size);
	printf("Number of elements in output: %lu, or %zu bytes of memory in floats.\n", num_elem_output, output_size);

	/* print total memory needed */
	double up3zeros = 1024.0;
	long double input_size_mb = (num_elem_input*sizeof(float)) / (up3zeros*up3zeros); 
	long double output_size_gb = (num_elem_output*sizeof(float)) / (up3zeros*up3zeros*up3zeros); 
	printf("Size of the input is %LFMB\n", input_size_mb);
	printf("Size of output cc-matrix is %LFGB\n", output_size_gb);

	/* see if we can divide up the job into parts if the job is too big */
	unsigned long long system_memory = getTotalSystemMemory();
	if ( input_size*3+output_size > system_memory*0.9 ) 
	{
		printf("The job is too big to fit onto the system memory,\n");
		printf("please upgrade your potato computer and try again!\n");
		return -1;
	}
	size_t device_free, device_total;
	cudaMemGetInfo(&device_free, &device_total);
	printf("GPU free mem: %lu, Total mem: %lu\n", device_free, device_total);
	if ( device_free < (input_size+output_size) ) 
	{
		printf("There isn't sufficient GPU memory to calculate everything in one go,\n");
		printf("please consider splitting the job into parts.\n");
		return -1;
	}

	/* allocate memeory on host */
	float *h_data = (float *)malloc(input_size);
	float *gpu_ccmat = (float *)malloc(output_size);
	if ( h_data == NULL || gpu_ccmat == NULL ) 
	{
		fprintf(stderr, "Failed to allocate memory on host!\n");
		exit(EXIT_FAILURE);
	}

    /* Parsing input data into appropriate format on host */
    inputData4DtoArray(data, h_data, mask, nx_data, ny_data, nz_data, nt_data);
	
	/* Allocating the host data input matrix to GPU */
	printf("Allocating device memory ...\n");
	float *d_data = NULL;
	float *d_result = NULL;
	err = cudaMalloc((void **)&d_data, input_size);
	if (err != cudaSuccess) 
	{
		printf("cudaMalloc d_data returned error %s (code %d, line(%d)\n", cudaGetErrorString(err), err, __LINE__);
		exit(EXIT_FAILURE);
	}
	err = cudaMalloc((void **)&d_result, output_size);
	if (err != cudaSuccess) 
	{
		printf("cudaMalloc d_result returned error %s (code %d, line(%d)\n", cudaGetErrorString(err), err, __LINE__);
		exit(EXIT_FAILURE);
	}
	printf("... done! \n");
	
	/* copy the host input data to device input data */
	printf("Copy input data from the host memory to the CUDA device\n");
	err = cudaMemcpy(d_data, h_data, input_size, cudaMemcpyHostToDevice);
	if ( err != cudaSuccess ) 
	{
		printf("Failed to copy input data from host to device (error code %s)!\n", cudaGetErrorString(err));
	}
	printf(" ... done!\n");


	/* -------------------------------------------------------
     * Setting up, and performing calculations using CUDA
	 * -------------------------------------------------------
	 */

	dim3 dimBlock(16, 16);
	dim3 dimGrid;
	dimGrid.x = (gV + dimBlock.x - 1) / dimBlock.x;
	dimGrid.y = (gV + dimBlock.y - 1) / dimBlock.y;

	clock_t start_gpu = clock(), diff;
	
	/* running the cuda kernel */
	printf("Initializing CUDA kernel ...\n");
	CUDA_CC<<<dimGrid, dimBlock>>>(d_data, d_result, gV, nt_data);

	/* copy the results form the device back to host */
	err = cudaMemcpy(gpu_ccmat, d_result, output_size, cudaMemcpyDeviceToHost);
	if ( err != cudaSuccess ) 
	{
		printf("Failed to copy result from device to host (error code %s)!\n", cudaGetErrorString(err));
	}

	diff = clock() - start_gpu;
	int msec = diff * 1000 / CLOCKS_PER_SEC;

	printf(" ... done! \n ");
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

	/* -------------------------------------------------------
	 * saving output
	 * -------------------------------------------------------
	 */
	printf("Saving output ...  \n");
    //float **pmat = createMatrix(gV, gV);
    //upperTriangularToFull(gpu_ccmat, pmat, gV);
    saveToText_TRIA(gpu_ccmat, gV, outputFile);
	printf(" ... done! \n ");

	/* -------------------------------------------------------
	 * Doing some cleanup ...
	 * -------------------------------------------------------
	 */
	err = cudaFree(d_data);
	err = cudaFree(d_result);
	if ( err != cudaSuccess ) 
	{
		fprintf(stderr, "Failed to free device data (error code %s)!\n", cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
	free(h_data);
	free(gpu_ccmat);
	
	return 0;
}

