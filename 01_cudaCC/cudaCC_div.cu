// standard includes
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

// nifticlib includes
#include <nifti1.h>
#include <fslio.h>
#include <nifti1_io.h>
#include <math.h>

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
    size_t c = 0;
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

void inputData4DtoArrayWithRanges(double ****input4D, float *h_arr, double ***mask, 
                                  int dimx, int dimy, int dimz, int dimt, 
                                  size_t start_idx, size_t end_idx)
{
    int i, j, k, t;
    size_t c = 0; // index for output "current" index
    size_t d = 0; // index for "n-th voxel within mask"
    for (k=0; k<dimz; ++k) 
    {
        for (j=0; j<dimy; ++j) 
        {
            for (i=0; i<dimx; ++i) 
            {
                if (mask[k][j][i] != 0) 
                {
                    if (d >= start_idx && d <= end_idx) 
                    {
                        for (t=0; t<dimt; ++t) 
                        {
                            h_arr[c] = (float)input4D[t][k][j][i];
                            c++;
                        }
                    }
                    d++;
                }
            }
        }
    }
}

int inputData4DtoArrayWith2Ranges(double ****input4D, float *h_arr, double ***mask, 
                                  int dimx, int dimy, int dimz, int dimt, 
                                  size_t start_idx1, size_t end_idx1, 
                                  size_t start_idx2, size_t end_idx2)
{
    int i, j, k, t;
    size_t c = 0; // index for output "current" index
    size_t d = 0; // index for "n-th voxel within mask"
    int div;
    for (k=0; k<dimz; ++k) 
    {
        for (j=0; j<dimy; ++j) 
        {
            for (i=0; i<dimx; ++i) 
            {
                if (mask[k][j][i] != 0) 
                {
                    if ((d >= start_idx1 && d <= end_idx1) || (d >= start_idx2 && d <= end_idx2)) 
                    {
                        if (d == start_idx2) {
                            div = c;
                        }
                        for (t=0; t<dimt; ++t) 
                        {
                            h_arr[c] = (float)input4D[t][k][j][i];
                            c++;
                        }
                    }
                    d++;
                }
            }
        }
    }
    return div;
}

/* save (2D full) cc-matrix to text file */
void saveToText_RECT(float *outputData, size_t nrows, size_t ncols, char *fileName)
{
    int i, j;
    FILE *output = fopen(fileName, "w");
    for( i=0; i<nrows; ++i ) {
        for( j=0; j<ncols; ++j ) {
            fprintf(output, "%lf ", outputData[i*ncols+j]);
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
                output[i][j] = 1.0;
            } else {
                output[i][j] = input[c];
                output[j][i] = input[c];
                c++;
            }
        }
    }
}

/* retrieves n-th voxel from the from data of dimensions (nV x time) */
void getNthSeries(float *data, float *output, int n, size_t gV, int nt) 
{
    int c;
    int i;

	c = 0;
    for (i=n*nt; i<n*nt+nt; ++i) 
	{
        output[c] = data[i];
        c++; // all good cpp code needs at least one "c++"
    }
}

/*
 * calculate the cross-correlation matrix of data of size (gV x nt)
 * at n-th and m-th indices.
 * 
 * implementated for ease of reading and omp extensions.
 * 
 * there is a faster (fft-based) implementation, but this is only to check with
 * cuda results so, good enough?
*/
float calcCrossCorr_omp(float *data, int n, int m, size_t gV, int nt)
{
    float t1[nt]; 
    float t2[nt];
    getNthSeries(data, t1, n, gV, nt);
    getNthSeries(data, t2, m, gV, nt);

    float m1 = 0; 
    float m2 = 0;
    for (int i=0; i<nt; ++i) 
	{
        m1 += t1[i];
        m2 += t2[i];
    }
    m1 /= nt; 
    m2 /= nt;    

    float nom = 0; 
	float de1 = 0; 
	float de2 = 0;
    for (int i=0 ; i<nt ; ++i) 
	{
        nom += (t1[i] - m1) * (t2[i] - m2);     
        de1 += (t1[i] - m1) * (t1[i] - m1);  
        de2 += (t2[i] - m2) * (t2[i] - m2);
    }

    float output;
    output = nom / ( sqrt(de1*de2) ); 

    return output;
}

void calcCrossCorr_DIAG(float *data, float *result, size_t gV, int nt)
{
    printf("calculating diagonal block of size (%zu, %zu)\n", gV, gV); 
    #pragma omp parallel for
    for ( int i=0; i<gV; ++i ) 
	{
        for ( int j=i+1; j<gV; ++j ) 
		{
		    size_t k = (gV*(gV-1)/2)-(gV-i)*((gV-i)-1)/2+j-i-1;
            result[k] = calcCrossCorr_omp(data, i, j, gV, nt);  
        }
    }
}

void calcCrossCorr_OFFD(float *data, float *result, size_t gV, int nt, size_t div)
{
    size_t rows, cols;
    rows = div;
    cols = gV-div;
    printf("calculating off-diagonal block of size (%zu, %zu)\n", rows, cols); 
    #pragma omp parallel for
    for (size_t i=0; i<rows; ++i ) 
	{
        for (size_t j=0; j<cols; ++j ) 
		{
            result[(i*cols)+j] = calcCrossCorr_omp(data, i, j+rows, gV, nt);  
        }
    }
}

int sanityCheck(float *input, size_t input_size) 
{
    for(int i=0;i<input_size;++i)
    {
        if (input[i] == 0 || input[i] >= 1 || input[i] <= -1) 
        {
            printf("sanity check failed, at index %d with %f\n", i, input[i]);
            return 0;
        }
    }
    return 1;
}

/* 2D CUDA Kernel for calculating the cross-correlation matrix 
 * of data with size (gV x nt). 
 */
__global__ void CUDA_CC_DIAG(float *data, float *answer, size_t gV, int nt) 
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
	answer[k] = nom / ( sqrt(de1*de2) ); 
}

/* 2D CUDA Kernel for calculating the cross-correlation matrix 
 * of data with size (gV x nt). 
 */
__global__ void CUDA_CC_OFFD(float *data, float *answer, size_t gV, int nt, size_t div) 
{
	float mean_t1 = 0;
	float mean_t2 = 0;
	float nom = 0;
	float de1 = 0;
	float de2 = 0;

	int n = blockIdx.x * blockDim.x + threadIdx.x;
	int m = blockIdx.y * blockDim.y + threadIdx.y;

	if ( n >= gV || m >= gV ) return;

    // boundary conditions:
    // n(i) = [0;div)
    // m(j) = [div;gV)
    if ( n >= div || m < div ) return;

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

	answer[(n*(gV-div))+(m-div)] = nom / ( sqrt(de1*de2) ); 
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
	size_t up3zeros = 1024;
	long double input_size_mb = (num_elem_input*sizeof(float)) / (up3zeros*up3zeros); 
	long double output_size_gb = (num_elem_output*sizeof(float)) / (up3zeros*up3zeros*up3zeros); 
	printf("Size of the input is %LFMB\n", input_size_mb);
	printf("Size of output cc-matrix is %LFGB\n", output_size_gb);

    /* total memory we have */
    printf("Getting system memory availability:\n");
    printf("===================================\n");
    /* host memory */
	unsigned long long system_memory = getTotalSystemMemory();
    /* device memory */
	size_t device_free_mem, device_total_mem;
	cudaMemGetInfo(&device_free_mem, &device_total_mem);
    printf("CPU available memory: %llu bytes\n", system_memory);
	printf("GPU free mem: %lu, Total mem: %lu\n", device_free_mem, device_total_mem);

    /* If we have enough host memory and attempt to allocate output buffer,
     * If we do not have enough system memory, complain and exit since the
     * program will not make things faster. */

    /* buffer output image */ 
    printf("allocate buffer memory on host ...\n");
    float *output_buffer = (float*)malloc(output_size);
    if ( output_buffer == NULL ) 
    {
        fprintf(stderr, "Failed to allocate output buffer on host memory!\n");
        printf("The job is too big to fit onto the system (host) memory,\n");
        printf("please upgrade your potato computer and try again!\n");
        return EXIT_FAILURE;
    }

    /* For device memory:
     * We need approximately input_size + output_size 
     * to compulte everything in one go */
    size_t device_req_mem = (size_t)(input_size+output_size);

    float num_divisions = 1.0;
    float num_runs = 1.0;
	if ( device_req_mem > device_free_mem ) 
	{
		printf("The job is too big to fit onto the device memory all at once ...\n");
		printf("... so we need to divide the computation into blocks.\n");

        size_t max_blk_size = (size_t)((sqrtf(device_free_mem-input_size))/sizeof(float));
        num_divisions = ceil((float)gV/(float)max_blk_size);
        num_runs = (num_divisions*(num_divisions+1))/2;
        printf("so you need at least %d runs to complete the entire job.\n", (int)num_runs);
        printf("You can process, at most, %zu indices at once.\n", max_blk_size);
        printf("We will split the voxels into %d divisions.\n", (int)num_divisions);
        printf("So you need %d runs to complete the entire job.\n", (int)num_runs);

        int start_idx[(int)num_divisions];
        int end_idx[(int)num_divisions];
        for (int i=0; i<num_divisions; ++i)
        {
            start_idx[i] = i*max_blk_size;
            end_idx[i] = ((i*max_blk_size)+max_blk_size)-1;
            if (end_idx[i] >= gV) 
            {
                end_idx[i] = gV - 1;
                continue;
            } 
        }
        printf("[DBG]: Printing out the division indices:\n");
        for (int i=0; i<num_divisions; ++i)
        {
            printf("division %d, [%d, %d]\n", i, start_idx[i], end_idx[i]);
        }

        /* ========================= */
        /* performing the block runs */
        /* ========================= */

        /* allocate variables needed for block runs */
        float *h_data;
        float *gpu_ccmat;
        size_t Fi, Fj, Fk, bk;
        size_t blk_size_x, blk_size_y, blk_size;
        size_t range_x, range_y, range_n;
        int current_run = 0;
        char s[26]; // for outputing blocks into files
        clock_t start_gpu = clock(), diff; // start the clock!
        for (int blk_x = 0; blk_x < num_divisions; ++blk_x) 
        {
            for (int blk_y = blk_x; blk_y < num_divisions; ++blk_y) 
            {
                if (blk_x == blk_y) 
                {
                    printf("Current run number: %d\n", current_run);
                    printf("Block indices (blk_x, blk_y) = (%d, %d)\n", blk_x, blk_y);

                    blk_size = end_idx[blk_x] - start_idx[blk_x] + 1;
                    printf("This is a diagonal block with block size: %zu\n", blk_size);

                    input_size = blk_size*nt_data*sizeof(float);
                    output_size = ((blk_size*(blk_size-1))/2)*sizeof(float);

                    /* allocate memory on host */
                    h_data = (float *)malloc(input_size);
                    gpu_ccmat = (float *)malloc(output_size);

                    printf("Parsing input data into appropriate format on host\n");
                    inputData4DtoArrayWithRanges(data, h_data, mask, 
                                                 nx_data, ny_data, nz_data, nt_data, 
                                                 start_idx[blk_x], end_idx[blk_x]);
                    if ( h_data == NULL || gpu_ccmat == NULL ) 
                    {
                        fprintf(stderr, "Failed to allocate memory on host!\n");
                        return EXIT_FAILURE;
                    }
                    printf("...done!\n");

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
                    /* running the cuda kernel */
                    printf("Initializing CUDA (DIAG) kernel ...\n");
                    CUDA_CC_DIAG<<<dimGrid, dimBlock>>>(d_data, d_result, blk_size, nt_data);

                    /* copy the results form the device back to host */
                    err = cudaMemcpy(gpu_ccmat, d_result, output_size, cudaMemcpyDeviceToHost);
                    if ( err != cudaSuccess ) 
                    {
                        printf("Failed to copy result from device to host (error code %s)!\n", cudaGetErrorString(err));
                    }
                    printf(" ... done! \n ");


                    printf("Piece together the calculated results into full output matrix.\n");
                    if ( blk_size != max_blk_size ) {
                        range_n = blk_size;
                    } 
                    else {
                        range_n = max_blk_size;
                    }
                    for (int bi = 0; bi < range_n; ++bi) 
                    {
                       for (int bj = bi+1; bj < range_n; ++bj) 
                       {
                            bk = (range_n*(range_n-1)/2)-(range_n-bi)*((range_n-bi)-1)/2+bj-bi-1;

                            Fi = (blk_x*max_blk_size)+bi;
                            Fj = (blk_y*max_blk_size)+bj;
                            Fk = (gV*(gV-1)/2)-(gV-Fi)*((gV-Fi)-1)/2+Fj-Fi-1;

                            output_buffer[Fk] = gpu_ccmat[bk];
                       }
                    }
                    printf("...done!\n");

                    //printf("saving block to txt file ...\n"); 
                    //if (current_run < 10) { 
                    //    sprintf(s, "dbgblk_0%d.txt", current_run);
                    //} else {
                    //    sprintf(s, "dbgblk_%d.txt", current_run);
                    //}
                    //printf("to file: %s", s);
                    //saveToText_TRIA(gpu_ccmat, blk_size, s);
                    //printf("...done!\n");

                    printf("Freeing block-temporary memory\n");
	                err = cudaFree(d_data);
	                err = cudaFree(d_result);
	                if ( err != cudaSuccess ) 
	                {
                        fprintf(stderr, "Failed to free device data (error code %s)!\n", cudaGetErrorString(err));
		                exit(EXIT_FAILURE);
	                }
                    free(h_data);
                    free(gpu_ccmat);
                    printf("...done!\n");
                }
                else 
                {
                    printf("Current run number: %d\n", current_run);
                    printf("Block indices (blk_x, blk_y) = (%d, %d)\n", blk_x, blk_y);

                    blk_size_x = (end_idx[blk_x]-start_idx[blk_x]+1);
                    blk_size_y = (end_idx[blk_y]-start_idx[blk_y]+1);
                    printf("This is a off-diagonal block with block size (%zu, %zu)\n", blk_size_x, blk_size_y);
                    input_size = (blk_size_x+blk_size_y)*nt_data*sizeof(float);
                    output_size = (blk_size_x*blk_size_y)*sizeof(float);

                    /* allocate memory on host */
                    h_data = (float *)malloc(input_size);
                    gpu_ccmat = (float *)malloc(output_size);
                    printf("Parsing input data into appropriate format on host\n");
                    inputData4DtoArrayWith2Ranges(data, h_data, mask, nx_data, ny_data, nz_data, nt_data, 
                                                  start_idx[blk_y], end_idx[blk_y], 
                                                  start_idx[blk_x], end_idx[blk_x]);

                    if ( h_data == NULL || gpu_ccmat == NULL ) 
                    {
                        fprintf(stderr, "Failed to allocate memory on host!\n");
                        return EXIT_FAILURE;
                    }
                    printf("...done!\n");

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
                    /* running the cuda kernel */
                    printf("Initializing CUDA (OFFD) kernel ...\n");
                    CUDA_CC_OFFD<<<dimGrid, dimBlock>>>(d_data, d_result, blk_size_x+blk_size_y, nt_data, blk_size_x);

                    /* copy the results form the device back to host */
                    err = cudaMemcpy(gpu_ccmat, d_result, output_size, cudaMemcpyDeviceToHost);
                    if ( err != cudaSuccess ) 
                    {
                        printf("Failed to copy result from device to host (error code %s)!\n", cudaGetErrorString(err));
                    }
                    printf(" ... done! \n ");


                    printf("Piece together the calculated results into full output matrix.\n");
                    if ( blk_size_x != blk_size_y ) {
                        range_x = blk_size_x;
                        range_y = blk_size_y;
                    } 
                    else {
                        range_x = max_blk_size;
                        range_y = max_blk_size;
                    }
                    for (int bi = 0; bi < blk_size_x; ++bi) 
                    {
                        for (int bj = 0; bj < blk_size_y; ++bj)
                        {
                            Fi = (blk_x*max_blk_size)+bi;
                            Fj = (blk_y*max_blk_size)+bj;
		                    Fk = (gV*(gV-1)/2)-(gV-Fi)*((gV-Fi)-1)/2+Fj-Fi-1;

                            bk = (bi*range_y)+bj;
                            output_buffer[Fk] = gpu_ccmat[bk];
                        }
                    }
                    printf("...done!\n");

                    //printf("saving block to txt file ...\n"); 
                    //if (current_run < 10) {
                    //    sprintf(s, "testfile_0%d.txt", current_run);
                    //} else {
                    //    sprintf(s, "testfile_%d.txt", current_run);
                    //}
                    //printf("to file: %s", s);
                    //saveToText_RECT(gpu_ccmat, blk_size_x, blk_size_y, s);
                    //printf("...done!\n");

                    printf("Freeing block-temporary memory\n");
	                err = cudaFree(d_data);
	                err = cudaFree(d_result);
	                if ( err != cudaSuccess ) 
	                {
                        fprintf(stderr, "Failed to free device data (error code %s)!\n", cudaGetErrorString(err));
		                exit(EXIT_FAILURE);
	                }
                    free(h_data);
                    free(gpu_ccmat);
                    printf("...done!\n");
                }
                current_run++;
            }
        }
        /* end clock and report time */
        diff = clock() - start_gpu;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    }
    else
    {
	    printf("allocate memory on host ...\n");
	    float *h_data = (float *)malloc(input_size);
	    if ( h_data == NULL ) 
        {
		    fprintf(stderr, "Failed to allocate memory on host!\n");
		    exit(EXIT_FAILURE);
	    }

        printf("Parsing input data into appropriate format on host ...\n");
        inputData4DtoArray(data, h_data, mask, nx_data, ny_data, nz_data, nt_data);
        printf("... immediately after, free up a little memory by freeing the 4D data ...\n");
        free(data);
        printf("...done!\n");


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
        clock_t start_gpu = clock(), diff; // start the clock!
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

        /* running the cuda kernel */
        printf("Initializing CUDA kernel ...\n");
        CUDA_CC_DIAG<<<dimGrid, dimBlock>>>(d_data, d_result, gV, nt_data);

        /* copy the results form the device back to host */
        err = cudaMemcpy(output_buffer, d_result, output_size, cudaMemcpyDeviceToHost);
        if ( err != cudaSuccess ) 
        {
            printf("Failed to copy result from device to host (error code %s)!\n", cudaGetErrorString(err));
        }
        printf("..done!\n");

        /* end clock and report time */
        diff = clock() - start_gpu;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

        /* free things we have malloced: 
           h_data, d_data, and d_result */
        free(h_data);
        err = cudaFree(d_data);
        err = cudaFree(d_result);
        if ( err != cudaSuccess ) 
        {
            fprintf(stderr, "Failed to free device data (error code %s)!\n", cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }
    }
	
	/* -------------------------------------------------------
	 * saving output
	 * -------------------------------------------------------
	 */

    clock_t start_save = clock(), diff_save; // start the clock!
    printf("Saving output ...  \n");
    //float **pmat = createMatrix(gV, gV);
    //upperTriangularToFull(cpu_ccmat, pmat, gV);
    saveToText_TRIA(output_buffer, gV, outputFile);
    printf(" ... done! \n ");
    diff_save = clock() - start_save;
    int msec_save = diff_save * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds\n", msec_save/1000, msec_save%1000);

    free(output_buffer);
	return 0;
}

