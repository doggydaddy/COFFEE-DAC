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

// for openMP parallel computing
#include <omp.h>

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
            // linear indexing here will not work, because i and j is not linear!
            // if we keep calcCrossCorr_omp as-is, we need to manipulate the i,j indices, specifically index j 
            //if ( (i*cols)+j >= rows*cols ) 
            //{
            //    printf("result index (%zu, %zu) = (%zu) exceeded size of result (%zu) array\n", i, j, (size_t)(i*cols)+j, (size_t)rows*cols);
            //    
            //}
            //if ( j+rows > gV ) 
            //{
            //    printf("index j exceeds size size of data array\n");
            //}
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

    /* buffer output image */ 
    printf("allocate buffer memory on host ...\n");
    float *output_buffer = (float*)malloc(output_size);
    if ( output_buffer == NULL ) 
    {
        fprintf(stderr, "Failed to allocate output buffer on host memory!\n");
        printf("The job is too big to fit onto the system (host) memory,\n");
        printf("please upgrade your potato computer and try again!\n");
        exit(EXIT_FAILURE);
    }

    printf("See how much memory we have: \n");
	unsigned long long system_memory = getTotalSystemMemory();
    // we required (about) input_size*3 plus the output_size to be able to
    // compute everything in one go 
    size_t required_device_memory = (size_t)(input_size+(output_size*2));
    printf("We actually have %llu bytes ...\n", system_memory);
    // 4GB virtual memory cap, with 90% of that being the practical amount we can use
    size_t virtual_device_memory = (size_t)(1*up3zeros*up3zeros*up3zeros);
    printf("... but let us pretend to have only %zu.\n", virtual_device_memory);

    float num_divisions = 1.0;
    float num_runs = 1.0;
	if ( required_device_memory > virtual_device_memory ) 
	{
		printf("The job is too big to fit onto the device memory,\n");

        size_t index_cap = (size_t)((sqrtf(virtual_device_memory-input_size))/(2*sizeof(float)));
        printf("You can process at most %zu indices at a time.\n", index_cap);

        num_divisions = ceil((float)gV/(float)index_cap);
        printf("so you need at least to make %d divisions.\n", (int)num_divisions);
        int start_idx[(int)num_divisions];
        int end_idx[(int)num_divisions];
        for (int i=0; i<num_divisions; ++i)
        {
            start_idx[i] = i*index_cap;
            end_idx[i] = ((i*index_cap)+index_cap)-1;
            if (end_idx[i] >= gV) 
            {
                end_idx[i] = gV - 1;
                continue;
            } 
        }
        printf("DBG: printing out the indices: \n");
        for (int i=0; i<num_divisions; ++i)
        {
            printf("division %d, [%d, %d]\n", i, start_idx[i], end_idx[i]);
        }

        num_runs = (num_divisions*(num_divisions+1))/2;
        printf("so you need at least %d runs to complete the entire job.\n", (int)num_runs);

        size_t run_output_size = index_cap*index_cap*sizeof(float); // device can only allocate so much memory at a time
        output_size = run_output_size;


        /* performing the block runs*/
        int current_run = 0; // for naming purposes
        float *h_data;
        float *cpu_ccmat;
        size_t Fi, Fj, Fk, bk;
        size_t blk_size_x, blk_size_y, blk_size;
        char s[26];
        size_t range_n, range_x, range_y;
        for (int blk_x = 0; blk_x < num_divisions; ++blk_x) 
        {
            for (int blk_y = blk_x; blk_y < num_divisions; ++blk_y) 
            {
                if (blk_x == blk_y) 
                {
                    printf("Current run number: %d\n", current_run);
                    printf("Block indices (blk_x, blk_y) = (%d, %d)\n", blk_x, blk_y);
                    printf("this is a diagonal block, run diagonal code.\n");

                    printf("allocate memory on host ...\n");
                    blk_size = end_idx[blk_x] - start_idx[blk_x] + 1;

                    input_size = blk_size*nt_data*sizeof(float);
                    output_size = ((blk_size*(blk_size-1))/2)*sizeof(float);
                    h_data = (float *)malloc(input_size);
                    cpu_ccmat = (float *)malloc(output_size);

                    printf("Parsing input data into appropriate format on host\n");
                    inputData4DtoArrayWithRanges(data, h_data, mask, 
                                                 nx_data, ny_data, nz_data, nt_data, 
                                                 start_idx[blk_x], end_idx[blk_x]);
                    printf("...done!\n");

                    /* double-check if we succeeded or not */
                    if ( h_data == NULL || cpu_ccmat == NULL ) 
                    {
                        fprintf(stderr, "Failed to allocate memory on host!\n");
                        exit(EXIT_FAILURE);
                    }

                    printf("Performing calcCrossCorr_DIAG\n");
                    calcCrossCorr_DIAG(h_data, cpu_ccmat, blk_size, nt_data);
                    printf("...done!\n");
                    printf("sanity check cpu_ccmat:\n");
                    if ( sanityCheck(cpu_ccmat, (size_t)(output_size/sizeof(float))) == 0 ) {
                        printf("sanity check failed! exiting \n");
                        return -1;
                    }
                    printf("...done!\n");

                    printf("Piece together the calculated results into full output matrix.\n");
                    if ( blk_size != index_cap ) {
                        range_n = blk_size;
                    } 
                    else {
                        range_n = index_cap;
                    }
                    for (int bi = 0; bi < range_n; ++bi) 
                    {
                       for (int bj = bi+1; bj < range_n; ++bj) 
                       {
                            bk = (range_n*(range_n-1)/2)-(range_n-bi)*((range_n-bi)-1)/2+bj-bi-1;

                            Fi = (blk_x*index_cap)+bi;
                            Fj = (blk_y*index_cap)+bj;
                            Fk = (gV*(gV-1)/2)-(gV-Fi)*((gV-Fi)-1)/2+Fj-Fi-1;

                            output_buffer[Fk] = cpu_ccmat[bk];
                       }
                    }
                    printf("...done!\n");

                    printf("saving block to txt file ...\n"); 
                    if (current_run < 10) {
                        sprintf(s, "testfile_0%d.txt", current_run);
                    } else {
                        sprintf(s, "testfile_%d.txt", current_run);
                    }
                    printf("to file: %s", s);
                    saveToText_TRIA(cpu_ccmat, blk_size, s);
                    printf("...done!\n");

                    printf("Freeing block-temporary memory\n");
                    free(h_data);
                    free(cpu_ccmat);
                    printf("...done!\n");
                }
                else 
                {
                    printf("Current run number: %d\n", current_run);
                    printf("Block indices (blk_x, blk_y) = (%d, %d)\n", blk_x, blk_y);
                    printf("this is a off-diagonal block, run off-diagonal code.\n");

                    printf("allocate memory on host ...\n");
                    blk_size_x = (end_idx[blk_x]-start_idx[blk_x]+1);
                    blk_size_y = (end_idx[blk_y]-start_idx[blk_y]+1);
                    input_size = (blk_size_x+blk_size_y)*nt_data*sizeof(float);
                    output_size = (blk_size_x*blk_size_y)*sizeof(float);
                    h_data = (float *)malloc(input_size);
                    cpu_ccmat = (float *)malloc(output_size);

                    printf("Parsing input data into appropriate format on host\n");
                    inputData4DtoArrayWith2Ranges(data, h_data, mask, nx_data, ny_data, nz_data, nt_data, 
                                                  start_idx[blk_y], end_idx[blk_y], 
                                                  start_idx[blk_x], end_idx[blk_x]);
                    printf("...done!\n");

                    printf("Performing calcCrossCorr_OFFD\n");
                    calcCrossCorr_OFFD(h_data, cpu_ccmat, blk_size_x+blk_size_y, nt_data, blk_size_x);
                    printf("...done!\n");
                    printf("sanity check cpu_ccmat:\n");
                    if ( sanityCheck(cpu_ccmat, (size_t)(output_size/sizeof(float))) == 0 ) {
                        printf("sanity check failed! exiting \n");
                        return -1;
                    }
                    printf("...done!\n");

                    printf("Piece together the calculated results into full output matrix.\n");
                    if ( blk_size_x != blk_size_y ) {
                        range_x = blk_size_x;
                        range_y = blk_size_y;
                    } 
                    else {
                        range_x = index_cap;
                        range_y = index_cap;
                    }
                    for (int bi = 0; bi < blk_size_x; ++bi) 
                    {
                        for (int bj = 0; bj < blk_size_y; ++bj)
                        {
                            Fi = (blk_x*index_cap)+bi;
                            Fj = (blk_y*index_cap)+bj;
		                    Fk = (gV*(gV-1)/2)-(gV-Fi)*((gV-Fi)-1)/2+Fj-Fi-1;

                            bk = (bi*range_y)+bj;
                            output_buffer[Fk] = cpu_ccmat[bk];
                        }
                    }
                    printf("...done!\n");

                    printf("saving block to txt file ...\n"); 
                    if (current_run < 10) {
                        sprintf(s, "testfile_0%d.txt", current_run);
                    } else {
                        sprintf(s, "testfile_%d.txt", current_run);
                    }
                    printf("to file: %s", s);
                    saveToText_RECT(cpu_ccmat, blk_size_x, blk_size_y, s);
                    printf("...done!\n");

                    printf("Freeing block-temporary memory\n");
                    free(h_data);
                    free(cpu_ccmat);
                    printf("...done!\n");
                }
                current_run++;
            }
        }
    }
    else
    {
	    printf("allocate memory on host ...\n");
	    float *h_data = (float *)malloc(input_size);

        printf("Parsing input data into appropriate format on host\n");
        inputData4DtoArray(data, h_data, mask, nx_data, ny_data, nz_data, nt_data);
        printf("...done!\n");
        printf("immediately after, free up the little memory we have by freeing the 4D data\n");
        free(data);
        printf("...done!\n");

        /* double-check if we succeeded or not */
        if ( h_data == NULL ) 
        {
            fprintf(stderr, "Failed to allocate memory on host!\n");
            exit(EXIT_FAILURE);
        }

        /* performing the run */
        printf("Performing calcCrossCorr_DIAG\n");
        calcCrossCorr_DIAG(h_data, output_buffer, gV, nt_data);
        printf("...done!\n");

        free(h_data);
    }
	

	/* -------------------------------------------------------
	 * saving output
	 * -------------------------------------------------------
	 */
    printf("Saving output ...  \n");
    //float **pmat = createMatrix(gV, gV);
    //upperTriangularToFull(cpu_ccmat, pmat, gV);
    saveToText_TRIA(output_buffer, gV, outputFile);
    printf(" ... done! \n ");

    free(output_buffer);
	return 0;
}

