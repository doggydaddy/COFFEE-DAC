#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <omp.h>

///
/// @brief returns number of lines a file has
/// @param filename input txt file
/// @return number of lines 
/// 
/// used to get number of voxels (not connections!) a subject data file has.
///
size_t getNumberLines(char* filename)
{
    FILE *stream;
    char *line = NULL;
    size_t len = 0;
    size_t nlines = 0;

    stream = fopen(filename, "r");
    if (stream == NULL) 
    {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    while( getline(&line, &len, stream) != -1 ) 
    {
        nlines++;
    }

    fclose(stream);
    free(line);

    return(nlines);
}

///
/// @brief returns the number of values a file has in the first line
/// @param filename input txt file
/// @return number of values in the first line
///
/// used to get number of voxels (not connections!) a subject data file has.
/// Sanity check getNumberLines() too as proper subject connectivity data should
/// have the same number of lines as values in the first line!
///
size_t getNumberValsFirstLine(char* filename)
{
    FILE *stream;
    char *line = NULL;
    char *linebuff;
    size_t len = 0;

    stream = fopen(filename, "r");
    if (stream == NULL) 
    {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    size_t line_length = 0;
    while( getline(&line, &len, stream) != -1 ) 
    {
        line_length = 0;
        linebuff = strtok(line, " ");
        while ( linebuff != NULL ) 
        {
            linebuff = strtok(NULL, " ");
            line_length++;
        }
        break;
    }

    fclose(stream);
    free(line);

    return(line_length);
}

///
/// @brief peek into file list and grabs the number of lines each subject has
/// @param filelist input file list 
/// @return number of lines each subject has
///
/// convenience routine to grab subject dimensions using the file list
/// references only
///
size_t peekFileList(char* filelist)
{
    FILE *fl;
    char* fl_line = NULL;
    size_t fl_len = 0;

    size_t output;

    fl = fopen(filelist, "r");
    if (fl == NULL) 
    {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    while( getline(&fl_line, &fl_len, fl) != -1 )
    {
        fl_line[strcspn(fl_line, "\n")] = 0;
        output = getNumberLines(fl_line);
    } 

    fclose(fl);
    free(fl_line);

    return(output);
}

///
/// @brief reads permutations file
/// @param filename generated permutations file from generatePermutations.py program
/// @param buffer array of ints to store the parsed permutations
/// @param nr_subs number of subjects in total for the test
///
/// note the expected input (as obtained from the output of
/// generatePermutations.py) is indices of a group (group A), and NOT one-hot
///
/// the output buffer contains one-hot labels of the subject permutations
///
void parsePermutations(char* filename, int* buffer, size_t nr_subs)
{
    /* variables to load permutations file */
    FILE *pt;
    char *pt_line = NULL;
    size_t pt_len = 0;
    int pt_lines;

    /* variables to parse file */
    char* line_buff;
    float val_buff;
    size_t line_idx;

    // opening permutations file
    pt = fopen(filename, "r");
    if (pt == NULL) 
    {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    line_idx = 0;
    while ((pt_lines = getline(&pt_line, &pt_len, pt)) != -1) 
    {
        if (line_idx == 0)
        {
            line_buff = (char*)malloc(sizeof(char)*pt_lines);
        }

        line_buff = strtok(pt_line, " ");
        while ( line_buff != NULL ) 
        {
            if (*line_buff == '\n') 
            {
                line_buff = strtok(NULL, " ");
                continue;
            } 
            else
            {
                val_buff = atof(line_buff);
                buffer[(int)((line_idx*nr_subs)+val_buff)] = 1;
                line_buff = strtok(NULL, " ");
            }
        }
        line_idx++;
    } 

    // cleanup
    free(line_buff);
}

///
/// @brief reads file list and parses all subjects from index N to M
/// @param filename input file list
/// @param nr_sub number of subjects in total in file list
/// @param N from (and including) index 
/// @param M to (and including) index
/// @param buffer output buffer of subject connection values
///
/// Output buffer contains index N as first row to index M as last row of
/// values. For each row subject order is the same as specified in the file
/// list.
///
void parseFileListNtoM(char* filename, int nr_sub, int N, int M, float* buffer)
{
    /* vars */
    /* variables to load file list */
    FILE *fl;
    char* fl_line = NULL;
    size_t fl_len = 0;

    /* variables to load file */
    FILE *stream;
    char *line = NULL;
    size_t len = 0;
    int nread;

    /* variables to parse file */
    char* line_buff;
    float buff;
    size_t line_idx;
    size_t k;
    size_t line_length;

    /* variables for device buffer */
    size_t row_counter = 0;
    /* /vars */

    fl = fopen(filename, "r");
    if (fl == NULL) 
    {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    size_t sub_idx = 0;
    while( getline(&fl_line, &fl_len, fl) != -1 )
    {
        /* get filename (remove trailing \n)*/
        fl_line[strcspn(fl_line, "\n")] = 0;

        /* load file */
        stream = fopen(fl_line, "r");
        if (stream == NULL) 
        {
            perror("fopen");
            exit(EXIT_FAILURE);
        }

        /* read file line by line */
        line_idx = 0;
        k = 0;
        row_counter = 0;
        while ((nread = getline(&line, &len, stream)) != -1) 
        {
            if (line_idx == 0) // firstline
            {
                line_buff = (char*)malloc(sizeof(char)*nread);
            }

            line_length = 0;
            line_buff = strtok(line, " ");
            while ( line_buff != NULL ) 
            {
                if (*line_buff == '\n') 
                {
                    line_buff = strtok(NULL, " ");
                    continue;
                } 
                else
                {
                    buff = atof(line_buff);
                    if (k >= N && k <= M)
                    {
                        buffer[(row_counter*nr_sub)+sub_idx] = buff;
                        row_counter++;
                    }
                    ++k;
                    line_length++;
                    line_buff = strtok(NULL, " ");
                }
            }

            line_idx++;
        } // end of reading a file line by line

        sub_idx++;
    } // end of parsing file list  

    /* cleanup */
    fclose(fl);
    free(fl_line);
    
    fclose(stream);
    free(line);

    free(line_buff);

    printf("done!\n");
}

///
/// @brief save results as upper triangular format
/// @param outputData data to be saved to file
/// @param nrows number of rows in the output file/nr voxels the output should have
/// @param fileName output file name
///
/// note that nrows specifies the number of voxels, not connections each subject
/// has.
///
void saveResToText(float *outputData, size_t nrows, char *fileName)
{
    FILE *output = fopen(fileName, "w");
    
    size_t c = 0;
    for (size_t i=0; i<nrows; ++i)
    {
        for (size_t j=i+1; j<nrows; ++j) 
        {
            fprintf(output, "%f ", outputData[c]);
            c++;
        }
        fprintf(output, "\n");
    }

    printf("[DBG] saved a total of %zu values\n", c);
    fclose(output);
}

///
/// @brief CUDA kernel to permform permutation test
/// @param input input buffer (subject values)
/// @param onehot permutations buffer in one-hot format
/// @param nr_vals number of values to process (connections)
/// @param nr_sub number of subjects in total
/// @param nr_perm number of permutations
/// @param output output buffer
/// @param two_tailed is two-tailed test or not (1 = two-tailed, 1 = single-tailed)
///
/// note GPU memory usage is NOT optimized!
///
__global__ 
void CUDA_perm(float *input, int* onehot, 
               size_t nr_vals, size_t nr_sub, size_t nr_perm, 
               float *output, int two_tailed)
{
    int n = blockIdx.x;
    if (n > nr_vals) return;

    float t_obs = 0.;
    float p_val = 0.;

    float a_mean = 0.;
    float b_mean = 0.;
    float nA = 0;
    float nB = 0;
    float tstat = 0.0;

    for (int i=0; i<nr_perm; ++i) 
    {
        // t-stat
        a_mean = 0.;
        b_mean = 0.;
        nA = 0;
        nB = 0;
        tstat = 0;
        for (int j=0; j<nr_sub; ++j)
        {
            if (onehot[(i*nr_sub)+j] == 0) // i-th row in perm
            {
                b_mean += input[(n*nr_sub)+j];
                nB++;
            } 
            else 
            {
                a_mean += input[(n*nr_sub)+j]; // n-th row in subject data
                nA++;
            }
        }
        a_mean /= nA;
        b_mean /= nB;
        if (two_tailed == 1) 
        {
            tstat = fabs(a_mean - b_mean);
        }
        else
        {
            tstat = a_mean - b_mean;
        }
        // /t-stat

        if (i == 0) // first permutation, t_obs is tstat
        {
            t_obs = tstat;
        }

        if (tstat > t_obs) 
        {
            p_val++;
        }

    }
    p_val = p_val/(float)(nr_perm);
    output[n] = p_val;
}

int
main(int argc, char *argv[])
{
    /* read arguments */    
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <file list> <permutations file> <output file>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* parse input arguments */
    char* filelist = argv[1];
    char* permutations = argv[2];
    char* outfile = argv[3];

    /* get dimensions */
    ssize_t nr_r1vals = peekFileList(filelist);
    ssize_t nr_vals = (nr_r1vals*(nr_r1vals-1))/2;
    size_t nr_perm = getNumberLines(permutations);
    size_t nr_subs = getNumberLines(filelist);
    printf("Number of permutations: %zd\n", nr_perm);
    printf("Number of subjects: %zu\n", nr_subs);
    printf("Number of connections in each subject: %zu\n", nr_vals);

    /* calculate how many voxels values we can test at once with our available memory */
    size_t device_free_mem, device_total_mem;
	cudaMemGetInfo(&device_free_mem, &device_total_mem);
	printf("GPU free mem: %lu, Total mem: %lu\n", device_free_mem, device_total_mem);
    size_t nr_vals_max = 0.9*((device_free_mem-(sizeof(int)*nr_perm*nr_subs))/(sizeof(float)*(nr_subs+1)));
    printf("Number of values we can load at once into device memory: %zu\n", nr_vals_max);

    if (nr_vals <= nr_vals_max) 
    {
        printf("can load the entire thing in one go!\n");

        printf("allocating memory buffers ...\n");
        int* perm_buff = (int*)malloc(sizeof(int)*nr_perm*nr_subs);
        /* we have to zero out the permutation buffer */
        for (int i=0; i<nr_perm*nr_subs; ++i)
        {
            perm_buff[i] = 0; 
        }
        float* device_buff = (float*)malloc(sizeof(float)*nr_vals*nr_subs);
        
        printf("parsing input files ...\n");
        parsePermutations(permutations, perm_buff, nr_subs);
        parseFileListNtoM(filelist, nr_subs, 0, nr_vals-1, device_buff);

        float* perm_test_res;
        perm_test_res = (float*)malloc(sizeof(float)*nr_vals);

        printf("performing permutation tests ...\n");
        ///* allocating gpu mem */
        int *d_perm;
        float *d_input;
        float *d_output;
        cudaError_t err = cudaSuccess;;
        err = cudaMalloc((void **)&d_perm, sizeof(int)*nr_perm*nr_subs);
        err = cudaMalloc((void **)&d_input, sizeof(float)*nr_vals*nr_subs);
        err = cudaMalloc((void **)&d_output, sizeof(float)*nr_vals );
        if (err!=cudaSuccess)
        {
            printf("CUDA ERROR! Failed to allocate device memory! (error code %s)\n", cudaGetErrorString(err));
        }
        /* copy input data host -> device */
        err = cudaMemcpy(d_input, device_buff, sizeof(float)*nr_vals*nr_subs, cudaMemcpyHostToDevice);
        err = cudaMemcpy(d_perm, perm_buff, sizeof(int)*nr_perm*nr_subs, cudaMemcpyHostToDevice);
        if (err!=cudaSuccess)
        {
            printf("CUDA ERROR! Failed to copy memory from host to device (error code %s)\n", cudaGetErrorString(err));
        }

        /* setting up and perform calculations on gpu */
        CUDA_perm<<<nr_vals, 1>>>(d_input, d_perm, nr_vals, nr_subs, nr_perm, d_output, 0);
        /* copy results from the device back to host */
        err = cudaMemcpy(perm_test_res, d_output, sizeof(float)*nr_vals, cudaMemcpyDeviceToHost);
        if (err!=cudaSuccess)
        {
            printf("CUDA ERROR! Failed to copy memory from device to host (error code %s)\n", cudaGetErrorString(err));
        }

        printf("writing to file ...\n");
        saveResToText(perm_test_res, nr_r1vals, outfile);
        printf("done!\n");

        /* cleanup */
        free(device_buff);
        free(perm_buff);
        free(perm_test_res);

        err = cudaFree(d_input);
        err = cudaFree(d_output);
        err = cudaFree(d_perm);
        if (err!=cudaSuccess)
        {
            printf("CUDA ERROR! Failed to free device memory! (error code %s)\n", cudaGetErrorString(err));
        }
    }
    else 
    {
        printf("we cannot load the entire thing in one go ...\n");
        /* how many parts do we have to split stuff in? */
        int nr_parts = ceil(nr_vals/nr_vals_max);
        printf("... so we have to split the job into %i parts\n", nr_parts);

        /* calculate indices and sizes */
        size_t part_starts[nr_parts];
        size_t part_ends[nr_parts];
        size_t part_vals[nr_parts];
        for (int p=0; p<nr_parts; ++p)
        {
            if (p==0)
            {
                part_vals[p] = int(nr_vals/nr_parts);
                part_starts[p] = 0;
                part_ends[p] = part_vals[0]-1;
            }
            else 
            {
                part_starts[p] = part_ends[p-1]+1;
                part_vals[p] = part_vals[p-1];
                part_ends[p] = part_starts[p]+part_vals[p]-1;
            }    
        }
        int diffset;
        if (part_ends[nr_parts-1] < nr_vals-1) 
        {   /* correcting up */
            diffset = (nr_vals-1) - part_ends[nr_parts-1];
            part_ends[nr_parts-1] = nr_vals-1;
            part_vals[nr_parts-1] += diffset;
        } 
        else if (part_ends[nr_parts-1] > nr_vals-1) 
        {   /* correcting down */
            diffset = part_ends[nr_parts-1] - (nr_vals-1);
            part_ends[nr_parts-1] = nr_vals-1;
            part_vals[nr_parts-1] -= diffset;
        }

        /* dbg printout part divisions */ 
        printf("[DBG] part division indices:\n");
        for (int i=0; i<nr_parts; ++i) 
        {
            printf("[%zu, %zu] (%zu)\n", part_starts[i], part_ends[i], part_vals[i]);
        }

        /* allocate full result buffer */
        float* full_results = (float*)malloc(sizeof(float)*nr_vals);

        /* performing calculation in parts */

        // host buffers 
        int* perm_buff;
        float* device_buff;
        float* perm_test_res;

        // device buffers
        int *d_perm;
        float *d_input;
        float *d_output;

        // counters
        cudaError_t err = cudaSuccess;;
        size_t id = 0;

        for (int p=0; p<nr_parts; ++p)
        {
            printf("performing part %i of %i calculations ...\n", p, nr_parts);

            printf("allocating memory buffers ...\n");
            perm_buff = (int*)malloc(sizeof(int)*nr_perm*nr_subs);
            // we have to zero out the permutation buffer
            for (int i=0; i<nr_perm*nr_subs; ++i)
            {
                perm_buff[i] = 0; 
            }
            device_buff = (float*)malloc(sizeof(float)*part_vals[p]*nr_subs);
            
            printf("parsing input files ...\n");
            parsePermutations(permutations, perm_buff, nr_subs);
            parseFileListNtoM(filelist, nr_subs, part_starts[p], part_ends[p], device_buff);

            perm_test_res = (float*)malloc(sizeof(float)*part_vals[p]);

            printf("performing permutation tests ...\n");
            /* allocating gpu mem */
            err = cudaMalloc((void **)&d_perm, sizeof(int)*nr_perm*nr_subs);
            err = cudaMalloc((void **)&d_input, sizeof(float)*part_vals[p]*nr_subs);
            err = cudaMalloc((void **)&d_output, sizeof(float)*part_vals[p] );
            if (err!=cudaSuccess)
            {
                printf("CUDA ERROR! Failed to allocate device memory! (error code %s)\n", cudaGetErrorString(err));
            }
            /* copy input data host -> device */
            err = cudaMemcpy(d_input, device_buff, sizeof(float)*part_vals[p]*nr_subs, cudaMemcpyHostToDevice);
            err = cudaMemcpy(d_perm, perm_buff, sizeof(int)*nr_perm*nr_subs, cudaMemcpyHostToDevice);
            if (err!=cudaSuccess)
            {
                printf("CUDA ERROR! Failed to copy memory from host to device (error code %s)\n", cudaGetErrorString(err));
            }

            /* setting up and perform calculations on gpu */
            CUDA_perm<<<part_vals[p], 1>>>(d_input, d_perm, part_vals[p], nr_subs, nr_perm, d_output, 0);
            /* copy results from the device back to host */
            err = cudaMemcpy(perm_test_res, d_output, sizeof(float)*part_vals[p], cudaMemcpyDeviceToHost);
            if (err!=cudaSuccess)
            {
                printf("CUDA ERROR! Failed to copy memory from device to host (error code %s)\n", cudaGetErrorString(err));
            }

            /* writing part results to full result buffer */
            printf("writing to full results ...\n");
            id = 0;
            for (size_t c=part_starts[p]; c<=part_ends[p]; ++c)
            {
                full_results[c] = perm_test_res[id];
                id++;
            }

            // cleanup
            free(device_buff);
            free(perm_buff);
            free(perm_test_res);

            err = cudaFree(d_input);
            err = cudaFree(d_output);
            err = cudaFree(d_perm);
            if (err!=cudaSuccess)
            {
                printf("CUDA ERROR! Failed to free device memory! (error code %s)\n", cudaGetErrorString(err));
            }

        }

        printf("writing to file ...\n");
        saveResToText(full_results, nr_r1vals, outfile);
        printf("done!\n");

        free(full_results);
    }

    return(EXIT_SUCCESS);
}