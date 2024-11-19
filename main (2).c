#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include <string.h>
#include <time.h>

int main(int argc, char *argv[])
{
    
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    if (argc < 2 || argc > 5 || argc == 3)
    {
        fprintf(stderr, "Please use the correct format\n");
        return 1;
    }

    const char *filename = argv[1];
    printf("Trying to open file: %s\n", filename);

    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    if (argc == 2)
    {
        print_CSR_Matrix(&A);
        end_time = clock();
        cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
        printf("CPU time: %f seconds\n", cpu_time_used);
        freeCSR(&A);
        return 0;
    }
    else if (argc == 4 && (strcmp(argv[2], "transpose") == 0))
    {
        CSRMatrix T = transposeCSR(&A);
        if (atoi(argv[3]) == 1)
        {
            printf("Matrix from %s\n", filename);
            print_CSR_Matrix(&A);
            printf("\n");
            printf("Transpose of matrix from %s\n", filename);
            print_CSR_Matrix(&T);
            printf("\n");
        }
        freeCSR(&A);
        freeCSR(&T);
        end_time = clock();
        cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
        printf("CPU time: %f seconds\n", cpu_time_used);
        return 0;
    }

    const char *filename_2 = argv[2];
    CSRMatrix B;
    ReadMMtoCSR(filename_2, &B);
    CSRMatrix C;
    const char *calc = argv[3];

    if (argc == 5)
    {
        if (strcmp(calc, "addition") == 0)
        {
            C = addCSR(&A, &B);
        }
        else if (strcmp(calc, "subtract") == 0)
        {
            C = subtractCSR(&A, &B);
        }
        else if (strcmp(calc, "multiply") == 0)
        {
            C = multiplyCSR(&A, &B);
        }
        else
        {
            fprintf(stderr, "Please use one of the following when calculating: addition, subtract, multiply or transpose. \n");
            freeCSR(&A);
            freeCSR(&B);
            end_time = clock();
            cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
            printf("CPU time: %f seconds\n", cpu_time_used);
            return 1;
        }
    }

    if (atoi(argv[4]) == 1)
    {
        printf("Matrix A:\n");
        print_CSR_Matrix(&A);
        printf("\n");
        printf("Matrix B:\n");
        print_CSR_Matrix(&B);
        printf("\n");
        printf("Resultant Matrix C:\n");
        print_CSR_Matrix(&C);
        printf("\n");
    }

    freeCSR(&A);
    freeCSR(&B);
    freeCSR(&C);

    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time: %f seconds\n", cpu_time_used);

    return 0;
}


