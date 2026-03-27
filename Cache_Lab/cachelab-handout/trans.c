/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 * 
 * the cache has parameters that s = 5, E = 1, b = 5
 * 
 * evaluate matrix size:
 *  1. M = 32, N = 32 ; m<300  will be marked as correct
 *  2. M = 64, N = 64 ; m<1300 will be marked as correct
 *  3. M = 61, N = 67 ; m<2000 will be marked as correct
 * Ok to explicitly check the input size and implement different transpose strategies for different cases
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 * 
 * 最终miss分别为288,1228,1929
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int a0 , a1 , a2 , a3 , a4 , a5 , a6 , a7 , i , j , k ;
    if(M==32)
    {
        for(i=0 ; i<4 ; i++)
        for(j=0 ; j<4 ; j++)
            for(k=0 ; k<8 ; k++)
            {
                int a0 = A[i*8+k][j*8+0] ;
                int a1 = A[i*8+k][j*8+1] ;
                int a2 = A[i*8+k][j*8+2] ;
                int a3 = A[i*8+k][j*8+3] ;
                int a4 = A[i*8+k][j*8+4] ;
                int a5 = A[i*8+k][j*8+5] ;
                int a6 = A[i*8+k][j*8+6] ;
                int a7 = A[i*8+k][j*8+7] ;
                B[j*8+0][i*8+k] = a0 ;
                B[j*8+1][i*8+k] = a1 ;
                B[j*8+2][i*8+k] = a2 ;
                B[j*8+3][i*8+k] = a3 ;
                B[j*8+4][i*8+k] = a4 ;
                B[j*8+5][i*8+k] = a5 ;
                B[j*8+6][i*8+k] = a6 ;
                B[j*8+7][i*8+k] = a7 ;
            }
    }
    else if(M==64)
    {
        for(i=0 ; i<8 ; i++)
        for(j=0 ; j<8 ; j++)
        {
            //A左上->B左上,A右上临时存入B右上
            for(k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+k][j*8+0] ;
                a1 = A[i*8+k][j*8+1] ;
                a2 = A[i*8+k][j*8+2] ;
                a3 = A[i*8+k][j*8+3] ;
                
                a4 = A[i*8+k][j*8+4] ;
                a5 = A[i*8+k][j*8+5] ;
                a6 = A[i*8+k][j*8+6] ;
                a7 = A[i*8+k][j*8+7] ;

                B[j*8+0][i*8+k] = a0 ;
                B[j*8+1][i*8+k] = a1 ;
                B[j*8+2][i*8+k] = a2 ;
                B[j*8+3][i*8+k] = a3 ;
                
                B[j*8+0][i*8+4+k] = a4 ;
                B[j*8+1][i*8+4+k] = a5 ;
                B[j*8+2][i*8+4+k] = a6 ;
                B[j*8+3][i*8+4+k] = a7 ;
            }
            //B右上->B左下,A左下->A右上
            for(k=0 ; k<4 ; k++)
            {
                a0 = B[j*8+k][i*8+4] ;
                a1 = B[j*8+k][i*8+5] ;
                a2 = B[j*8+k][i*8+6] ;
                a3 = B[j*8+k][i*8+7] ;
                a4 = A[i*8+4][j*8+k] ;
                a5 = A[i*8+5][j*8+k] ;
                a6 = A[i*8+6][j*8+k] ;
                a7 = A[i*8+7][j*8+k] ;
                B[j*8+k][i*8+4] = a4 ;
                B[j*8+k][i*8+5] = a5 ;
                B[j*8+k][i*8+6] = a6 ;
                B[j*8+k][i*8+7] = a7 ;
                B[j*8+4+k][i*8+0] = a0 ;
                B[j*8+4+k][i*8+1] = a1 ;
                B[j*8+4+k][i*8+2] = a2 ;
                B[j*8+4+k][i*8+3] = a3 ;
            }
            //A右下->B右下
            for(k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+4+k][j*8+4] ;
                a1 = A[i*8+4+k][j*8+5] ;
                a2 = A[i*8+4+k][j*8+6] ;
                a3 = A[i*8+4+k][j*8+7] ;
                B[j*8+4][i*8+4+k] = a0 ;
                B[j*8+5][i*8+4+k] = a1 ;
                B[j*8+6][i*8+4+k] = a2 ;
                B[j*8+7][i*8+4+k] = a3 ;
            }
        }
    }
    else
    {
        for(int i=0 ; i<N ; i+=23)
        for(int j=0 ; j<M ; j+=23)
            for(a0=i ; a0<i+23&&a0<N ; a0++)
            for(a1=j ; a1<j+23&&a1<M ; a1++)
                B[a1][a0] = A[a0][a1] ;
    }
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

// M32 N32 
char trans_M32N32_desc[] = "Transpose for M=32 N=32" ;
void trans_M32N32(int M, int N, int A[N][M], int B[M][N])
{
    for(int i=0 ; i<4 ; i++)
        for(int j=0 ; j<4 ; j++)
            for(int k=0 ; k<8 ; k++)
            {
                int a0 = A[i*8+k][j*8+0] ;
                int a1 = A[i*8+k][j*8+1] ;
                int a2 = A[i*8+k][j*8+2] ;
                int a3 = A[i*8+k][j*8+3] ;
                int a4 = A[i*8+k][j*8+4] ;
                int a5 = A[i*8+k][j*8+5] ;
                int a6 = A[i*8+k][j*8+6] ;
                int a7 = A[i*8+k][j*8+7] ;
                B[j*8+0][i*8+k] = a0 ;
                B[j*8+1][i*8+k] = a1 ;
                B[j*8+2][i*8+k] = a2 ;
                B[j*8+3][i*8+k] = a3 ;
                B[j*8+4][i*8+k] = a4 ;
                B[j*8+5][i*8+k] = a5 ;
                B[j*8+6][i*8+k] = a6 ;
                B[j*8+7][i*8+k] = a7 ;
            }
}

// M64 N64 ; not fully marked (miss 1428)
char trans_M64N64_desc[] = "Transpose for M=64 N=64 ; Not fully marked" ;
void trans_M64N64(int M, int N, int A[N][M], int B[M][N])
{
    int a0 , a1 , a2 , a3 ;
    for(int i=0 ; i<8 ; i++)
        for(int j=0 ; j<8 ; j++)
        {
            //A左上,B左上
            for(int k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+k][j*8+0] ;
                a1 = A[i*8+k][j*8+1] ;
                a2 = A[i*8+k][j*8+2] ;
                a3 = A[i*8+k][j*8+3] ;
                B[j*8+0][i*8+k] = a0 ;
                B[j*8+1][i*8+k] = a1 ;
                B[j*8+2][i*8+k] = a2 ;
                B[j*8+3][i*8+k] = a3 ;
            }
            //A左下,B右上
            for(int k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+4+k][j*8+0] ;
                a1 = A[i*8+4+k][j*8+1] ;
                a2 = A[i*8+4+k][j*8+2] ;
                a3 = A[i*8+4+k][j*8+3] ;
                B[j*8+0][i*8+4+k] = a0 ;
                B[j*8+1][i*8+4+k] = a1 ;
                B[j*8+2][i*8+4+k] = a2 ;
                B[j*8+3][i*8+4+k] = a3 ;
            }
            //A右下,B右下
            for(int k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+4+k][j*8+4] ;
                a1 = A[i*8+4+k][j*8+5] ;
                a2 = A[i*8+4+k][j*8+6] ;
                a3 = A[i*8+4+k][j*8+7] ;
                B[j*8+4][i*8+4+k] = a0 ;
                B[j*8+5][i*8+4+k] = a1 ;
                B[j*8+6][i*8+4+k] = a2 ;
                B[j*8+7][i*8+4+k] = a3 ;
            }
            //A右上,B左下
            for(int k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+k][j*8+4] ;
                a1 = A[i*8+k][j*8+5] ;
                a2 = A[i*8+k][j*8+6] ;
                a3 = A[i*8+k][j*8+7] ;
                B[j*8+4][i*8+k] = a0 ;
                B[j*8+5][i*8+k] = a1 ;
                B[j*8+6][i*8+k] = a2 ;
                B[j*8+7][i*8+k] = a3 ;
            }
        }

}

// M64 N64 
char trans_M64N64_v2_desc[] = "Transpose for M=64 N=64 ; Not fully marked" ;
void trans_M64N64_v2(int M , int N , int A[N][M] , int B[M][N])
{
    int a0 , a1 , a2 , a3 , a4 , a5 , a6 , a7 ;
    for(int i=0 ; i<8 ; i++)
        for(int j=0 ; j<8 ; j++)
        {
            //A左上->B左上,A右上临时存入B右上
            for(int k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+k][j*8+0] ;
                a1 = A[i*8+k][j*8+1] ;
                a2 = A[i*8+k][j*8+2] ;
                a3 = A[i*8+k][j*8+3] ;
                
                a4 = A[i*8+k][j*8+4] ;
                a5 = A[i*8+k][j*8+5] ;
                a6 = A[i*8+k][j*8+6] ;
                a7 = A[i*8+k][j*8+7] ;

                B[j*8+0][i*8+k] = a0 ;
                B[j*8+1][i*8+k] = a1 ;
                B[j*8+2][i*8+k] = a2 ;
                B[j*8+3][i*8+k] = a3 ;
                
                B[j*8+0][i*8+4+k] = a4 ;
                B[j*8+1][i*8+4+k] = a5 ;
                B[j*8+2][i*8+4+k] = a6 ;
                B[j*8+3][i*8+4+k] = a7 ;
            }
            //B右上->B左下,A左下->A右上
            for(int k=0 ; k<4 ; k++)
            {
                a0 = B[j*8+k][i*8+4] ;
                a1 = B[j*8+k][i*8+5] ;
                a2 = B[j*8+k][i*8+6] ;
                a3 = B[j*8+k][i*8+7] ;
                a4 = A[i*8+4][j*8+k] ;
                a5 = A[i*8+5][j*8+k] ;
                a6 = A[i*8+6][j*8+k] ;
                a7 = A[i*8+7][j*8+k] ;
                B[j*8+k][i*8+4] = a4 ;
                B[j*8+k][i*8+5] = a5 ;
                B[j*8+k][i*8+6] = a6 ;
                B[j*8+k][i*8+7] = a7 ;
                B[j*8+4+k][i*8+0] = a0 ;
                B[j*8+4+k][i*8+1] = a1 ;
                B[j*8+4+k][i*8+2] = a2 ;
                B[j*8+4+k][i*8+3] = a3 ;
            }
            //A右下->B右下
            for(int k=0 ; k<4 ; k++)
            {
                a0 = A[i*8+4+k][j*8+4] ;
                a1 = A[i*8+4+k][j*8+5] ;
                a2 = A[i*8+4+k][j*8+6] ;
                a3 = A[i*8+4+k][j*8+7] ;
                B[j*8+4][i*8+4+k] = a0 ;
                B[j*8+5][i*8+4+k] = a1 ;
                B[j*8+6][i*8+4+k] = a2 ;
                B[j*8+7][i*8+4+k] = a3 ;
            }
        }

}

// M61 N67
char trans_M61N67_desc[] = "Transpose for M=61 N=67" ;
void trans_M61N67(int M , int N , int A[N][M] , int B[M][N])
{
    int blocksize = 23 ;
    for(int i=0 ; i<N ; i+=blocksize)
        for(int j=0 ; j<M ; j+=blocksize)
            for(int i2=i ; i2<i+blocksize&&i2<N ; i2++)
                for(int j2=j ; j2<j+blocksize&&j2<M ; j2++)
                    B[j2][i2] = A[i2][j2] ;
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 
    registerTransFunction(trans_M32N32, trans_M32N32_desc) ;
    registerTransFunction(trans_M64N64, trans_M64N64_desc) ;
    registerTransFunction(trans_M64N64_v2,trans_M64N64_v2_desc) ;
    registerTransFunction(trans_M61N67,trans_M61N67_desc) ;

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

