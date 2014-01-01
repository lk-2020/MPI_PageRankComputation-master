//
//  main.c
//  PageRank_Sparse_Input

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "math.h"
#include <time.h>

int N;
int NumOfLines = 0;
double OneByNTimesE;
long dim;


void matrix_mult(double ** matA,float matB[N],float matC[N]);
float vector_norm(float vector[N], int n);
int numprocs, rank, namelen;


int main(int argc, const char * argv[])
{

    // insert code here...
    //printf("Hello, World!\n");
    
    const char * filename = argv[1];
    N = atoi(argv[2]);
    printf("N = %d \n",N);
    
    double psi = 0.00001;
    double delta = 1;
    int i,k;
    int i_sparse = 0;
    
    int done = 1;
  
    double temp, temp1;
    double OneByN;
    OneByN = 1.0/(double)N;
    double vnorm;
    int iterations = 0;
    float* v;
    float* vold;
    float* vdiff;
    char ch;
    double** P1;

    double float_temp;
    
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    double time_initial,time_current,time;

    int greatest = 0;
    
    FILE * pFile;
    FILE * pFile1;
    pFile = fopen (filename,"r");
    pFile1 = pFile;
    fpos_t pos;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);
    MPI_Request request[numprocs];

    
    fgetpos (pFile,&pos);
    
    
    if(pFile==NULL)
    {
        printf("File Read Error \n");
        return 0;
    }
    
    printf("File Opened \n");
    
    while(!feof(pFile)){
        ch = fgetc(pFile);
        if( ch== '\n')
        NumOfLines++;
    }
    
    NumOfLines++;
    
    printf("Number of Lines = %d",NumOfLines);
    
    if(NumOfLines%numprocs)
    {
        NumOfLines += numprocs - (NumOfLines%numprocs);
    }
    
    
    fsetpos(pFile, &pos);
    
    dim = (long)(NumOfLines)/(long)numprocs;
       
    v = malloc(N*sizeof(float));
    vold = malloc(N*sizeof(float));
    vdiff = malloc(N*sizeof(float));
    printf("Array V has been initialised \n");

    
    
    /* Initial approximation of the vector */
    for(i = 0; i<N; i++)
    {
        v[i] = 1;
        vold[i] = 0;
        vdiff[i] = 0;
        
    }

    
    printf("Number of lines = %d \n",NumOfLines);
    P1 = (double**)malloc( NumOfLines*sizeof( double* ));

    
    for(i=0;i<NumOfLines;i++)
    {
        P1[i] = ( double* )malloc( 3* sizeof( double));
    }
    
    OneByNTimesE = OneByN * 0.15;
    
   
    i_sparse =0;
    while(!feof(pFile)){
        fscanf (pFile, "%lf", &temp);
        fscanf (pFile, "%lf", &temp1);
        fscanf (pFile, "%lf", &float_temp);
        
        //    printf("%f \t %f \t %f \n",temp,temp1,float_temp);
        
        P1[i_sparse][0]= temp;
        P1[i_sparse][1]= temp1;
        P1[i_sparse][2]= 0.85*float_temp;//+OneByNTimesE;
        
        i_sparse++;
    }
    
    printf("i_sparse = %d \n", i_sparse);
    i_sparse = 0;

    
    printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);
    
        if(rank == 0)
        {
           float delta1 = 1;
                       
            
            /*
             for(i=0;i<NumOfLines;i++)
             {
             
             printf("%f \t %f \t %f \n",P1[i][0],P1[i][1],P1[i][2]);
             }
             */
            
            clock_t start = clock() ;
            time_initial  = MPI_Wtime();

            
            for(i=1;i<numprocs;i++)
            {
            MPI_Send(&(P1[i*dim]), (int)dim, MPI_FLOAT, i,0, MPI_COMM_WORLD);
                MPI_Isend(&(P1[i*(int)dim]), (int)dim, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request);
            }
             //   MPI_Scatter(&P1, dim, MPI_FLOAT,
               //         &P2, dim, MPI_FLOAT, 0,
                 //       MPI_COMM_WORLD);
            printf("P1 sent \n");
            

           // delta = 1.0;
            printf("psi = %f, delta1 = %f", psi,delta1);
            
            while (psi < delta1) {
                greatest = 0;
                printf("iterarion = %d \n",iterations);
                
                for(i=0;i<N;i++)
                {
                    vold[i] = v[i];
                }
                
                // Matrix multiplcation
                matrix_mult(P1, vold, v);
                
                
                printf("Matrix Mutl Done \n");
                
                //Receive all elements of V from other processes
                for(i=1;i<numprocs;i++)
                {
                //MPI_Recv(&(v[i*dim]), (int)dim, MPI_FLOAT,
                  //       i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Irecv(&(v[i*dim]), (int)dim, MPI_FLOAT,
                              i, 2, MPI_COMM_WORLD, &request[i]);
                    MPI_Wait(&request[i],MPI_STATUS_IGNORE);
                    printf("received from node %d\n",i);


                }
                
                
                
                
               // MPI_Barrier(MPI_COMM_WORLD);

                
                printf("crossed barrier \n");
                
                vnorm = vector_norm(v, N);
                printf("vnorm = %f\n",vnorm);
                
                for(i=0;i<N;i++)
                {
                    //        printf("v = %f \t", v[i]);
                    v[i] = v[i] / vnorm;
                    //        printf("vnormalised = %f \t, vold = %f \n", v[i],vold[i]);
                    if(v[i]>greatest)
                    {
                        greatest = i;
                    }
                    
                    vdiff[i] = vold[i] - v[i];
                }
                
                
                //find the 2 norm and their difference
                
                delta1 = vector_norm(vdiff, N);
                //printf("psi = %f, \t delta1 = %f \n",psi,delta1);
                iterations++;
                
              //  MPI_Bcast(&v, N, MPI_FLOAT,
                //          0, MPI_COMM_WORLD);
                
                
                ////////
                /*
                for(i=1;i<numprocs;i++)
                {
                    MPI_Isend(&(v),(int)N, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request[i]);
                     MPI_Wait(&request[i],MPI_STATUS_IGNORE);
                    // MPI_ISend(&(done), 1, MPI_INT, i,1, MPI_COMM_WORLD);
                }
                 */   
               /////////
                
                if((numprocs>1)){
                    
                    if(psi>delta1)
                    done = 1;
                    
                    ///////
                    /*
                    for(i=1;i<numprocs;i++)
                    {
                         MPI_Isend(&(done),(int)dim, MPI_INT, i, 0, MPI_COMM_WORLD, &request[i]);
                   // MPI_ISend(&(done), 1, MPI_INT, i,1, MPI_COMM_WORLD);
                        MPI_Wait(&request[i],MPI_STATUS_IGNORE);
                    }
                    */
                    ///////
                
                }
                
            }
            
            printf("\n PageRank available, iterations = %d \n",iterations);
            printf("Node Id of page with highest pagerank = %d \n",greatest);
            
            /*
             for(i=0;i<N;i++)
             {
             printf("%f \t",v[i]);
             }
             */
            
            time_current  = MPI_Wtime();
            
            clock_t end = clock() ;
            double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
            
            printf("Time taken = %f \n",elapsed_time);
            
            time  = time_current - time_initial;
            
            
            
            printf("%.3f tid=%i : hello MPI user: machine=%s [NCPU=%i]\n",
                   time, rank, processor_name, numprocs);
            
            
            
            }
    
    
        else{
            
 //           printf("Receiving P2 \n");
 //           MPI_Recv(&(P1), (int)dim, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("P2 received");
            
            
            while(!done)
            {
                
                
                
   //         MPI_Barrier(MPI_COMM_WORLD);
    
            for(i=0;i<N;i++)
            {
                vold[i] = v[i];
            }
            
        
            printf("MatMult on %d\n", rank);
            matrix_mult(&P1[rank*dim], vold, v);
    
            MPI_Isend(&v[rank*dim], (int)dim, MPI_FLOAT, 0, 2, MPI_COMM_WORLD,&request[0]);
            MPI_Wait(&request[0],MPI_STATUS_IGNORE);
                

                // MPI_Barrier(MPI_COMM_WORLD);
                
              //  MPI_Bcast(&v, N, MPI_FLOAT,
                //          0, MPI_COMM_WORLD);
                
               // MPI_Recv(&(v), N, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

     ////////           MPI_Irecv(&(v), N, MPI_FLOAT,  0, MPI_ANY_TAG, MPI_COMM_WORLD, &request[0]);
     ////////           MPI_Wait(&request[0],MPI_STATUS_IGNORE);

                
                vnorm = vector_norm(v, N);
                
                for(i=0;i<N;i++)
                {
                    v[i] = v[i] / vnorm;
                }
                
                //MPI_Recv(&(done), 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     ///////           MPI_Irecv(&(done), 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request[0]);
                
     ///////           MPI_Wait(&request[0],MPI_STATUS_IGNORE);

           
            }
    
        }
  
    printf("Finalize %d", rank);
    MPI_Finalize();

    

    
}

/*Function to compute matrix multiplication*/

void matrix_mult(double ** matA,float matB[N],float matC[N])
{
    int i,k;
    double sum = 0;
    double vsum = 0;
    int index;
    
    for(i=0;i<N;i++)
    {
        vsum+= matB[i];
    }
    
    vsum = vsum * OneByNTimesE;
    
    for(i=0;i<N;i++)
    {
        matC[i] = vsum;
    }
    
    for(i=0;i<dim-1;i++)
    {
        matC[(int)matA[i][0]]+= (matA[i][2])*(matB[(int)matA[i][1]]);
    }
    
}


/*Function to Compute 2-norm of a vector */
float vector_norm(float vector[N], int n)
{
    int i;
    float sum = 0 ;
    
    for(i=0; i<n; i++)
    {
        sum+= vector[i] * vector[i];
    }
    
    return sqrtf(sum);
}


