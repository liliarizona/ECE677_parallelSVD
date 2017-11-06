#include "svd.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stddef.h>
#include <wctype.h>

#define TAG 13
void trim (char *s);
//extern char* strtok_r(char *str, const char *delim, char **nextp);

int main(int argc, char * argv[])
{
	// Initializes MPI
    int numNodes, myRank;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    long long int i,j,k;
    long long int dimensionK;
    long long int rowa, cola;
    char *inputFileDir;
    char *outLeftMatrixDir;
    char *outRightMatrixDir;
    char *outSingularValueDir;
    char *reducedRightDir;
    char *simMatrixDir;
    char *markUSDir;
    char *markSVDir;
        
    int file_size;
    if (argc>0 && myRank == 0)
    {
        inputFileDir=argv[1];
        outLeftMatrixDir=argv[2];
        outRightMatrixDir=argv[3];
        outSingularValueDir=argv[4];
        dimensionK=atoi(argv[5]);
        reducedRightDir=argv[6];
        simMatrixDir=argv[7];
        markUSDir=argv[8];
        markSVDir=argv[9];
        printf("input number: %d\n",argc);
    }
    
    double *a = NULL;

    /*  double a[4][3]={{1.0,1.0,-1.0},{2.0,1.0,0.0},{1.0,-1.0,0.0},{-1.0,2.0,1.0}};
        double b[3][4]={{1.0,1.0,-1.0,-1.0},{2.0,1.0,0.0,2.0},{1.0,-1.0,0.0,1.0}};
        double u[4][4],v[3][3],c[4][3],d[3][4];
        double a[12]={1.0,1.0,-1.0,2.0,1.0,0.0,1.0,-1.0,0.0,-1.0,2.0,1.0};
        double b[12]={1.0,1.0,-1.0,-1.0,2.0,1.0,0.0,2.0,1.0,-1.0,0.0,1.0};
        double u[16],v[9],c[12],d[12];
    */

        /////////////////
        //double * a;
    if (myRank == 0) 
    {
        long long int tempInt;
        char * pch;
        char * pchh;
        rowa=0;
        cola=0;
        FILE *fp;
        fp=fopen(inputFileDir,"rb");
        fseek(fp,0,SEEK_END);
        file_size = ftell(fp);
        //printf( "%d" , file_size );
        char *tmp;
        fseek( fp,0,SEEK_SET);
        tmp=(char *)malloc((file_size+1) * sizeof(char));
        tmp[file_size]='\0';
        fread(tmp,file_size,sizeof(char),fp);
        trim(tmp);

        double *tempArray = (double *) malloc (sizeof(double) * file_size);
        char *outer_ptr=NULL;
        char *inner_ptr=NULL;
        char *tempString;
        char *end;
        double temp;
        int index=0;
        pch = strtok_r(tmp,"\n",&outer_ptr);
        while (pch != NULL)
        {
            //printf("\n one line: \n");
            //printf ("%s\n",pch);
            cola=0;
            pchh = strtok_r(pch,"|",&inner_ptr);
            while (pchh != NULL)
            {
                //printf ("%s\n",pchh);
                tempString=pchh;
                temp=strtod(tempString,&end);
                //printf("temp: %f\n",temp);
                tempArray[index]=temp;
                index=index+1;
                pchh=strtok_r(NULL, "|",&inner_ptr);
                cola=cola+1;
            }
            rowa=rowa+1;
            pch = strtok_r(NULL,"\n",&outer_ptr);
        }
        tempInt=cola;
        cola=rowa;
        rowa=tempInt;

        fclose(fp);
       
        long long int tempCol, tempRow;
        //printf("malloc successed\n");
        a = (double *)malloc(sizeof(double) * (rowa * cola));
        for(i=0;i<rowa*cola;i++)
        {
            //a[i]=tempArray[i];
            tempCol=i%cola;
            tempRow=(i-tempCol)/cola;
            a[tempRow*cola+tempCol]=tempArray[i];
        }
        free(tempArray);
    }
        
    // Passing row a to all processes
    if (myRank == 0) {
        for (j = 1; j < numNodes; j++) {
            printf("Master is sending rowa to process %lld\n", j);
            MPI_Send(&rowa, 1, MPI_LONG_LONG_INT, j, TAG, MPI_COMM_WORLD);
            printf("Master has sent rowa to process %lld\n", j);
        }
    } else {
        MPI_Recv(&rowa, 1, MPI_LONG_LONG_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received rowa = %lld\n", myRank, rowa);
       
    }

    if (myRank == 0) {
        for (j = 1; j < numNodes; j++) {
            printf("Master is sending cola to process %lld\n", j);
            MPI_Send(&cola, 1, MPI_LONG_LONG_INT, j, TAG, MPI_COMM_WORLD);
            printf("Master has sent cola to process %lld\n", j);
        }
    } else {
        MPI_Recv(&cola, 1, MPI_LONG_LONG_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received cola = %lld\n", myRank, cola);
    }

    // MPI broadcast the reading from file to all the processes. 
    if (myRank != 0) {
        a = (double *)malloc(sizeof(double) * rowa * cola);
    }
    if (myRank == 0) {
        for (j = 1; j < numNodes; j++) {
            printf("Master is sending a to process %lld\n", j);
            MPI_Send(a, rowa * cola, MPI_DOUBLE, j, TAG, MPI_COMM_WORLD);
            printf("Master has sent cola to process %lld\n", j);
        }
    } else {
        printf("Process %d received a\n", myRank);
        MPI_Recv(a, rowa * cola, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

//	printf("i: %lld\n",i);
//	printf("a[lastElement]: %f\n",a[i-1]);
    double *u;
    double *v;
    double *c;
    v=(double *)malloc (sizeof(double)*(rowa * cola));
    u=(double *)malloc(sizeof(double)*(rowa * cola));
    c=(double *)malloc (sizeof(double)*(rowa * cola));
    double eps=0.0000001;
    long long int ks,singCount;
    if(cola>rowa)
    {
        ks=cola+1;
        singCount=rowa;
    }
    else
    {
        ks=rowa+1;
        singCount=cola;
    }

    if(dimensionK>singCount)
    {
        dimensionK=singCount;
    }

    if(u==NULL || v==NULL)
    {
        printf("fail malloc\n");
        return 1;
    }
    
    for(i=0;i<cola*cola;i++)
    {
        u[i]=0.0;
    }

    if (myRank == 0) 
    {
        printf("number of singular values: %lld\n", singCount);
    }
    printf("Process %d before calling dluav, cola = %lld, rowa = %lld, ks = %lld\n", myRank, cola, rowa, ks);
    i=dluav(numNodes, myRank, a, cola, rowa, u, v, eps, ks);
    
    //i=dluav(a,rowa,cola,u,v,eps,ks);
    //printf("lenght of leftMatrix: %d \n",(int)(sizeof(v)));
    
    printf("\nProcess %d: dluv Result: i=%lld\n", myRank, i);
    
/* Commented out so the output file writes are not counted. */
/*
	FILE *fpLeft;
	FILE *fpRight;
	FILE *fpSing;
	FILE *fpReRight;
	FILE *fpSimMat;


	fpLeft=fopen(outLeftMatrixDir,"w");
	//printf("\nMAT U Is:\n");
   for(i=0;i<rowa;i++)
	{
	    for(j=0;j<cola;j++)
	    {
            fprintf(fpLeft,"%f",v[i+j*rowa]);
            if(j<(cola-1))
            {
                fprintf(fpLeft,"|");
            }
	    }
	    if(i<(rowa-1))
	    {
            fprintf(fpLeft,"\n");
	    }
	}
	fclose(fpLeft);
	//printf("\n");
	//printf("MAT V IS:\n");
	fpRight=fopen(outRightMatrixDir,"w");
	for(i=0;i<cola*cola;i++)
	{
			//printf("%e ",v[i*3+j]);
			tempCol=i%cola;
			tempRow=(i-tempCol)/cola;
			fprintf(fpRight,"%f",u[i]);
			if(tempCol==cola-1)
			{
			    if(tempRow<(cola-1))
			    {
			        fprintf(fpRight,"\n");
			    }
			}
			else
			{
			    fprintf(fpRight,"|");
			}
	}
	fclose(fpRight);

    double *vk;
    vk=(double *)malloc(sizeof(double)*(cola*dimensionK));
    fpReRight=fopen(reducedRightDir,"w");
	for(i=0;i<cola*cola;i++)
	{
			//printf("%e ",v[i*3+j]);
			tempCol=i%cola;
			tempRow=(i-tempCol)/cola;
			// assign value write the dimension reduced right matrix
			if(tempCol<dimensionK)
			{
                vk[tempRow*dimensionK+tempCol]=u[i];
                fprintf(fpReRight,"%f",vk[tempRow*dimensionK+tempCol]);
                //fprintf(fpReRight,"%f",u[i]);
                if(tempCol==(dimensionK-1))
                {
                    if(tempRow<(cola-1))
                    {
                        fprintf(fpReRight,"\n");
                    }
                }
                else
                {
                    fprintf(fpReRight,"|");
                }
			}
	}
	fclose(fpReRight);

	//printf("\n");
	//printf("MAT A Is:\n");
	fpSing=fopen(outSingularValueDir,"w");
	for(i=0;i<singCount;i++)
	{
		//printf("%e ",a[i*3+j]);
		fprintf(fpSing,"%f",a[i*rowa+i]);
		if(i<singCount-1)
		{
			fprintf(fpSing,"|");
		}
	}
    fclose(fpSing);

    long long int numDocument;
    long long int numTerm;
    numDocument=cola;
    numTerm=rowa;

    double *simMatrix;
    simMatrix=(double *)malloc(sizeof(double)*(numDocument*numDocument));

    double saa,sbb,sab;
    for(i=0;i<numDocument;i++)
    {
        for(j=i;j<numDocument;j++)
        {
            if(j==i)
            {
                simMatrix[i*cola+i]=1;
            }
            else
            {
                saa=0.0;sbb=0.0;sab=0.0;
                for(k=0;k<dimensionK;k++)
                {
                    saa = saa + u[i*cola+k] * u[i*cola+k];
                    sbb = sbb + u[j*cola+k] * u[j*cola+k];
                    sab = sab + u[i*cola+k] * u[j*cola+k];
                }
                simMatrix[i*cola+j]=sab/(sqrt(saa)*sqrt(sbb));
                simMatrix[j*cola+i]=simMatrix[i*cola+j];
            }
        }
    }

    fpSimMat=fopen(simMatrixDir,"w");
	for(i=0;i<numDocument;i++)
	{
        for(j=0;j<numDocument;j++)
        {
            fprintf(fpSimMat,"%f",simMatrix[i*numDocument+j]);
            if(j<(numDocument-1))
            {
                fprintf(fpSimMat,"|");
            }
            else
            {
                if(i<(numDocument-1))
                {
                    fprintf(fpSimMat,"\n");
                }
            }
        }
	}
	fclose(fpSimMat);


    double *singMatrix;
    singMatrix=(double *)malloc(sizeof(double)*(singCount*singCount));
    for(i=0;i<singCount;i++)
    {
        for(j=0;j<singCount;j++)
        {
            if(i==j)
            {
                singMatrix[i*singCount+i]=a[i*rowa+i];
            }
            else
            {
                singMatrix[i*singCount+j]=0.0;
            }
        }
    }
    double *marksVectorUS;
    double *marksVectorSV;
    marksVectorUS=(double *)malloc(sizeof(double)*(numTerm));
    marksVectorSV=(double *)malloc(sizeof(double)*(numDocument));

    for(i=0;i<numTerm;i++)
    {
        temp=0;
        for(j=0;j<singCount;j++)
        {
            temp=temp+v[i+j*rowa]*singMatrix[j*singCount+j];
            printf("%f | %f \n",v[i+j*rowa],singMatrix[j*singCount+j]);
        }
        marksVectorUS[i]=temp;
    }

    for(i=0;i<numDocument;i++)
    {
        temp=0;
        for(j=0;j<singCount;j++)
        {
            for(k=0;k<singCount;k++)
            {
                temp=temp+singMatrix[j*singCount+k]*u[i*cola+k];
            }
        }
        marksVectorSV[i]=temp;
    }

    FILE *fpMarkUS;
    fpMarkUS=fopen(markUSDir,"w");
    for(i=0;i<numTerm;i++)
    {
        fprintf(fpMarkUS,"%f",marksVectorUS[i]);
        fprintf(fpMarkUS,"\n");

    }
    fclose(fpMarkUS);

    FILE *fpMarkSV;
    fpMarkSV=fopen(markSVDir,"w");
    for(i=0;i<numDocument;i++)
    {
        fprintf(fpMarkSV,"%f",marksVectorSV[i]);
        fprintf(fpMarkSV,"\n");
    }
    fclose(fpMarkSV);

	free(a);
	free(u);
	free(v);
	free(c);
	free(vk);
	free(singMatrix);
	free(marksVectorUS);
	free(marksVectorSV);*/
    free(a);
    free(u);
    free(v);
    free(c);
    MPI_Finalize();
    return 0;
}

void trim (char *s)
{
    int i;

    while (isspace (*s)) s++;   // skip left side white spaces
    for (i = strlen (s) - 1; (isspace (s[i])); i--) ;   // skip right side white spaces
    s[i + 1] = '\0';
   // printf ("%s\n", s);
}

