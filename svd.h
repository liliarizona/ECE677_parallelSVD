/*************************************************************
 * 一般实矩阵的奇异值分解，参见《c 常用算法程序集》徐世良P169
 *************************************************************/
#ifndef _SVD_H
#define _SVD_H
#define MAX_ITERA 60
#define MIN_DOUBLE (1e-30)

#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include  <stdio.h>

int dluav(int numNodes, int myRank, double a[], long long int m, long long int n, double u[], double v[], double eps, long long int ka);
//int dluav(double *a,long long int m,long long int n,double *u,double *v,double eps,long long int ka);
void damul(double a[],double b[],int m,int n,int k,double c[]);
void ppp(double a[],double e[],double s[],double v[],int m,int n);
void sss(double fg[2],double cs[2]);

#endif
