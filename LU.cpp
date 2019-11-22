/*
 * LU.cpp
 *
 *  Created on: Nov 21, 2019
 *      Author: siv0021
 */




/************** LU Decomposition for solving linear equations ***********/
#include "LU.h"
#include "main.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <bits/stdc++.h>
#include <vector>
using namespace std;



void LU_solver(vector<double> &aa,vector<double> &bb,vector<double> &solution)
{
   vector<double> LU; LU.resize(aa.size());
   Crout( bb.size(),aa,LU);
   solveCrout(bb.size(),LU,bb,solution);
   cout << "LU solver"<< endl;
}

void Crout(double d,vector<double> &aa,vector<double> &LU){
	   for(int k=0;k<d;++k){
	      for(int i=k;i<d;++i){
	         double sum=0.;
	         for(int p=0;p<k;++p)sum+=LU[i*d+p]*LU[p*d+k];
	         LU[i*d+k]=aa[i*d+k]-sum; // not dividing by diagonals
	      }
	      for(int j=k+1;j<d;++j){
	         double sum=0.;
	         for(int p=0;p<k;++p)sum+=LU[k*d+p]*LU[p*d+j];
	         LU[k*d+j]=(aa[k*d+j]-sum)/LU[k*d+k];
	      }
	   }
	}

void solveCrout(int d,vector<double> &LU,vector<double> &b,vector<double> &x){
	   double y[d];
	   for(int i=0;i<d;++i){
	      double sum=0.;
	      for(int k=0;k<i;++k)sum+=LU[i*d+k]*y[k];
	      y[i]=(b[i]-sum)/LU[i*d+i];
	   }
	   for(int i=d-1;i>=0;--i){
	      double sum=0.;
	      for(int k=i+1;k<d;++k)sum+=LU[i*d+k]*x[k];
	      x[i]=(y[i]-sum); // not dividing by diagonals
	   }
	}
