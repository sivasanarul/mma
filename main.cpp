/*
 * main.cpp
 *
 *  Created on: Nov 19, 2019
 *      Author: siv0021
 */

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

int main() {

	int m =2;
	int n = 3;
	double epsimin = 0.0000001;

    vector<double> aa;aa.resize(12);
	vector<double> eeen = { 1.0,1.0,1.0 };
	vector<double> eeem = { 1.0,1.0 };
	vector<double> zeron = { 0.0,0.0,0.0 };
	vector<double> zerom = { 0.0,0.0,0.0 };

	vector<double> xval = { 4,3,2 };

	vector<double> xold1 = xval;
	vector<double> xold2 = xval;

	vector<double> xmin = zeron;
	vector<double> xmax = {5.0,5.0,5.0};
	vector<double> low = zeron;
	vector<double> upp = xmax;

	vector<double> c_MMA = { 1000.0,1000.0 };
	vector<double> d     = { 1000.0,1000.0 };

	vector<double> f0val;
	vector<double>df0dx;
	vector<double>fval;
	vector<double>dfdx;

	double a0 = 1;
	vector<double> a = zerom;
	double outeriter = 0.0;
	double maxoutit  = 11;
	double kktol     = 0.0;
	double move = 1.0;
	f0val.resize(1);
	df0dx.resize(3);
	fval.resize(2);


	toy2(xval,f0val,df0dx,fval,dfdx );

	 mma_main1( n,m,outeriter,xval,xmin,xmax,xold1,
				xold2,f0val,df0dx,fval,dfdx,low,
				upp,a0,a,c_MMA,d,move);


cout<< "MMA done";

	return 0;
}
