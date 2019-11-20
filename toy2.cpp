/*
 * toy2.cpp
 *
 *  Created on: Nov 19, 2019
 *      Author: siv0021
 */




#include "toy2.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <bits/stdc++.h>
#include <vector>
#include <numeric>
using namespace std;

void toy2(vector<double> &xval,vector<double> &f0val,vector<double> &df0dx, vector<double> &fval,vector<double> &dfdx){

	f0val[0] = xval[0]*xval[0]  + xval[1]*xval[1] + xval[2]*xval[2] ;
	//auto df0dx = 2 * xval ;

	vector <double> fconst1 = {5,2,1};
	vector <double> fconst2 = {3,4,3};
	double fval1,fval2 = 0;

	vector <double> dfdx1,dfdx2;
	for (int i=0;i<xval.size();i++){
		fval1+=(xval[i] - fconst1[i])*(xval[i] - fconst1[i]);
		fval2+=(xval[i] - fconst2[i])*(xval[i] - fconst2[i]);
        dfdx1.push_back(2*(xval[i] - fconst1[i]));
        dfdx2.push_back(2*(xval[i] - fconst2[i]));
        df0dx[i] = 2*xval[i];
	}
	fval1 = fval1 - 9;
	fval2 = fval2 - 9;

	fval = {fval1,fval2};

	dfdx.insert( dfdx.end(), dfdx1.begin(), dfdx1.end() );
	dfdx.insert( dfdx.end(), dfdx2.begin(), dfdx2.end() );



	//cout<< std::inner_product( xval.begin(), xval.end(), fconst1.begin(), 0 )<<endl;
	cout << "evaluating function"<< endl;




}
