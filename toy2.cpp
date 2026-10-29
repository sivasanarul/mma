// SPDX-License-Identifier: GPL-3.0-or-later
// C++ translation of Krister Svanberg's canonical toy2.m example.




#include "toy2.h"
#include <cmath>
#include <vector>
using namespace std;

void toy2(vector<double> &xval,vector<double> &f0val,vector<double> &df0dx, vector<double> &fval,vector<double> &dfdx){

	f0val[0] = xval[0]*xval[0]  + xval[1]*xval[1] + xval[2]*xval[2] ;
	//auto df0dx = 2 * xval ;

	vector <double> fconst1 = {5,2,1};
	vector <double> fconst2 = {3,4,3};
	double fval1=0,fval2 = 0;

	vector <double> dfdx1,dfdx2;
	for (std::size_t i=0;i<xval.size();i++){
		fval1+=(xval[i] - fconst1[i])*(xval[i] - fconst1[i]);
		fval2+=(xval[i] - fconst2[i])*(xval[i] - fconst2[i]);
        dfdx1.push_back(2*(xval[i] - fconst1[i]));
        dfdx2.push_back(2*(xval[i] - fconst2[i]));
        df0dx[i] = 2*xval[i];
	}
	fval1 = fval1 - 9;
	fval2 = fval2 - 9;

	fval = {fval1,fval2};

	dfdx.resize(6);
	for(std::size_t i=0;i<dfdx1.size();i++){
		dfdx[2*i] = dfdx1[i];
		dfdx[2*i + 1] =  dfdx2[i];
	}



}
