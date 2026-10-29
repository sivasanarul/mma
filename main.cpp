// SPDX-License-Identifier: GPL-3.0-or-later
// Canonical MMA toy problem originally published by Krister Svanberg.

#include "mma.h"
#include "toy2.h"

#include <iostream>
#include <vector>
using namespace std;

int main() {

	int m =2;
	int n = 3;
	vector<double> zeron = { 0.0,0.0,0.0 };
	vector<double> zerom = { 0.0,0.0 };

	vector<double> xval = { 4,3,2 };

	vector<double> xold1 = xval;
	vector<double> xold2 = xval;

	vector<double> xmin = zeron;
	vector<double> xmax = {5.0,5.0,5.0};
	vector<double> low = zeron;
	vector<double> upp = xmax;

	vector<double> c_MMA = { 1000.0,1000.0 };
	vector<double> d     = { 1.0,1.0 };

	vector<double> f0val;
	vector<double>df0dx;
	vector<double>fval;
	vector<double>dfdx;

	double a0 = 1;
	vector<double> a = zerom;
	int outeriter = 0;
	const int maxoutit = 100;
	const double kkttol = 1.0e-7;
	double move = 1.0;
	f0val.resize(1);
	df0dx.resize(3);
	fval.resize(2);
	vector<double> xmma,ymma,out_lam,out_xsi,out_eta,out_mu,out_s,out_low,out_upp={};
	double zmma = 0.0;
    double out_zet = 0.0;
	toy2(xval,f0val,df0dx,fval,dfdx );

double kktnorm = kkttol + 10.0;
int outit = 0;
while (kktnorm>kkttol && outit < maxoutit){
	outit++;
	outeriter++;
	mma_main1( n,m,outeriter,xval,xmin,xmax,xold1,
					xold2,f0val,df0dx,fval,dfdx,low,
					upp,a0,a,c_MMA,d,move,
					xmma,ymma,zmma,out_lam,out_xsi,out_eta,out_mu,out_zet,out_s,out_low,out_upp);
	xold2=xold1;
	xold1=xval;
	xval=xmma;
	toy2(xval,f0val,df0dx,fval,dfdx );
	kktnorm = kktcheck(m,n,xmma,ymma,zmma,out_lam,out_xsi,out_eta,out_mu,out_zet,out_s,xmin,xmax,df0dx,fval,dfdx,a0,a,c_MMA,d);
}
if (kktnorm > kkttol) {
	cerr << "MMA failed to converge after " << outit << " iterations\n";
	return 1;
}
cout << "MMA done after " << outit << " iterations\n";
cout << "KKT norm: " << kktnorm << "\n";
cout << xmma[0] << "\t" << xmma[1] << "\t" << xmma[2] << "\n";
	return 0;
}
