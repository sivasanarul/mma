// SPDX-License-Identifier: GPL-3.0-or-later
// See NOTICE for provenance and copyright information.

#ifndef MMA_H_
#define MMA_H_


#include <vector>

using namespace std;

void mma_main1(int n,int m, int outeriter, vector<double> &xval, vector<double> &xmin, vector<double> &xmax, vector<double> &xold1,
		vector<double> &xold2, vector<double> &f0val, vector<double> &df0dx, vector<double> &fval, vector<double> &dfdx, vector<double> &low,
		vector<double> &upp, double a0, vector<double> &a, vector<double> &c_MMA, vector<double> &d, double move,
		vector<double> &xmma,vector<double> &ymma,double &zmma,vector<double> &lamma,vector<double> &xsimma,
		vector<double> &etamma,vector<double> &mumma,double &zetmma,vector<double> &smma,vector<double> &out_low,vector<double> &out_upp);
void mma_main();



double kktcheck(int m,int n,vector<double> xmma,vector<double> ymma, double zmma, vector<double>out_lam,
		vector<double> out_xsi, vector<double> out_eta,vector<double> out_mu,double out_zet,vector<double> out_s,
		vector<double> xmin,vector<double> xmax,vector<double> df0dx, vector<double> fval,vector<double> dfdx,
		double a0,vector<double> a,vector<double> c_MMA,vector<double> d);


#endif /* MMA_H_ */



