//============================================================================
// Name        : MMA.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <bits/stdc++.h>
#include <vector>
#include "mma.h"


using namespace std;

void mma_main1(int n,int m, int outeriter, vector<double> &xval, vector<double> &xmin, vector<double> &xmax, vector<double> &xold1,
		vector<double> &xold2, vector<double> &f0val, vector<double> &df0dx, vector<double> &fval, vector<double> &dfdx, vector<double> &low,
		vector<double> &upp, double a0, vector<double> &a, vector<double> &c_MMA, vector<double> &d, double move){
	cout << "Hello MMA1234!!!" << endl; // prints Hello World!!!


	double xmamieps = 0.00001;
	double epsimin = std::sqrt(n + m) * 1e-9;//.0000001
	double raa0 = 0.000001;
	double albefa = 0.1;
	double asyminit = 0.5;
	double asymdec = 0.65;
	double asyminc = 1.05;

	vector<double> b = {}; b.resize(m);



	vector<double> y = {};y.resize(m);
	vector<double> lam = {}; lam.resize(m);
	vector<double> mu = {}; mu.resize(m);
	vector<double> s= {}; s.resize(m);

	vector<double> alpha = {}; alpha.resize(n);
	vector<double> beta = {}; beta.resize(n);
	vector<double> p0 = {}; p0.resize(n);
	vector<double> q0 = {}; q0.resize(n);
	vector<double> pij = {}; pij.resize(m*n);
	vector<double> qij = {}; qij.resize(m*n);
	vector<double> grad = {}; grad.resize(m);
	vector<double> hess = {}; hess.resize(m*m);

	vector<double> gx = {}; gx.resize(4);
	double sum_volume = 0;
	vector<double> dely, delz, dellam, diagx, diagxinv, diagy, diagyinv, diaglam, diaglamyi = {};
	dely.resize(n); delz.resize(n); dellam.resize(n); diagx.resize(n); diagxinv.resize(n); diagy.resize(m); diagyinv.resize(m);
	diaglam.resize(m);diaglamyi.resize(m);

	for (int i = 0; i<m; i++) {
		gx[i] = fval[i];


	}



	if (outeriter < 3) {

		for (int i = 0; i < n; i++) {
			low[i] = xval[i] - asyminit * (xmax[i] - xmin[i]);
			upp[i] = xval[i] + asyminit * (xmax[i] - xmin[i]);
		}
	}
	else {

		for (int i = 0; i < n; i++) {
			double zzz = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]);
			double gamma;
			if (zzz < 0.0) {
				gamma = asymdec;
			}
			else if (zzz > 0.0) {
				gamma = asyminc;
			}
			else {
				gamma = 1.0;
			}
			low[i] = xval[i] - gamma * (xold1[i] - low[i]);
			upp[i] = xval[i] + gamma * (upp[i] - xold1[i]);

			double xmami = std::max(xmamieps, xmax[i] - xmin[i]);
			// double xmami = xmax[i][i] - xmin[i];
			low[i] = std::max(low[i], xval[i] - 10.0 * xmami);
			low[i] = std::min(low[i], xval[i] - 1.0e-2 * xmami);
			upp[i] = std::max(upp[i], xval[i] + 1.0e-2 * xmami);
			upp[i] = std::min(upp[i], xval[i] + 10.0 * xmami);

			double xmi = xmin[i] - 1.0e-6;
			double xma = xmax[i] + 1.0e-6;
			if (xval[i] < xmi) {
				low[i] = xval[i] - (xma - xval[i]) / 0.9;
				upp[i] = xval[i] + (xma - xval[i]) / 0.9;
			}
			if (xval[i] > xma) {
				low[i] = xval[i] - (xval[i] - xmi) / 0.9;
				upp[i] = xval[i] + (xval[i] - xmi) / 0.9;
			}
		}
	}



	for (int i = 0; i < n; ++i) {
		// Compute bounds alpha and beta
		alpha[i] = std::max(xmin[i], low[i] + albefa * (xval[i] - low[i]));
		alpha[i] = std::max(alpha[i], xval[i] - move * (xmax[i] - xmin[i]));
		alpha[i] = std::min(alpha[i], xmax[i]);
		beta[i] = std::min(xmax[i], upp[i] - albefa * (upp[i] - xval[i]));
		beta[i] = std::min(beta[i], xval[i] + move * (xmax[i] - xmin[i]));
		beta[i] = std::max(beta[i], xmin[i]);

		// Objective function
		{
			double dfdxp = std::max(0.0, df0dx[i]);
			double dfdxm = std::max(0.0, -1.0 * df0dx[i]);
			double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
			double pq = 0.001 * (dfdxp + dfdxm) + raa0 * xmamiinv;
			p0[i] = std::pow(upp[i] - xval[i], 2.0) * (dfdxp + pq);
			q0[i] = std::pow(xval[i] - low[i], 2.0) * (dfdxm + pq);
		}

		// Constraints
		for (int j = 0; j < m; j++) {
			double dgdxp = std::max(0.0, dfdx[i * m + j]);
			double dgdxm = std::max(0.0, -1.0 * dfdx[i * m + j]);
			double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
			double pq = 0.001 * (dgdxp + dgdxm) + raa0 * xmamiinv;
			pij[i * m + j] = std::pow(upp[i] - xval[i], 2.0) * (dgdxp + pq);
			qij[i * m + j] = std::pow(xval[i] - low[i], 2.0) * (dgdxm + pq);
		}
	}

	// The constant for the constraints
	for (int j = 0; j < m; j++) {
		b[j] = -gx[j];
		for (int i = 0; i < n; i++) {
			b[j] += pij[i * m + j] / (upp[i] - xval[i]) + qij[i * m + j] / (xval[i] - low[i]);
		}
	}


	////////////////////////////////////   DIP  ///////////////////////

	double x_intermediate[3] = {};
	double xsi[3] = {};
	double eta[3] = {};
	double epsvecn[3] = {};
	double epsvecm[4] = {};
	double pjlam[3], qjlam[3], gvec[3] = {};
	double rex[3], rey[4], rez, relam[4], rexsi[3], reeta[3], remu[4], rezet, res[4] = {};
	double norm_rex, norm_rey, norm_rez, norm_relam, norm_rexsi, norm_reeta, norm_remu, norm_rezet, norm_res, norm_sum = 0.0;
	double zet = 1.0, z = 1.0;
	double dpsidx[3], delx[3], gg[3 * 4] = {};

	double blam[m],bb[n],Alam[m*m],AAr1[m*m + m],AAr2[m+1],AA,solution[m+1] = {};
	for (int i = 0; i<n; i++) {
		x_intermediate[i] = (alpha[i] + beta[i]) / 2;
		xsi[i] = std::max(1 / (x_intermediate[i] - alpha[i]), 1.0);
		eta[i] = std::max(1 / (beta[i] - x_intermediate[i]), 1.0);

	}
	for (int i = 0; i<m; i++) {
		y[i] = 1.0;
		mu[i] = std::max(1.0, 0.5*c_MMA[i]);
		s[i] = 1.0;
		lam[i] = 1.0;
	}
	double tol = epsimin; // 1.0e-9*sqrt(m+n);
	double epsi = 1.0;
	double err = 1.0;
	int loop;

	while (epsi > tol) {
		for (int i = 0; i<n; i++) {
			epsvecn[i] = 1.0*epsi;
		}
		for (int i = 0; i<m; i++) {
			epsvecm[i] = 1.0*epsi;
		}
		for (int i = 0; i < n; i++) {
			pjlam[i] = p0[i];
			qjlam[i] = q0[i];
			for (int j = 0; j < m; j++) {
				pjlam[i] += pij[i * m + j] * lam[j];
				qjlam[i] += qij[i * m + j] * lam[j];
			}
		}
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < n; i++) {
				grad[j] += pij[i * m + j] / (upp[i] - x_intermediate[i]) + qij[i * m + j] / (x_intermediate[i] - low[i]);
			}

		}
		for (int i = 0; i < n; i++) {
			dpsidx[i] = pjlam[i] / (std::pow(upp[i] - x_intermediate[i], 2.0)) - qjlam[i] / (std::pow(x_intermediate[i] - low[i], 2.0));
		}


		for (int i = 0; i < n; i++) {
			rex[i] = dpsidx[i] - xsi[i] + eta[i];
			rexsi[i] = xsi[i] * (x_intermediate[i] - alpha[i]) - epsvecn[i];
			reeta[i] = eta[i] * (beta[i] - xval[i]) - epsvecn[i];

			norm_rex += rex[i] * rex[i];
			norm_rexsi += rexsi[i] * rexsi[i];
			norm_reeta += reeta[i] * reeta[i];
		}
		for (int i = 0; i < m; i++) {
			rey[i] = c_MMA[i] + d[i] - mu[i] - lam[i];
			relam[i] = grad[i] - 1 - y[i] + s[i] - b[i];
			remu[i] = mu[i] * y[i] - epsvecm[i];
			res[i] = lam[i] * s[i] - epsvecm[i];

			norm_rey += rey[i] * rey[i];
			norm_relam += relam[i] * relam[i];
			norm_remu += remu[i] * remu[i];
			norm_res += res[i] * res[i];

		}
		double a_lam = 0.0;
		for (int i = 0; i<m; i++) {
			a_lam += a[i]*lam[i];
		}

		rez = a0 - zet - a_lam;
		rezet = zet*z - epsi;
		norm_rez = rez*rez;
		norm_rezet = rezet*rezet;
		norm_sum = norm_rex + norm_rey + norm_rez + norm_relam + norm_rexsi + norm_reeta + norm_remu + norm_rezet + norm_res;
		int ittt = 0;
		while (ittt < 200) {
			ittt++;
			for (int i = 0; i < n; i++) {
				pjlam[i] = p0[i];
				qjlam[i] = q0[i];
				for (int j = 0; j < m; j++) {
					pjlam[i] += pij[i * m + j] * lam[j];
					qjlam[i] += qij[i * m + j] * lam[j];
				}
			}
			for (int j = 0; j < m; j++) {
				for (int i = 0; i < n; i++) {
					grad[j] += pij[i * m + j] / (upp[i] - xval[i]) + qij[i * m + j] / (xval[i] - low[i]);
				}

			}
			for (int i = 0; i < n; i++) {
				dpsidx[i] = pjlam[i] / (std::pow(upp[i] - xval[i], 2.0)) - qjlam[i] / (std::pow(xval[i] - low[i], 2.0));
				delx[i] = dpsidx[i] - epsvecn[i] / (x_intermediate[i] - alpha[i]) + epsvecn[i] / (beta[i] - x_intermediate[i]);
			}

			double sum_lam = 0.0;
			for (int i = 0; i<m; i++) {
				sum_lam = sum_lam + lam[i];
			}
			for (int j = 0; j < m; j++) {
				dely[j] = c_MMA[j] + d[j] * y[j] - lam[j] - epsvecm[j] / y[j];
				delz[j] = a0 - sum_lam - epsi / z;
				dellam[j] = grad[j] - a[j] * z - y[j] - b[j] + epsvecm[j] / lam[j];

			}
			for (int i = 0; i < n; i++) {
				diagx[i] = pjlam[i] / (std::pow(upp[i] - xval[i], 3.0)) + qjlam[i] / (std::pow(xval[i] - low[i], 3.0));
				diagx[i] = diagx[i] * 2 + xsi[i] / (xval[i] - alpha[i]) + eta[i] / (beta[i] - xval[i]);
				diagxinv[i] = 1 / diagx[i];
			}

			for (int j = 0; j < m; j++) {
				diagy[j] = d[j] + mu[j] / y[j];
				diagyinv[j] = 1 / diagy[j];
				diaglam[j] = s[j] / lam[j];
				diaglamyi[j] = diaglam[j] + diagyinv[j];

			}

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					gg[i*n + j] = pij[i*n + j] * (1 / (std::pow(upp[i] - xval[i], 2.0))) - qij[i*n + j] * (1 / (std::pow(upp[i] - xval[i], 2.0)));
				}
			}

		    if (m<n){






		    }

		    /*for (int i=0;i<m;i++){
		    	dy[i] = -dely[i]/diagy[i] + dlam[i]/diagy[i];
		    	dmu[i] = -mu[i] + epsvecm[i]/y[i] - (mu[i]*dy[i])/y[i];
		    	ds[i]  = -s[i] + epsvecm[i]/lam[i] - (s[i]*dlam[i])/lam[i];
		        }

		    for (int i=0;i<n;i++){
		   		    	dxsi[i] = -xsi[i] + epsvecn[i]/(xval[i]-alpha[i]) - (xsi[i]*dx[i])/(xval[i]-alpha[i]);
		   		    	deta[i] = -eta[i] + epsvecn[i]/(beta[i]-xval[i])  - (eta[i]*dx[i])/(beta[i]-xval[i]);


		   		}

		    dzet = -zet + epsi/z - zet*dz/z; */




		}
		epsi = 0;
	}



	cout << "Hello MMA!!!" << endl;






}

void mma_main() {
	cout << "Hello MMA!!!" << endl; // prints Hello World!!!


	cout<<"Hello World";
	vector<double> aaaa = { 1, 45, 54, 71, 76, 12 };


	cout<< aaaa[3]+1;
	      cout<<*max_element(aaaa.begin(),aaaa.end());

	      vector<double> aa;aa.resize(12);

	         for (int i = 0; i<12;i++){
	             aa[i] = i*100;
	         }
	int m = 4;
	int n = 3;
	int iter = 0;
	int iteration = 0;
	double  xval[3]= { 1.0,1.0,1.0 };
	double xmin[3] = { 0.0010,0.0010,0.0010 };
	double xmax[3] = { 3.0,3.0,3.0 };
	double xold1[3] = { 1.0,1.0,1.0 };
	double xold2[3] = { 1.0,1.0,1.0 };
	double f0val = 0.0;
	double df0dx[3] = { 0.0,0.0,0.0 };
	double fval[4] = { 0.75,1.0,0.75,0.0 };
	double dgdx[12] = { -.5625,-.25,-.0625,1.0,-.125,-.50,-.125,1.00,-.0625,-.25,-.5625,1.0 };
	double low[3] = { 0.0010,0.0010,0.0010 };
	double upp[3] = { 3.0,3.0,3.0 };
	double a0 = 1.0;
	double a[4] = { 1.0,1.0,1.0,0 };
	double c_MMA = 1000;
	double d[4] = { 0.0,0.0,0.0,0.0 };
	double move = 1.0;

	double xmamieps = 0.00001;
	double epsimin = std::sqrt(n + m) * 1e-9;//.0000001
	double raa0 = 0.000001;
	double albefa = 0.1;
	double asyminit = 0.5;
	double asymdec = 0.65;
	double asyminc = 1.05;

	double b[4] = {};



	double y[4] = {};
	double lam[4] = {};
	double mu[4] = {};
	double s[4] = {};

	double alpha[3] = {};
	double beta[3] = {};
	double p0[3] = {};
	double q0[3] = {};
	double pij[3 * 4] = {};
	double qij[3 * 4] = {};
	double grad[4] = {};
	double hess[4 * 4] = {};

	double gx[4] = {};
	double sum_volume = 0;
	double dely[3], delz[3], dellam[3], diagx[3], diagxinv[3], diagy[4], diagyinv[4], diaglam[4], diaglamyi[4] = {};
	for (int i = 0; i<m; i++) {
		gx[i] = fval[i];


	}



	if (iteration < 3) {

		for (int i = 0; i < n; i++) {
			low[i] = xval[i] - asyminit * (xmax[i] - xmin[i]);
			upp[i] = xval[i] + asyminit * (xmax[i] - xmin[i]);
		}
	}
	else {

		for (int i = 0; i < n; i++) {
			double zzz = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]);
			double gamma;
			if (zzz < 0.0) {
				gamma = asymdec;
			}
			else if (zzz > 0.0) {
				gamma = asyminc;
			}
			else {
				gamma = 1.0;
			}
			low[i] = xval[i] - gamma * (xold1[i] - low[i]);
			upp[i] = xval[i] + gamma * (upp[i] - xold1[i]);

			double xmami = std::max(xmamieps, xmax[i] - xmin[i]);
			// double xmami = xmax[i][i] - xmin[i];
			low[i] = std::max(low[i], xval[i] - 10.0 * xmami);
			low[i] = std::min(low[i], xval[i] - 1.0e-2 * xmami);
			upp[i] = std::max(upp[i], xval[i] + 1.0e-2 * xmami);
			upp[i] = std::min(upp[i], xval[i] + 10.0 * xmami);

			double xmi = xmin[i] - 1.0e-6;
			double xma = xmax[i] + 1.0e-6;
			if (xval[i] < xmi) {
				low[i] = xval[i] - (xma - xval[i]) / 0.9;
				upp[i] = xval[i] + (xma - xval[i]) / 0.9;
			}
			if (xval[i] > xma) {
				low[i] = xval[i] - (xval[i] - xmi) / 0.9;
				upp[i] = xval[i] + (xval[i] - xmi) / 0.9;
			}
		}
	}



	for (int i = 0; i < n; ++i) {
		// Compute bounds alpha and beta
		alpha[i] = std::max(xmin[i], low[i] + albefa * (xval[i] - low[i]));
		alpha[i] = std::max(alpha[i], xval[i] - move * (xmax[i] - xmin[i]));
		alpha[i] = std::min(alpha[i], xmax[i]);
		beta[i] = std::min(xmax[i], upp[i] - albefa * (upp[i] - xval[i]));
		beta[i] = std::min(beta[i], xval[i] + move * (xmax[i] - xmin[i]));
		beta[i] = std::max(beta[i], xmin[i]);

		// Objective function
		{
			double dfdxp = std::max(0.0, df0dx[i]);
			double dfdxm = std::max(0.0, -1.0 * df0dx[i]);
			double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
			double pq = 0.001 * (dfdxp + dfdxm) + raa0 * xmamiinv;
			p0[i] = std::pow(upp[i] - xval[i], 2.0) * (dfdxp + pq);
			q0[i] = std::pow(xval[i] - low[i], 2.0) * (dfdxm + pq);
		}

		// Constraints
		for (int j = 0; j < m; j++) {
			double dgdxp = std::max(0.0, dgdx[i * m + j]);
			double dgdxm = std::max(0.0, -1.0 * dgdx[i * m + j]);
			double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
			double pq = 0.001 * (dgdxp + dgdxm) + raa0 * xmamiinv;
			pij[i * m + j] = std::pow(upp[i] - xval[i], 2.0) * (dgdxp + pq);
			qij[i * m + j] = std::pow(xval[i] - low[i], 2.0) * (dgdxm + pq);
		}
	}

	// The constant for the constraints
	for (int j = 0; j < m; j++) {
		b[j] = -gx[j];
		for (int i = 0; i < n; i++) {
			b[j] += pij[i * m + j] / (upp[i] - xval[i]) + qij[i * m + j] / (xval[i] - low[i]);
		}
	}


	////////////////////////////////////   DIP  ///////////////////////

	double x_intermediate[3] = {};
	double xsi[3] = {};
	double eta[3] = {};
	double epsvecn[3] = {};
	double epsvecm[4] = {};
	double pjlam[3], qjlam[3], gvec[3] = {};
	double rex[3], rey[4], rez, relam[4], rexsi[3], reeta[3], remu[4], rezet, res[4] = {};
	double norm_rex, norm_rey, norm_rez, norm_relam, norm_rexsi, norm_reeta, norm_remu, norm_rezet, norm_res, norm_sum = 0.0;
	double zet = 1.0, z = 1.0;
	double dpsidx[3], delx[3], gg[3 * 4] = {};

	double blam[m],bb[n],Alam[m*m],AAr1[m*m + m],AAr2[m+1],AA,solution[m+1] = {};
	for (int i = 0; i<n; i++) {
		x_intermediate[i] = (alpha[i] + beta[i]) / 2;
		xsi[i] = std::max(1 / (x_intermediate[i] - alpha[i]), 1.0);
		eta[i] = std::max(1 / (beta[i] - x_intermediate[i]), 1.0);

	}
	for (int i = 0; i<m; i++) {
		y[i] = 1.0;
		mu[i] = std::max(1.0, 0.5*c_MMA);
		s[i] = 1.0;
		lam[i] = 1.0;
	}
	double tol = epsimin; // 1.0e-9*sqrt(m+n);
	double epsi = 1.0;
	double err = 1.0;
	int loop;

	while (epsi > tol) {
		for (int i = 0; i<n; i++) {
			epsvecn[i] = 1.0*epsi;
		}
		for (int i = 0; i<m; i++) {
			epsvecm[i] = 1.0*epsi;
		}
		for (int i = 0; i < n; i++) {
			pjlam[i] = p0[i];
			qjlam[i] = q0[i];
			for (int j = 0; j < m; j++) {
				pjlam[i] += pij[i * m + j] * lam[j];
				qjlam[i] += qij[i * m + j] * lam[j];
			}
		}
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < n; i++) {
				grad[j] += pij[i * m + j] / (upp[i] - x_intermediate[i]) + qij[i * m + j] / (x_intermediate[i] - low[i]);
			}

		}
		for (int i = 0; i < n; i++) {
			dpsidx[i] = pjlam[i] / (std::pow(upp[i] - x_intermediate[i], 2.0)) - qjlam[i] / (std::pow(x_intermediate[i] - low[i], 2.0));
		}


		for (int i = 0; i < n; i++) {
			rex[i] = dpsidx[i] - xsi[i] + eta[i];
			rexsi[i] = xsi[i] * (x_intermediate[i] - alpha[i]) - epsvecn[i];
			reeta[i] = eta[i] * (beta[i] - xval[i]) - epsvecn[i];

			norm_rex += rex[i] * rex[i];
			norm_rexsi += rexsi[i] * rexsi[i];
			norm_reeta += reeta[i] * reeta[i];
		}
		for (int i = 0; i < m; i++) {
			rey[i] = c_MMA + d[i] - mu[i] - lam[i];
			relam[i] = grad[i] - 1 - y[i] + s[i] - b[i];
			remu[i] = mu[i] * y[i] - epsvecm[i];
			res[i] = lam[i] * s[i] - epsvecm[i];

			norm_rey += rey[i] * rey[i];
			norm_relam += relam[i] * relam[i];
			norm_remu += remu[i] * remu[i];
			norm_res += res[i] * res[i];

		}
		double a_lam = 0.0;
		for (int i = 0; i<m; i++) {
			a_lam += a[i]*lam[i];
		}

		rez = a0 - zet - a_lam;
		rezet = zet*z - epsi;
		norm_rez = rez*rez;
		norm_rezet = rezet*rezet;
		norm_sum = norm_rex + norm_rey + norm_rez + norm_relam + norm_rexsi + norm_reeta + norm_remu + norm_rezet + norm_res;
		int ittt = 0;
		while (ittt < 200) {
			ittt++;
			for (int i = 0; i < n; i++) {
				pjlam[i] = p0[i];
				qjlam[i] = q0[i];
				for (int j = 0; j < m; j++) {
					pjlam[i] += pij[i * m + j] * lam[j];
					qjlam[i] += qij[i * m + j] * lam[j];
				}
			}
			for (int j = 0; j < m; j++) {
				for (int i = 0; i < n; i++) {
					grad[j] += pij[i * m + j] / (upp[i] - xval[i]) + qij[i * m + j] / (xval[i] - low[i]);
				}

			}
			for (int i = 0; i < n; i++) {
				dpsidx[i] = pjlam[i] / (std::pow(upp[i] - xval[i], 2.0)) - qjlam[i] / (std::pow(xval[i] - low[i], 2.0));
				delx[i] = dpsidx[i] - epsvecn[i] / (x_intermediate[i] - alpha[i]) + epsvecn[i] / (beta[i] - x_intermediate[i]);
			}

			double sum_lam = 0.0;
			for (int i = 0; i<m; i++) {
				sum_lam = sum_lam + lam[i];
			}
			for (int j = 0; j < m; j++) {
				dely[j] = c_MMA + d[j] * y[j] - lam[j] - epsvecm[j] / y[j];
				delz[j] = a0 - sum_lam - epsi / z;
				dellam[j] = grad[j] - a[j] * z - y[j] - b[j] + epsvecm[j] / lam[j];

			}
			for (int i = 0; i < n; i++) {
				diagx[i] = pjlam[i] / (std::pow(upp[i] - xval[i], 3.0)) + qjlam[i] / (std::pow(xval[i] - low[i], 3.0));
				diagx[i] = diagx[i] * 2 + xsi[i] / (xval[i] - alpha[i]) + eta[i] / (beta[i] - xval[i]);
				diagxinv[i] = 1 / diagx[i];
			}

			for (int j = 0; j < m; j++) {
				diagy[j] = d[j] + mu[j] / y[j];
				diagyinv[j] = 1 / diagy[j];
				diaglam[j] = s[j] / lam[j];
				diaglamyi[j] = diaglam[j] + diagyinv[j];

			}

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) {
					gg[i*n + j] = pij[i*n + j] * (1 / (std::pow(upp[i] - xval[i], 2.0))) - qij[i*n + j] * (1 / (std::pow(upp[i] - xval[i], 2.0)));
				}
			}

		    if (m<n){






		    }

		    /*for (int i=0;i<m;i++){
		    	dy[i] = -dely[i]/diagy[i] + dlam[i]/diagy[i];
		    	dmu[i] = -mu[i] + epsvecm[i]/y[i] - (mu[i]*dy[i])/y[i];
		    	ds[i]  = -s[i] + epsvecm[i]/lam[i] - (s[i]*dlam[i])/lam[i];
		        }

		    for (int i=0;i<n;i++){
		   		    	dxsi[i] = -xsi[i] + epsvecn[i]/(xval[i]-alpha[i]) - (xsi[i]*dx[i])/(xval[i]-alpha[i]);
		   		    	deta[i] = -eta[i] + epsvecn[i]/(beta[i]-xval[i])  - (eta[i]*dx[i])/(beta[i]-xval[i]);


		   		}

		    dzet = -zet + epsi/z - zet*dz/z; */




		}
		epsi = 0;
	}



	cout << "Hello MMA!!!" << endl;

}
