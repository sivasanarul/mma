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
	double sum_volume,delz = 0;
	vector<double> dely,  dellam, diagx, diagxinv, diagy, diagyinv, diaglam, diaglamyi = {};
	dely.resize(m); dellam.resize(n); diagx.resize(n); diagxinv.resize(n); diagy.resize(m); diagyinv.resize(m);
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

	vector<double> x_intermediate = {}; x_intermediate.resize(n);
	vector<double> xsi = {}; xsi.resize(n);
	vector<double> eta = {}; eta.resize(n);
	vector<double> epsvecn = {}; epsvecn.resize(n);
	vector<double> epsvecm = {}; epsvecm.resize(n);
	vector<double> norm_max; norm_max.resize(9);
	vector<double> pjlam, qjlam, gvec; pjlam.resize(n);qjlam.resize(n);gvec.resize(n);
	vector<double> rex, rey , relam, rexsi, reeta, remu, res = {};
	rex.resize(n);rey.resize(m);relam.resize(m);rexsi.resize(n);reeta.resize(n);remu.resize(m);res.resize(m);
	double  rez,rezet  =0.0;
	double norm_rex, norm_rey, norm_rez, norm_relam, norm_rexsi, norm_reeta, norm_remu, norm_rezet, norm_res, norm_sum = 0.0;
	double zet = 1.0, z = 1.0;
	vector<double> dpsidx, delx, gg ;
	dpsidx.resize(n);delx.resize(n);gg.resize(n*m);

	vector<double> blam,bb,Alam,AAr1,AAr2,AA,solution = {};
	vector<double> stepalfa, stepbeta = {};stepalfa.resize(n);stepbeta.resize(n);
	vector<double> dlam, dx,dxsi,deta,dy,ds,dmu; double dz = 0.0; double dzet=0.0;
	dlam.resize(m);dx.resize(n);dxsi.resize(n);deta.resize(n);dy.resize(m);ds.resize(m);dmu.resize(m);
	blam.resize(m);bb.resize(n);Alam.resize(m*m);AAr1.resize(m*m + m);AAr2.resize(m+1);solution.resize(m+1);
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
		std::fill (grad.begin(),grad.end(),0);
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < n; i++) {
				grad[j] += pij[i * m + j] / (upp[i] - x_intermediate[i]) + qij[i * m + j] / (x_intermediate[i] - low[i]);
			}

		}
		for (int i = 0; i < n; i++) {
			dpsidx[i] = pjlam[i] / (std::pow(upp[i] - x_intermediate[i], 2.0)) - qjlam[i] / (std::pow(x_intermediate[i] - low[i], 2.0));
		}

		norm_rex = 0.0; norm_rey = 0.0; norm_rez = 0.0; norm_relam = 0.0;norm_rexsi = 0.0; norm_reeta = 0.0; norm_remu = 0.0; norm_rezet = 0.0; norm_res = 0.0; norm_sum = 0.0;
		for (int i = 0; i < n; i++) {
			rex[i] =std::abs( dpsidx[i] - xsi[i] + eta[i]);
			rexsi[i] = std::abs(xsi[i] * (x_intermediate[i] - alpha[i]) - epsvecn[i]);
			reeta[i] =std::abs( eta[i] * (beta[i] - x_intermediate[i]) - epsvecn[i]);

			norm_rex += rex[i] * rex[i];
			norm_rexsi += rexsi[i] * rexsi[i];
			norm_reeta += reeta[i] * reeta[i];
		}
		for (int i = 0; i < m; i++) {
			rey[i] = std::abs(c_MMA[i] + d[i] - mu[i] - lam[i]);
			relam[i] =std::abs( grad[i] - a[i]*z - y[i] + s[i] - b[i]);
			remu[i] = std::abs(mu[i] * y[i] - epsvecm[i]);
			res[i] = std::abs(lam[i] * s[i] - epsvecm[i]);

			norm_rey += rey[i] * rey[i];
			norm_relam += relam[i] * relam[i];
			norm_remu += remu[i] * remu[i];
			norm_res += res[i] * res[i];

		}
		double a_lam = 0.0;
		for (int i = 0; i<m; i++) {
			a_lam += a[i]*lam[i];
		}

		rez = std::abs(a0 - zet - a_lam);
		rezet = std::abs(zet*z - epsi);
		norm_rez = rez*rez;
		norm_rezet = rezet*rezet;

		norm_sum = norm_rex + norm_rey + norm_rez + norm_relam + norm_rexsi + norm_reeta + norm_remu + norm_rezet + norm_res;
		norm_sum = std::pow(norm_sum, 1.0/2.0);
		norm_max[0]= *max_element(rex.begin(),rex.end());
		norm_max[1]= *max_element(rexsi.begin(),rexsi.end());
		norm_max[2]= *max_element(reeta.begin(),reeta.end());
		norm_max[3]= *max_element(rey.begin(),rey.end());
		norm_max[4]= *max_element(relam.begin(),relam.end());
		norm_max[5]= *max_element(remu.begin(),remu.end());
		norm_max[6]= *max_element(res.begin(),res.end());
		norm_max[7]= rez;
		norm_max[8]= rezet;
		double residmax = *max_element(norm_max.begin(),norm_max.end());

		int ittt = 0;
		while (ittt < 200 && residmax>0.9*epsi) {
			ittt++;
			for (int i = 0; i < n; i++) {
				pjlam[i] = p0[i];
				qjlam[i] = q0[i];
				for (int j = 0; j < m; j++) {
					pjlam[i] += pij[i * m + j] * lam[j];
					qjlam[i] += qij[i * m + j] * lam[j];
				}
			}
			std::fill (grad.begin(),grad.end(),0);
			for (int j = 0; j < m; j++) {
				for (int i = 0; i < n; i++) {
					grad[j] += pij[i * m + j] / (upp[i] - x_intermediate[i]) + qij[i * m + j] / (x_intermediate[i] - low[i]);
				}

			}
			for (int i = 0; i < n; i++) {
				dpsidx[i] = pjlam[i] / (std::pow(upp[i] - x_intermediate[i], 2.0)) - qjlam[i] / (std::pow(x_intermediate[i] - low[i], 2.0));
				delx[i] = dpsidx[i] - epsvecn[i] / (x_intermediate[i] - alpha[i]) + epsvecn[i] / (beta[i] - x_intermediate[i]);
			}

			double sum_lam = 0.0;
			for (int i = 0; i<m; i++) {
				sum_lam = sum_lam + a[i]*lam[i];
			}
			delz = a0 - sum_lam - epsi / z;

			for (int j = 0; j < m; j++) {
				dely[j] = c_MMA[j] + d[j] * y[j] - lam[j] - epsvecm[j] / y[j];
				dellam[j] = grad[j] - a[j] * z - y[j] - b[j] + epsvecm[j] / lam[j];

			}
			for (int i = 0; i < n; i++) {
				diagx[i] = pjlam[i] / (std::pow(upp[i] - x_intermediate[i], 3.0)) + qjlam[i] / (std::pow(x_intermediate[i] - low[i], 3.0));
				diagx[i] = diagx[i] * 2 + xsi[i] / (x_intermediate[i] - alpha[i]) + eta[i] / (beta[i] - x_intermediate[i]);
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
					gg[i*m + j] = pij[i*m + j] * (1 / (std::pow(upp[i] - x_intermediate[i], 2.0))) - qij[i*m + j] * (1 / (std::pow(x_intermediate[i] - low[i], 2.0)));
				}
			}

			for (int j = 0; j < m; j++) {
				blam[j] = dellam[j]+ dely[j]/diagy[j];
				for (int i = 0; i < n; i++) {
					blam[j] = blam[j] - gg[i*m+j]  *delx[i]/diagx[i];
							}
						}


			//vector<vector<int> > vec( n , vector<int> (m, 0));
			vector<double> bb = {};
			bb.insert(bb.end(),blam.begin(),blam.end());
			bb.push_back(delz);

			vector<double> AAr123 = {}; AAr123.resize(n*m);


			for (int i=0;i<n;i++){
				for (int j=0;j<m;j++){
					AAr123[j*n+i] = diagxinv[i]*gg[i*m+j];
				}
			}
			vector<double> AA = {}; AA.resize(m*m);

			double sum = 0;
			for (int i = 0; i < m; i++)
			  {
				for (int j = 0; j < m; j++)
				{
					sum = 0;
					for (int k = 0; k < n; k++)
					{
						sum = sum + AAr123[n*i+k] * gg[m*k+j];
					}
					AA[i+j*m] += sum;
				}
				AA[(m+1)*i] += diaglamyi[i];
			}
			AA.insert(AA.end(),a.begin(),a.end());
			vector<double> AA1 = {}; AA1.insert(AA1.end(),a.begin(),a.end());
			AA1.push_back(-zet/z);
			for (int i =0;i<m+1;i++){
				AA.emplace(AA.begin()+(i+1)*m+i,AA1[i]);
			}
			vector<double> AAr1; vector<double> AAr2;
			vector<double> solution; solution.resize(bb.size());
			LU_solver(AA,bb,solution);
			vector<double> ggdotdlam = {}; ggdotdlam.resize(3);
			std::copy(solution.begin(),solution.begin() +m,dlam.begin());
			dz = solution[m];

			for(int i= 0;i<n;i++){
				for(int j=0;j<m;j++){
				ggdotdlam[i] += gg[i*m+j]*dlam[j];}
			}
			for(int i=0;i<n;i++){
				dx[i] = -delx[i]/diagx[i]- ggdotdlam[i]/diagx[i];
				dxsi[i] = -xsi[i] + 1/(x_intermediate[i]-alpha[i]) - (xsi[i]*dx[i])/(x_intermediate[i]-alpha[i]);////-
				deta[i] = -eta[i] + 1/(beta[i]-x_intermediate[i]) + (eta[i]*dx[i])/(beta[i]-x_intermediate[i]);

			}
			dzet = -zet + epsi/z - dz/z;
			for(int i=0;i<m;i++){
			   dy[i] = -dely[i]/diagy[i]+ dlam[i]/diagy[i];
			   ds[i]   = -s[i] + 1/lam[i] - (s[i]*dlam[i])/lam[i];
			   dmu[i]  = -mu[i] + 1/y[i] - (mu[i]*dy[i])/y[i];
			 }


			for (int i=0;i<n;i++){
				stepalfa[i] = -1.01*dx[i]/(x_intermediate[i]-alpha[i]);
				stepbeta[i] = 1.01*dx[i]/(beta[i]-x_intermediate[i]);
			}

			vector<double> stepd; stepd.resize(8);

			vector<double> stepdy;stepdy.resize(m);
			std::transform (dy.begin(), dy.end(), y.begin(), stepdy.begin(), std::divides<double>());
			std::for_each(stepdy.begin(), stepdy.end(), [](double &el){el *= -1.01; });
			stepd[0]  = *max_element(stepdy.begin(),stepdy.end());
			stepd[1] = -1.01*dz/z;

			vector<double> stepdlam;stepdlam.resize(m);
			std::transform (dlam.begin(), dlam.end(), lam.begin(), stepdlam.begin(), std::divides<double>());
			std::for_each(stepdlam.begin(), stepdlam.end(), [](double &el){el *= -1.01; });
			stepd[2]  = *max_element(stepdlam.begin(),stepdlam.end());

			vector<double> stepdxsi;stepdxsi.resize(n);
			std::transform (dxsi.begin(), dxsi.end(), xsi.begin(), stepdxsi.begin(), std::divides<double>());
			std::for_each(stepdxsi.begin(), stepdxsi.end(), [](double &el){el *= -1.01; });
			stepd[3]  = *max_element(stepdxsi.begin(),stepdxsi.end());

			vector<double> stepdeta;stepdeta.resize(n);
			std::transform (deta.begin(), deta.end(), eta.begin(), stepdeta.begin(), std::divides<double>());
			std::for_each(stepdeta.begin(), stepdeta.end(), [](double &el){el *= -1.01; });
			stepd[4]  = *max_element(stepdeta.begin(),stepdeta.end());

			vector<double> stepdmu;stepdmu.resize(m);
			std::transform (dmu.begin(), dmu.end(), mu.begin(), stepdmu.begin(), std::divides<double>());
			std::for_each(stepdmu.begin(), stepdmu.end(), [](double &el){el *= -1.01; });
			stepd[5]  = *max_element(stepdmu.begin(),stepdmu.end());

			stepd[6]  = 1.01*dzet/zet;
			vector<double> stepds;stepds.resize(m);
			std::transform (ds.begin(), ds.end(), s.begin(), stepds.begin(), std::divides<double>());
			std::for_each(stepds.begin(), stepds.end(), [](double &el){el *= -1.01; });
			stepd[7]  = *max_element(stepds.begin(),stepds.end());

			double stmalfa = *max_element(stepalfa.begin(),stepalfa.end());
			double stmbeta = *max_element(stepbeta.begin(),stepbeta.end());
			double stamlbe = max(stmalfa,stmbeta);
			double stmxx =  *max_element(stepd.begin(),stepd.end());
			double stmalbexx = std::max(stamlbe,stmxx);
			double stminv = std::max(stmalbexx,1.0);
			double steg   = 1.0/stminv;

			vector<double> xold(x_intermediate);
			vector<double> yold(y);
			double zold = z;
			vector<double> lamold(lam);
			vector<double> xsiold(xsi);
			vector<double> muold(mu);
			vector<double> etaold(eta);
			double zetold  = zet;
			vector<double> sold(s);
			//double stmalbexx = max(stamlbe,stmxx);
			//double stminv =  std::max(stmalbexx,1.0);
			//double steg = 1.0/stminv;
			int itto = 0;
			double resinew = 2*norm_sum;
			while (resinew > norm_sum && itto<50){
				itto = itto +1;
				for (int i=0;i<m;i++){
					y[i] = yold[i] + steg*dy[i];
					lam[i] = lamold[i] + steg*dlam[i];
					mu[i]  = muold[i] + steg*dmu[i];
					s[i] = sold[i] + steg*ds[i];
				}
				for (int i=0;i<n;i++){
					xsi[i] = xsiold[i] + steg*dxsi[i];
					x_intermediate[i]   = x_intermediate[i] + steg*dx[i];
					eta[i]   =  etaold[i] + steg*deta[i];
				}
				z = zold + steg*dz;
				zet = zetold + steg*dzet;
				for (int i = 0; i < n; i++) {
					pjlam[i] = p0[i];
					qjlam[i] = q0[i];
					for (int j = 0; j < m; j++) {
						pjlam[i] += pij[i * m + j] * lam[j];
						qjlam[i] += qij[i * m + j] * lam[j];
					}
				}
				std::fill (grad.begin(),grad.end(),0);
				for (int j = 0; j < m; j++) {
					for (int i = 0; i < n; i++) {
						grad[j] += pij[i * m + j] / (upp[i] - x_intermediate[i]) + qij[i * m + j] / (x_intermediate[i] - low[i]);
					}

				}
				for (int i = 0; i < n; i++) {
					dpsidx[i] = pjlam[i] / (std::pow(upp[i] - x_intermediate[i], 2.0)) - qjlam[i] / (std::pow(x_intermediate[i] - low[i], 2.0));
				}

				norm_rex = 0.0; norm_rey = 0.0; norm_rez = 0.0; norm_relam = 0.0; norm_rexsi = 0.0; norm_reeta = 0.0; norm_remu = 0.0; norm_rezet = 0.0; norm_res = 0.0;
				for (int i = 0; i < n; i++) {
					rex[i] =std::abs( dpsidx[i] - xsi[i] + eta[i]);
					rexsi[i] = std::abs(xsi[i] * (x_intermediate[i] - alpha[i]) - epsvecn[i]);
					reeta[i] =std::abs( eta[i] * (beta[i] - x_intermediate[i]) - epsvecn[i]);

					norm_rex += rex[i] * rex[i];
					norm_rexsi += rexsi[i] * rexsi[i];
					norm_reeta += reeta[i] * reeta[i];
				}
				for (int i = 0; i < m; i++) {
					rey[i] = std::abs(c_MMA[i] + d[i] - mu[i] - lam[i]);
					relam[i] =std::abs( grad[i] - a[i]*z - y[i] + s[i] - b[i]);
					remu[i] = std::abs(mu[i] * y[i] - epsvecm[i]);
					res[i] = std::abs(lam[i] * s[i] - epsvecm[i]);

					norm_rey += rey[i] * rey[i];
					norm_relam += relam[i] * relam[i];
					norm_remu += remu[i] * remu[i];
					norm_res += res[i] * res[i];

				}
				double a_lam = 0.0;
				for (int i = 0; i<m; i++) {
					a_lam += a[i]*lam[i];
				}

				rez = std::abs(a0 - zet - a_lam);
				rezet = std::abs(zet*z - epsi);
				norm_rez = rez*rez;
				norm_rezet = rezet*rezet;
				resinew = norm_rex + norm_rey + norm_rez + norm_relam + norm_rexsi + norm_reeta + norm_remu + norm_rezet + norm_res;
				resinew = std::pow(resinew, 1.0/2.0);
				steg = steg/2;
			    }
			norm_sum = resinew;
			steg = 2*steg;
			norm_max[0]= *max_element(rex.begin(),rex.end());
			norm_max[1]= *max_element(rexsi.begin(),rexsi.end());
			norm_max[2]= *max_element(reeta.begin(),reeta.end());
			norm_max[3]= *max_element(rey.begin(),rey.end());
			norm_max[4]= *max_element(relam.begin(),relam.end());
			norm_max[5]= *max_element(remu.begin(),remu.end());
			norm_max[6]= *max_element(res.begin(),res.end());
			norm_max[7]= rez;
			norm_max[8]= rezet;
			residmax = *max_element(norm_max.begin(),norm_max.end());
			cout << "so far so good";


		}
		epsi = 0.1*epsi;
	}
	vector<double> xmma(x_intermediate);
	vector<double> ymma(y);
	vector<double> zmma(z);
	vector<double> lamma(lam);
	vector<double> xsimma(xsi);
	vector<double> etamma(eta);
	vector<double> mumma(mu);
	vector<double> zetmma(zet);
	vector<double> smma(s);




	cout << "Hello MMA!!!" << endl;






}

