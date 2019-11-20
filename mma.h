/*
 * mma.h
 *
 *  Created on: Nov 19, 2019
 *      Author: siv0021
 */

#ifndef MMA_H_
#define MMA_H_


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <bits/stdc++.h>
#include <vector>
using namespace std;

#include "main.h"
void mma_main1(int n,int m, int outeriter, vector<double> &xval, vector<double> &xmin, vector<double> &xmax, vector<double> &xold1,
		vector<double> &xold2, vector<double> &f0val, vector<double> &df0dx, vector<double> &fval, vector<double> &dfdx, vector<double> &low,
		vector<double> &upp, double a0, vector<double> &a, vector<double> &c_MMA, vector<double> &d, double move);
void mma_main();

#endif /* MMA_H_ */



