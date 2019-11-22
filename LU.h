/*
 * LU.h
 *
 *  Created on: Nov 21, 2019
 *      Author: siv0021
 */

#ifndef LU_H_
#define LU_H_


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <bits/stdc++.h>
#include <vector>

using namespace std;

void LU_solver(vector<double> &aa,vector<double> &bb,vector<double> &solution);
void Crout(double d,vector<double> &aa,vector<double> &LU);
void solveCrout(int d,vector<double> &LU,vector<double> &b,vector<double> &x);

#endif /* LU_H_ */
