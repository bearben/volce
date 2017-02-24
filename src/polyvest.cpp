/*  polyvest.cpp
 *
 *  Copyright (C) 2015-2017 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of PolyVest.
 *  See COPYING for more information on using this software.
 */


#include "polyvest.h"
#include "glpk.h"
#include <boost/math/distributions/normal.hpp>

#define PI 3.1415926536

using namespace std;
using namespace polyvest;
using namespace arma;


// auxiliary functions
//double abs(double x){
	//return (x > 0) ? x : -x;
//}

double uballVol(int n){
	double vol = 1;
	if (n % 2 == 1){
		int k = (n - 1) / 2;
		vol *= pow(2, n);
		for (int i = 1; i < k + 1; i++) vol *= PI * i;
		for (int i = 1; i < n + 1; i++) vol /= i;
	}else{
		int k = n / 2;
		for (int i = 1; i < k + 1; i++) vol *= PI / i;
	}
	return vol;
}

double zvalue(double x) {

	double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

	double sign = 1;
    
    if (x < 0) {
        sign = -1;
	}
    x = abs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 1.0 - sign*y;

}

/*********** Delete Redundent Hyperplanes ***********/
void polytope::checkHPs(){

	//init GLPK
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, m);
	glp_add_cols(lp, n);

	//disable msg output
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;

	//load constraints
	int *ind = new int[n + 1];
	double *val = new double[n + 1];
	for (int i = 1; i < m + 1; i++){
		for (int j = 1; j < n + 1; j++){
			ind[j] = j, val[j] = A(i - 1, j - 1);
		}
		glp_set_mat_row(lp, i, n, ind, val);
		glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
	}
	delete []ind, delete []val;
	for (int i = 1; i < n + 1; i++)
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

	//feasiblity check
	int num[2];
	for (int i = 1; i < m + 1;){
		glp_set_row_bnds(lp, i, GLP_LO, b(i - 1) + 0.00001, 0);
		glp_set_obj_coef(lp, 1, 1);
		for (int j = 1; j < n; j++)
			glp_set_obj_coef(lp, j + 1, 0);
		glp_simplex(lp, &parm);
		if (glp_get_status(lp) == GLP_NOFEAS){
			num[1] = i;
			glp_del_rows(lp, 1, num);
			A.shed_row(i - 1);
			b.shed_row(i - 1);
			m--;
		}else{
			glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
			i++;
		}
	}
	cout << "Hyperplanes Left: " << m << endl;
	glp_delete_prob(lp);
}

void polytope::genInitE(double &R2, vec &Ori){
	R2 = 0, Ori.zeros();

	//init GLPK
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, m);
	glp_add_cols(lp, n);

	//disable msg output
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;


	//load constraints
	int *ind = new int[n + 1];
	double *val = new double[n + 1];
	for (int i = 1; i < m + 1; i++){
		for (int j = 1; j < n + 1; j++){
			ind[j] = j, val[j] = A(i - 1, j - 1);
		}
		glp_set_mat_row(lp, i, n, ind, val);
		glp_set_row_bnds(lp, i, GLP_UP, 0, b(i - 1));
	}
	delete []ind, delete []val;
	for (int i = 1; i < n + 1; i++)
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

	//get bounds
	for (int i = 0; i < n; i++){
		double max, min;
		for (int j = 0; j < n; j++)
			glp_set_obj_coef(lp, j + 1, 0);

		glp_set_obj_coef(lp, i + 1, 1);
		glp_simplex(lp, &parm);
		max = glp_get_obj_val(lp);
		for (int j = 0; j < n; j++)
			Ori(j) += glp_get_col_prim(lp, j + 1);

		glp_set_obj_coef(lp, i + 1, -1);
		glp_simplex(lp, &parm);
		min = -glp_get_obj_val(lp);
		for (int j = 0; j < n; j++)
			Ori(j) += glp_get_col_prim(lp, j + 1);

		R2 += (max - min) * (max - min);
	}
	Ori = Ori / (2 * n);
	
	glp_delete_prob(lp);
}

bool polytope::Preprocess(){

	if (!check_planes_off) checkHPs();

	double c1 = (2 * pow(n, 2) + pow(1 - n / beta_r, 2)) * (1 - 1.0 / pow(beta_r, 2)) / (2 * pow(n, 2) - 2);
	//double c1 = pow(n, 2) * (1 - 1.0 / pow(beta_r, 2)) / (pow(n, 2) - 1);
	double c2 = (1 - n / beta_r) / (n + 1);
	double c3 = beta_r * beta_r;
	double c4 = 2 * c2 / (1 - 1.0 / beta_r);

	//init E(R2I, 0), T = R2I, ori = 0.
	mat T;
	vec ori(n);
	double R2;
	genInitE(R2, ori);
	T.eye(n, n);
	T = R2 * T;

	vec distance = zeros<vec>(m);
	vec tm = zeros<vec>(m);

	int counter = 0;
	while (++counter > 0){
		int i;
		
		//check if ori in polytope
		distance = b - A * ori;
		for (i = 0; i < m; i++)
			if (distance(i) < 0){
				tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
				break;
			}
		if (i == m){
			//check if small ellipsoid contained in polytope
			for (i = 0; i < m; i++){
				tm(i) = as_scalar(A.row(i) * T * A.row(i).t());
				if (c3 * distance(i) * distance(i) - tm(i) < 0) break;
			}
		}
		
		//terminate if E satisfies two criterions
		if (i == m) break;
		
		vec t = T * A.row(i).t() / sqrt(tm(i));
		ori = ori - t * c2;
		T = c1 * (T - c4 * t * t.t());
	}
	
	if (!msg_off){ 
		cout << "R^2: " << R2 << endl << "Origin: " << endl;
		ori.print();
	}

	//apply affine transformation
	mat Trans;
	
	try{
		Trans = chol(T);
	}catch (const std::runtime_error& ex){
		cout << "\033[2A\r                                               \r";
		return false;		
	}
	
	//cout << Trans << endl;
	b = beta_r * (b - A * ori);
	A = A * Trans.t();

	if (!msg_off) cout << "The number of iterations: " << counter << endl;

	rowvec exp(n);
	exp.ones();
	for (int i = 0; i < n; i++){
		B[i] = b / A.col(i);
		Ai[i] = A / (A.col(i) * exp);
	}
	
	determinant = det(Trans) / pow(beta_r, n);
	
	return true;
}

double polytope::EstimateVol(double epsilon, double delta, double coef = 1.0){
	int k, i, j;

	boost::math::normal dist(0.0, 1.0);
	const double z = boost::math::quantile(dist, 1.0-delta/2);
	const long stepsz = coef * pow((z * l / log(1+epsilon) + z), 2) + 1; //size of sampling
	long counter = 0;

	double *alpha = new double[l];
	long *volK = new long[l];
	memset(alpha, 0, l * sizeof(double));
	memset(volK, 0, l * sizeof(long));
	
	x.zeros();
	for (k = l - 2; k >= 0; k--){
		for (i = volK[k + 1]; i < stepsz; i++){
			counter++;
			double m = 1;
			for (j = 0; j < n; j++) m = walk(k);
			
			if (m < r2[0]) volK[0]++;
			else if (m < r2[k])
				volK[(int)trunc(n * log(m) / (log((double)2) * 2)) + 1]++;
		}
		for (i = 0; i < k; i++){
			volK[k] += volK[i];
		}
		if (volK[k] < stepsz){
			alpha[k] = (double)(stepsz) / volK[k];
			x = x / pow((double)2, (double)1 / n);
		}else alpha[k] = 1;
	}
	
	ofstream outfile("number_of_points.out", ios::app);
	outfile << counter << endl;
	outfile.close();

	vol = uballVol(n) * determinant;
	if (!msg_off) cout << "k\tr^2\t\tvol(k+1)/vol(k)" << endl;

	for (i = 0; alpha[i] > 1 && i < l - 1; i++){
		if (!msg_off) cout << i << "\t" << r2[i] << "\t\t" << alpha[i] << endl;
		vol *= alpha[i];
	}

	delete []alpha;
	delete []volK;

	return vol;
}

double polytope::walk(int k){
	double r, max, min, C = 0;
	int dir = randi(n);
	vec::col_iterator it, itA, end;

	for (it = x.begin(), end = x.end(); it != end; ++it) C += *it * *it;
	C -= x(dir) * x(dir);
	r = sqrt(r2[k + 1] - C);
	max = r - x(dir), min = -r - x(dir);

//A(x + t v) <= b
//Av t <= b - Ax
	vec bound = B[dir] - Ai[dir] * x;
	for (it = bound.begin(), end = bound.end(), itA = A.begin_col(dir); it != end; ++it, ++itA){
		if (*itA > 0){
			if (*it < max) max = *it;
		}else if (*itA < 0)
			if (*it > min) min = *it; 
	}

	double t = x(dir) + randd(max - min) + min;
	x(dir) = t;

	return (C + t * t);
}

