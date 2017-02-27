/*  vol.cpp
 *
 *  Copyright (C) 2016-2017 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of VolCE.
 *  See COPYING for more information on using this software.
 */

#include "solver.h"
#include "glpk.h"
#include <limits>

bool check_all_zeros(arma::rowvec r) {
	for (arma::rowvec::iterator it = r.begin(); it != r.end(); it++) {
		if (*it != 0) return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////
//// Initialization //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void volce::solver::vol_init(){

	nVars = vnum_list.size();
	nFormulas = ineq_list.size();
	
	//linear constraints	
	bigA.zeros(nFormulas, nVars);
	bigb.zeros(nFormulas);
	bigop = new int[nFormulas];
	
	for (unsigned int i = 0; i < nFormulas; i++) {		
		//retrieve ineq
		volce::ineqc ie = ineq_list[i];
		
		//op: eq 0, le -10, lt -1, gt 1, ge 10
		if (ie.iseq()) bigop[i] = 0;
		else bigop[i] = -10;
		
		//terms
		for (unsigned int j = 0; j < ie.size(); j++) {
			bigA(i, ie[j].id) = ie[j].m;
		}
		bigb(i) = ie.get_const_r();
	}

}


//////////////////////////////////////////////////////////////////////
//// Factorization ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const unsigned int volce::solver::get_decided_vars(int *bools, std::vector<int> &vars){

	for (unsigned int i = 0; i < nVars; i++) {
		bool undecided = true;
		for (unsigned int j = 0; j < nFormulas; j++) {
			if (bools[j] < 0) continue;
			if (bigA(j, i) != 0) {
				undecided = false;
				break;
			}
		}
		if (!undecided) {
			vars.push_back(i);
		}
	}
	
	assert(vars.size() <= nVars);
	
	return vars.size();
}

//merge target solution into source solution
const bool volce::solver::merge_sols(int *source, int *target){

	bool need_merge = false;
	for (unsigned int i = 0; i < nFormulas; i++){
		if (source[i] < 0) continue;
		else if(source[i] == target[i]){
			need_merge = true;
			break;
		}
	}
	if (need_merge)
		for (unsigned int i = 0; i < nFormulas; i++) 
			source[i] = (target[i] == -1) ? source[i] : target[i];
	
	return need_merge;
}

const unsigned int volce::solver::factorize_bsol(int *bools, std::vector<int*> &pbools){
	std::vector<int> constraints;
	std::vector<int> vars;
	unsigned int nc, nv;
	
	for (unsigned int i = 0; i < nFormulas; i++){
		if (bools[i] < 0) continue;
		constraints.push_back(i);
	}
	
	nc = constraints.size();
	nv = get_decided_vars(bools, vars);

	for (unsigned int i = 0; i < nv; i++){
		int *pbsol = new int[nFormulas];
		for (unsigned int j = 0; j < nFormulas; j++) pbsol[j] = -1;
		for (unsigned int j = 0; j < nc; j++)
			if (bigA(constraints[j], vars[i]) != 0) 
				pbsol[constraints[j]] = bools[constraints[j]];
		pbools.push_back(pbsol);
	}

	for (unsigned int i = 0; i < pbools.size(); i++){
		unsigned int j = i + 1;
		while (j < pbools.size()){
			if (merge_sols(pbools[i], pbools[j])){
				delete[] pbools[j];
				pbools.erase(pbools.begin() + j);
				j = i + 1;
			}else
				j++;
		}
	}
	
	return pbools.size();
}

//bound checking for each call
//employing linear programming
const bool volce::solver::bound_checking(int *bools, unsigned int nRows, std::vector<int> vars) {

	unsigned int nVars = vars.size();
	unsigned int counter = 0;

	//init GLPK
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, nRows);
	glp_add_cols(lp, nVars);

	//disable msg output
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;

	//load constraints
	int *ind = new int[nVars + 1];
	double *val = new double[nVars + 1];
	
	for(unsigned int i = 0; i < nFormulas; i++){
		if (bools[i] < 0) continue;
		else counter++;
		
		int cmp = bigop[i];
		assert(cmp != 0);

		//insert one row		
		if ((cmp > 0 && bools[i] == 1) || (cmp < 0 && bools[i] == 0)) {
			for (unsigned int j = 0; j < nVars; j++) {
				ind[j + 1] = j + 1;
				val[j + 1] = bigA(i, vars[j]);
				//std::cout << bigA(i, vars[j]) << ' ';
			}
			glp_set_mat_row(lp, counter, nVars, ind, val);
			glp_set_row_bnds(lp, counter, GLP_LO, bigb[i], 0);
			//std::cout << "LO ";
			//std::cout << bigb[i] << std::endl;
		} else {
			for (unsigned int j = 0; j < nVars; j++) {
				ind[j + 1] = j + 1;
				val[j + 1] = bigA(i, vars[j]);
				//std::cout << bigA(i, vars[j]) << ' ';
			}
			glp_set_mat_row(lp, counter, nVars, ind, val);
			glp_set_row_bnds(lp, counter, GLP_UP, 0, bigb[i]);
			//std::cout << "UP ";
			//std::cout << bigb[i] << std::endl;
		}
	}
	assert(counter == nRows);
	delete []ind, delete []val;
	for (unsigned int i = 1; i < nVars + 1; i++)
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);

	//set obj and apply simplex method
	for (unsigned int j = 1; j < nVars + 1; j++) 
		glp_set_obj_coef(lp, j, 1);
	glp_simplex(lp, &parm);
	
	if (glp_get_status(lp) == GLP_UNBND){
		//no upper bound
		glp_delete_prob(lp);
		return false;
	}

	for (unsigned int j = 1; j < nVars + 1; j++) 
		glp_set_obj_coef(lp, j, -1);
	glp_simplex(lp, &parm);
	
	if (glp_get_status(lp) == GLP_UNBND){
		//no lower bound
		glp_delete_prob(lp);
		return false;
	}

	glp_delete_prob(lp);
	return true;

}

//////////////////////////////////////////////////////////////////////
//// Volume Estimation ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const double volce::solver::volume_estimation_basic(int *bools, unsigned int nRows, std::vector<int> vars, 
		double epsilon, double delta, double coef){
	//nVars: the number of "decided" numeric variables
	//nRows: the number of "decided" linear formulas
	//nFormulas: the number of linear formulas
	//bools[nFormulas], vars[nVars]
	unsigned int nVars = vars.size();
	
	//bound checking
	if (wordlength == 0 && !bound_checking(bools, nRows, vars)) {
		err_unbounded_polytope();
	}
	
	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims += nVars;

	//estimating
	if (wordlength > 0) nRows += 2 * nVars;
	
	polyvest::polytope p(nRows, nVars);
	
	p.msg_off = true;
	p.check_planes_off = true;
	
	unsigned int counter = 0;
	
	//add boundaries
	//Disable bounds when wordlength == 0
	if (wordlength > 0){
		for(unsigned int i = 0; i < nVars; i++){
			//Xi <= 2^(wordlength-1) - 1;
			p.b(counter) = pow(2, wordlength - 1) - 1;
			p.A.row(counter) = arma::zeros<arma::rowvec>(nVars);
			p.A(counter, i) = -1;
			counter++;
			//-Xi <= 2^(wordlength-1)
			p.b(counter) = pow(2, wordlength - 1);
			p.A.row(counter) = arma::zeros<arma::rowvec>(nVars);
			p.A(counter, i) = 1;
			counter++;
		}
	}
	
	for(unsigned int i = 0; i < nFormulas; i++){
		if (bools[i] < 0) continue;
		
		int cmp = bigop[i];
		assert(cmp != 0);

		//insert one row		
		if ((cmp > 0 && bools[i] == 1) || (cmp < 0 && bools[i] == 0)){
			p.b(counter) = -bigb[i];
			for (unsigned int j = 0; j < nVars; j++)
				p.A(counter, j) = -bigA(i, vars[j]);
		}else{
			p.b(counter) = bigb[i];
			for (unsigned int j = 0; j < nVars; j++)
				p.A(counter, j) = bigA(i, vars[j]);
		}
		
		for (unsigned int j = 0; j < counter; j++){
			if (check_all_zeros(p.A.row(counter) + p.A.row(j)) &&
				p.b(counter) == (-1) * p.b(j)){
				//degenerate
				//std::cout << "DEGENERATE" << std::endl;
				return 0;
			}
		}
		
		counter++;
	}

	if (p.Preprocess()){
		p.EstimateVol(epsilon, delta, coef);
		return p.Volume();
	}else{
		//degenerate
		//cout << "DEGENERATE" << std::endl;
		return 0;
	}
}

const double volce::solver::volume_estimation(int *boolsol, double epsilon, double delta, double coef){
	std::vector<int> vars;
	unsigned int nRows = 0;
	double vol = 1;
	
	for (unsigned int i = 0; i < nFormulas; i++){
		if (boolsol[i] < 0) continue;
		if (bigop[i] == 0){
			//equations imply volume of 0
			if (boolsol[i] == 1)
				return 0;
		}else{
			//cmp == 1 || cmp == 10 || cmp == -1 || cmp == -10
			nRows++;
		}
	}
	if (nRows == 0) return 0;

	unsigned int nVars_decided_total = get_decided_vars(boolsol, vars);

	if (!enable_fact){
		unsigned int nVars_undecided = nVars - nVars_decided_total;
		
		//volume of cube consisted of undecided variables
		double cube_vol = pow(pow(2, wordlength), nVars_undecided);
		if (wordlength == 0 && nVars_undecided > 0) {
			//unbounded
			err_unbounded_polytope();
		}
		
		if (nVars_decided_total == 0)
			return cube_vol;
		else if (nVars_decided_total == 1)
			return volume_computation_light(boolsol, vars.back()) * cube_vol;
		else				
			return volume_estimation_basic(boolsol, nRows, vars, epsilon, delta, coef) * cube_vol;
	}
	
	//factorization
	std::vector<int*> pbools;
	unsigned int nVars_tmp = 0;
	unsigned int npbools = factorize_bsol(boolsol, pbools);

	//count the bunch factorized
	if (npbools > 1) stats_fact_bunches++;
	
	//compute each subproblem
	for (unsigned int i = 0; i < npbools; i++){
		vars.clear();
		nRows = 0;
		for (unsigned int j = 0; j < nFormulas; j++){
			if (pbools[i][j] < 0) continue;
			nRows++;
		}
	
		unsigned int nVars_decided = get_decided_vars(pbools[i], vars);
		nVars_tmp += nVars_decided;
		
		if (nVars_decided == 1)
			vol *= volume_computation_light(pbools[i], vars.back());
		else {
			// increase coef while partitions into some pieces
			vol *= volume_estimation_basic(pbools[i], nRows, vars, epsilon, delta, coef * npbools);
		}
	}
	
	assert(nVars_decided_total <= nVars);
	assert(nVars_decided_total == nVars_tmp);
	
	if (wordlength == 0 && nVars - nVars_decided_total > 0) {
		//unbounded
		err_unbounded_polytope();
	} 
	
	if (nVars_decided_total == 0) {
		return pow(pow(2, wordlength), nVars);
	} else
		return vol * pow(pow(2, wordlength), nVars - nVars_decided_total);
}


//////////////////////////////////////////////////////////////////////
//// Volume Computation //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const double volce::solver::volume_computation_light(int *bools, int var){

	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims++;

	double max, min;

	if (wordlength == 0){
		max = std::numeric_limits<double>::max();
		min = -std::numeric_limits<double>::max();
	}else{
		max = pow(2, wordlength - 1) - 1;
		min = -pow(2, wordlength - 1);
	}
	for (unsigned int i = 0; i < nFormulas; i++){
		if(bools[i] < 0) continue;
		
		double v = bigb[i] / bigA(i, var);		
		int cmp = bigop[i];
		
		if (bigA(i, var) < 0) cmp = -cmp;

		if ((cmp > 0 && bools[i] == 1) || (cmp < 0 && bools[i] == 0)){
			if (v > min) min = v;
		}else{
			if (v < max) max = v;
		}
	}
	
	if (max == std::numeric_limits<double>::max() || 
		min == -std::numeric_limits<double>::max()) {
		err_unbounded_polytope();
	} 
	
	if ((max - min) < 0) return 0;
	else return (max - min);
}

const double volce::solver::volume_computation_basic(int *bools, unsigned int nRows, std::vector<int> vars){
	//nVars: the number of "decided" numeric variables
	//nRows: the number of "decided" linear formulas
	//nFormulas: the number of linear formulas
	//bools[nFormulas], vars[nVars]
	unsigned int nVars = vars.size();
	
	//bound checking
	if (wordlength == 0 && !bound_checking(bools, nRows, vars)) {
		err_unbounded_polytope();
	}
	
	//search previous computation result for reusing
	std::vector<int> bools_vec;
	for (unsigned int i = 0; i < nFormulas; i++) 
		bools_vec.push_back(bools[i]);
	std::map<std::vector<int>, double>::iterator vol_map_iter = vol_map.find(bools_vec);
	if (vol_map_iter != vol_map.end()) {
		//result exist
		stats_vol_reuses++;
		return vol_map_iter->second;
	}
	
	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims += nVars;
	
	// compute
	std::string filename = tooldir + "/vinci_input_tmp"; //.ine

	if (wordlength > 0) nRows += 2 * nVars;
	
	std::ofstream ofile;
	
	ofile.open(filename + ".ine");
	if (!ofile.is_open()) {
		err_open_file(filename);
	}

	ofile << "H-representation" << std::endl;
	ofile << "begin" << std::endl;
	ofile << nRows << " " << nVars + 1 << " real" << std::endl;
			
	//add boundaries
	//Disable bounds when wordlength == 0
	if (wordlength > 0) {
		for (unsigned int i = 0; i < nVars; i++){
			ofile << pow(2, wordlength - 1) - 1 << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << ((i == j) ? "-1 " : "0 ");
			ofile << std::endl;
			ofile << pow(2, wordlength - 1) << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << ((i == j) ? "1 " : "0 ");
			ofile << std::endl;
		}
	}
	
	for(unsigned int i = 0; i < nFormulas; i++) {
		if(bools[i] < 0) continue;

		//insert one row
		int cmp = bigop[i];
		assert(cmp != 0);	//if (cmp == 0) continue;
		
		if ((cmp > 0 && bools[i] == 1) || (cmp < 0 && bools[i] == 0)){
			ofile << -bigb[i] << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}else{
			ofile << bigb[i] << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << -bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}
	}
	
	ofile << "end" << std::endl;
	
	ofile.close();
	
	//execute vinci
	std::string cmd = tooldir + "/vinci " + filename + " >/dev/null";
	int proc = system(cmd.c_str());
	
	//read result
	std::ifstream ifile;
	double vol = 0;
	filename = resultdir + "/vinci.result";

	ifile.open(filename);
	if (!ifile.is_open()) {
		err_open_file(filename);
	}
	
	ifile >> vol;
	
	ifile.close();	
	
	//new entry
	vol_map.insert(std::pair<std::vector<int>, double>(bools_vec, vol));
	
	return vol;
	
}

const double volce::solver::volume_computation(int *boolsol){

	std::vector<int> vars;
	unsigned int nRows = 0;
	double vol = 1;
	
	for (unsigned int i = 0; i < nFormulas; i++){
		if (boolsol[i] < 0) continue;
		if (bigop[i] == 0){
			//equations imply volume of 0
			if (boolsol[i] == 1)
				return 0;
		}else{
			//cmp == 1 || cmp == 10 || cmp == -1 || cmp == -10
			nRows++;
		}
	}
	if (nRows == 0) return 0;
	
	unsigned int nVars_decided_total = get_decided_vars(boolsol, vars);

	if (!enable_fact){
		unsigned int nVars_undecided = nVars - nVars_decided_total;
		
		//volume of cube consisted of undecided variables
		double cube_vol = pow(pow(2, wordlength), nVars_undecided);
		if (wordlength == 0 && nVars_undecided > 0) {
			//unbounded
			err_unbounded_polytope();
		}
		
		if (nVars_decided_total == 0)
			return cube_vol;
		else if (nVars_decided_total == 1)
			return volume_computation_light(boolsol, vars.back()) * cube_vol;
		else
			return volume_computation_basic(boolsol, nRows, vars) * cube_vol;
	}
	
	//factorization
	std::vector<int*> pbools;
	unsigned int nVars_tmp = 0;
	unsigned int npbools = factorize_bsol(boolsol, pbools);
	
	//count the bunch factorized
	if (npbools > 1) stats_fact_bunches++;
	
	//compute each subproblem
	for (unsigned int i = 0; i < npbools; i++){
		vars.clear();
		nRows = 0;
		for (unsigned int j = 0; j < nFormulas; j++){
			if (pbools[i][j] < 0) continue;
			nRows++;
		}
	
		unsigned int nVars_decided = get_decided_vars(pbools[i], vars);
		nVars_tmp += nVars_decided;
		
		if (nVars_decided == 1)
			vol *= volume_computation_light(pbools[i], vars.back());
		else
			vol *= volume_computation_basic(pbools[i], nRows, vars);
	}
	
	assert(nVars_decided_total <= nVars);
	assert(nVars_decided_total == nVars_tmp);
	
	if (wordlength == 0 && nVars - nVars_decided_total > 0) {
		//unbounded
		err_unbounded_polytope();
	}
	
	if (nVars_decided_total == 0)
		return pow(pow(2, wordlength), nVars);
	else
		return vol * pow(pow(2, wordlength), nVars - nVars_decided_total);
}

//////////////////////////////////////////////////////////////////////
//// Lattice Counting ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
const double volce::solver::lattice_counting_light(int *bools, int var){

	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims++;

	long max, min;

	if (wordlength == 0){
		max = std::numeric_limits<long>::max();
		min = std::numeric_limits<long>::min();
	}else{
		max = pow(2, wordlength - 1) - 1;
		min = -pow(2, wordlength - 1);
	}
	for (unsigned int i = 0; i < nFormulas; i++){
		if(bools[i] < 0) continue;
		
		double v = bigb[i] / bigA(i, var);		
		int cmp = bigop[i];
		
		if (bigA(i, var) < 0) cmp = -cmp;
		
		if ((cmp == 1 && bools[i] == 1) || (cmp == -10 && bools[i] == 0)){
			//GT
			if (v >= min){
				if (v == ceil(v))
					min = v + 1;
				else
					min = ceil(v);
			}
		}else if ((cmp == 1 && bools[i] == 0) || (cmp == -10 && bools[i] == 1)){
			//LE
			if (v <= max) max = floor(v);
		}else if ((cmp == 10 && bools[i] == 1) || (cmp == -1 && bools[i] == 0)){
			//GE
			if (v >= min) min = ceil(v);
		}else if ((cmp == 10 && bools[i] == 0) || (cmp == -1 && bools[i] == 1)){
			//LT
			if (v <= max){
				if (v == floor(v))
					max = v - 1;
				else
					max = floor(v);
			}
		}else if (cmp == 0){
			//EQ
			if (v >= min) min = ceil(v);
			if (v <= max) max = floor(v);
		}
	}
	
	if (max == std::numeric_limits<long>::max() ||
		min == std::numeric_limits<long>::min()) {
		err_unbounded_polytope();
	}
	
	if ((max - min) < 0) return 0;
	else return (max - min + 1);
}

const double volce::solver::lattice_counting_basic(int *bools, unsigned int nRows, std::vector<int> vars){
	//nVars: the number of "decided" numeric variables
	//nRows: the number of "decided" linear formulas
	//nFormulas: the number of linear formulas
	//bools[nFormulas], vars[nVars]
	unsigned int nVars = vars.size();
	
	//bound checking
	if (wordlength == 0 && !bound_checking(bools, nRows, vars)) {
		err_unbounded_polytope();
	}
	
	//search previous counting result for reusing
	std::vector<int> bools_vec;
	for (unsigned int i = 0; i < nFormulas; i++) 
		bools_vec.push_back(bools[i]);
	std::map<std::vector<int>, double>::iterator vol_map_iter = vol_map.find(bools_vec);
	if (vol_map_iter != vol_map.end()) {
		//result exist
		stats_vol_reuses++;
		return vol_map_iter->second;
	}

	//update stats of vol calls
	stats_vol_calls++;
	stats_total_dims += nVars;
	
	//counting
	std::string filename = tooldir + "/latte_input_tmp"; //.ine

	if (wordlength > 0) nRows += 2 * nVars;
	
	std::ofstream ofile;
	
	ofile.open(filename);
	if (!ofile.is_open()){
		err_open_file(filename);
	}

	ofile << nRows << " " << nVars + 1 << std::endl;
			
	//add boundaries
	//Disable bounds when wordlength == 0
	if (wordlength > 0){
		for (unsigned int i = 0; i < nVars; i++){
			ofile << pow(2, wordlength - 1) - 1 << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << ((i == j) ? "1 " : "0 ");
			ofile << std::endl;
			ofile << pow(2, wordlength - 1) << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << ((i == j) ? "-1 " : "0 ");
			ofile << std::endl;
		}
	}

	for(unsigned int i = 0; i < nFormulas; i++) {
		if(bools[i] < 0) continue;

		//insert one row
		int cmp = bigop[i];
		
		if ((cmp == 1 && bools[i] == 1) || (cmp == -10 && bools[i] == 0)){
			//GT
			ofile << (-1) * (bigb[i] + 1) << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << (-1) * bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}else if ((cmp == 1 && bools[i] == 0) || (cmp == -10 && bools[i] == 1)){
			//LE
			ofile << bigb[i] << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}else if ((cmp == 10 && bools[i] == 1) || (cmp == -1 && bools[i] == 0)){
			//GE
			ofile << (-1) * bigb[i] << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << (-1) * bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}else if ((cmp == 10 && bools[i] == 0) || (cmp == -1 && bools[i] == 1)){
			//LT
			ofile << bigb[i] - 1 << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}else if (cmp == 0){
			//EQ = LE + GE
			ofile << (-1) * bigb[i] << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << (-1) * bigA(i, vars[j]) << " ";
			ofile << std::endl;
			ofile << bigb[i] << " ";
			for (unsigned int j = 0; j < nVars; j++)
				ofile << bigA(i, vars[j]) << " ";
			ofile << std::endl;
		}
	}
	
	ofile.close();

	//execute latte
	std::string cmd = tooldir + "/count --homog " + filename + " >/dev/null 2>/dev/null";
	int proc = system(cmd.c_str());
	
	//read result
	std::ifstream ifile;
	double count = 0;
	filename = resultdir + "/numOfLatticePoints";

	ifile.open(filename);
	if (!ifile.is_open()) {
		err_open_file(filename);
	}
	
	ifile >> count;
	
	ifile.close();
	
	//new entry
	vol_map.insert(std::pair<std::vector<int>, double>(bools_vec, count));
	
	return count;
}

const double volce::solver::lattice_counting(int *boolsol){
	std::vector<int> vars;
	unsigned int nRows = 0;
	double count = 1;
	
	for (unsigned int i = 0; i < nFormulas; i++){
		if (boolsol[i] < 0) continue;
		nRows++;
	}
	if (nRows == 0) return 0;

	unsigned int nVars_decided_total = get_decided_vars(boolsol, vars);

	if (!enable_fact){
		unsigned int nVars_undecided = nVars - nVars_decided_total;
		//cout << "The number of decided variables: " << vars.size() << std::endl;
		
		//volume of cube consisted of undecided variables
		double cube_count = pow(pow(2, wordlength), nVars_undecided);
		if (wordlength == 0 && nVars_undecided > 0) {
			//unbounded
			err_unbounded_polytope();
		}
		
		if (nVars_decided_total == 0)
			return cube_count;
		else if (nVars_decided_total == 1)
			return lattice_counting_light(boolsol, vars.back()) * cube_count;
		else
			return lattice_counting_basic(boolsol, nRows, vars) * cube_count;
	}
	
	//factorization
	std::vector<int*> pbools;
	unsigned int nVars_tmp = 0;
	unsigned int npbools = factorize_bsol(boolsol, pbools);
	
	//count the bunch factorized
	if (npbools > 1) stats_fact_bunches++;
	
	//compute each subproblem
	for (unsigned int i = 0; i < npbools; i++){
		vars.clear();
		nRows = 0;
		for (unsigned int j = 0; j < nFormulas; j++){
			if (pbools[i][j] < 0) continue;
			nRows++;
		}
	
		unsigned int nVars_decided = get_decided_vars(pbools[i], vars);
		nVars_tmp += nVars_decided;
		
		if (nVars_decided == 1)
			count *= lattice_counting_light(pbools[i], vars.back());
		else
			count *= lattice_counting_basic(pbools[i], nRows, vars);
		
	}
	
	assert(nVars_decided_total <= nVars);
	assert(nVars_decided_total == nVars_tmp);
	
	if (wordlength == 0 && nVars - nVars_decided_total > 0) {
		//unbounded
		err_unbounded_polytope();
	}
	
	if (nVars_decided_total == 0)
		return pow(pow(2, wordlength), nVars);
	else
		return count * pow(pow(2, wordlength), nVars - nVars_decided_total);
}




