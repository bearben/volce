#include <global.h>
#include <ineq.h>
#include <z3++.h>
#include <polyvest.h>

#define NDEBUG

//internal string length, inequalities' names and constants' names.
#define STRLEN 30

#ifndef SOLVER_HEADER
#define SOLVER_HEADER

namespace volce {

// enum for parser
enum ENUM_LOGIC {
	UNKNOWN_LOGIC,
	QF_LIA,
	QF_LRA
};

enum SCAN_MODE {
	SM_COMMON,
	SM_SYMBOL,
	SM_COMP_SYM,
	SM_COMMENT,
	SM_STRING
};


enum CMD_TYPE {
	CT_UNKNOWN, CT_EOF,
	//COMMANDS
	CT_ASSERT, CT_CHECK_SAT, CT_CHECK_SAT_ASSUMING,
	CT_DECLARE_CONST, CT_DECLARE_FUN, CT_DECLARE_SORT, 
	CT_DEFINE_FUN, CT_DEFINE_FUN_REC, CT_DEFINE_FUNS_REC, CT_DEFINE_SORT,
	CT_ECHO, CT_EXIT,
	CT_GET_ASSERTIONS, CT_GET_ASSIGNMENT, CT_GET_INFO, 
	CT_GET_MODEL, CT_GET_OPTION, CT_GET_PROOF, 
	CT_GET_UNSAT_ASSUMPTIONS, CT_GET_UNSAT_CORE, CT_GET_VALUE,
	CT_POP, CT_PUSH, CT_RESET, CT_RESET_ASSERTIONS, 
	CT_SET_INFO, CT_SET_LOGIC, CT_SET_OPTION, 
};

enum ERROR_TYPE {
	ERR_UNEXP_EOF, 
	ERR_SYM_MIS, 
	ERR_UNKWN_SYM, 
	ERR_PARAM_MIS,
	ERR_PARAM_NBOOL, 
	ERR_PARAM_NNUM, 
	ERR_PARAM_NSAME, 
	ERR_LOGIC, 
	ERR_MUL_DECL, 
	ERR_MUL_DEF, 
	ERR_NLINEAR, 
	ERR_ZERO_DIVISOR
};


enum NODE_TYPE {
	NT_UNKNOWN, NT_ERROR,
	//CONST
	NT_CONST_BOOL, NT_CONST_NUM, 
	//VAR
	NT_VBOOL, NT_VNUM, NT_INEQ,
	//CORE
	NT_AND, NT_OR, NT_EQ, NT_ITE_BOOL, NT_ITE_NUM,
	//LA, mul and div are temp nodes 
	NT_ADD, NT_MUL, NT_DIV
};



/*
	main class --- solver
*/

class solver {
public:


	//	DAG core 
	//	<t, id, v ,m> --- <type, index, value, multiplier>
	//	const:	id is disabled
	//	var:	id is the index of vbool_list or vnum_list, v is disabled
	//	ineq:	id is the index of ineq_list, v is disabled
	//	oper:	id is the index of op_list, v is enabled only for addition
	class dagc{
	public:
		NODE_TYPE		t;
		unsigned int	id;
		double			v;
		double			m;	//for bool, 1 = not negated, -1 = negated
	
		dagc(NODE_TYPE type = NT_UNKNOWN, unsigned int index = 0, double value = 0, double multiplier = 1) :
				t(type), id(index), v(value), m(multiplier) {};
		~dagc() {};
	
		//overload operator ==
		//only compare type and index, work for boolean-variables and operators
		const bool operator==(const dagc elem) const {	return (t == elem.t && id == elem.id); };

		//check const
		const bool iscbool() 	const { return (t == NT_CONST_BOOL) ? true : false; };
		const bool iscnum() 	const { return (t == NT_CONST_NUM) ? true : false; };
		const bool isconst() 	const { return (iscbool() || iscnum()) ? true : false; };
	
		//check var & ineq
		const bool isvbool() 	const { return (t == NT_VBOOL) ? true : false; };
		const bool isvnum() 	const { return (t == NT_VNUM) ? true : false; };
		const bool isineq()		const { return (t == NT_INEQ) ? true : false; };
		const bool isvar() 		const { return (isvbool() || isvnum() || isineq()) ? true : false; };
	
		//check '='
		const bool iscomp() 	const { return (t == NT_EQ) ? true : false; };
	
		//check ite
		const bool isitebool() 	const { return (t == NT_ITE_BOOL) ? true : false; };
		const bool isitenum() 	const { return (t == NT_ITE_NUM) ? true : false; };
		const bool isite() 		const { return (isitebool() || isitenum()) ? true : false; };
		
		//check other opers
		const bool iserr() 		const { return (t == NT_ERROR) ? true : false; };
		const bool isand() 		const { return (t == NT_AND) ? true : false; };
		const bool isor() 		const { return (t == NT_OR) ? true : false; };
		const bool isadd() 		const { return (t == NT_ADD) ? true : false; };
		const bool ismul() 		const { return (t == NT_MUL) ? true : false; };
		const bool isdiv() 		const { return (t == NT_DIV) ? true : false; };
		
		//check type of operators
		const bool isop()		const { return !(isconst() || isvar()); };
		const bool isnumop()	const { return (isadd() || ismul() || isdiv()) ? true : false; };
		const bool isboolop()	const { return (isand() || isor() || iscomp() || isitebool()) ? true : false; };
		
		//check returned type
		const bool isbool() 	const { return (iscbool() || isvbool() || isineq() || isboolop()) ? true : false; };
		
		const bool isneg() 		const { return !(m > 0); };
	
		const bool bval() 		const { return isneg() ? !(bool)v : (bool)v; };
		const double nval() 	const { return v * m; };

		void negate(){ m = -m; };
		void print() const { std::cout << t << ' ' << id << ' ' << v << ' ' << m << std::endl; };
	
	private:
	};
	
	// value of node
	class dagv {
	public:
		dagv() : v(-1) {};
		dagv(bool val) : v(val) {};
		~dagv() {};
		
		const bool is_unknown() const { return (v == -1) ? true : false; };
		const bool is_true() const { return (!is_unknown() && v) ? true : false; };
		const bool is_false() const { return (!is_unknown() && !v) ? true : false; };
		const dagv operator!() const { return (is_unknown()) ? dagv() : dagv(!v); };
		void negate() { if (!is_unknown()) v = !v; };
	private:
		int v;
	};
	
	// vector of variables
	// values, names and a map of names
	class var_vec {
	public:	
	
		var_vec() {};
		~var_vec() {};
		
		//get & set value
		dagv &val(const unsigned int index) { return vals[index]; };
		dagv &operator()(const unsigned int index) { return vals[index]; };
		
		//get name 
		const std::string name(const unsigned int index) const { return names[index]; };
		const std::string operator[](const unsigned int index) const { return names[index]; };
		
		const unsigned int find(const std::string name) { 
			boost::unordered_map<std::string, unsigned int>::iterator it = map.find(name);
			return (it != map.end()) ? it->second : size(); };
		
		//push back element
		bool push_back(const std::string name) { 
			std::pair<boost::unordered_map<std::string, unsigned int>::iterator, bool>
				p = map.insert(std::pair<std::string, unsigned int>(name, size()));
			if (p.second) {
				//new element
				names.push_back(name);
				return true;
			} else return false; };
		
		//size
		const unsigned int size() const { return names.size(); };
		
		void init_vals() { vals.assign(size(), dagv()); };
		
	protected:
		std::vector<dagv> vals;
		std::vector<std::string> names;
		boost::unordered_map<std::string, unsigned int>	map;
	};
	
	// vector of inequalities
	// values, names, a map of names, ineqs and a map of ineqs
	class ineq_vec : public var_vec {
	public:
		ineq_vec() : var_vec() {};
	
		//retrieve inequalities
		ineqc &ie(const unsigned int index) { return ineqs[index]; };
		ineqc &operator[](const unsigned int index) { return ineqs[index]; };
			
		//push back element
		unsigned int push_back(const ineqc ie) { 
			std::pair<boost::unordered_map<std::vector<double>, unsigned int>::iterator, bool>
				p = ineq_map.insert(std::pair<std::vector<double>, unsigned int>(ie.get_key(), size()));
			if (p.second) {
				//new inequalities
				std::string name = mk_ineq_name(size());
				map.insert(std::pair<std::string, unsigned int>(name, size()));
				names.push_back(name);
				ineqs.push_back(ie); 
				return size() - 1;
			} else return p.first->second; };

	private:
		std::vector<ineqc> ineqs;
		boost::unordered_map<std::vector<double>, unsigned int> ineq_map;
		
		std::string mk_ineq_name(unsigned int n) { char name[STRLEN]; sprintf(name, "_ie%d", n); return name; };
	};
	
	// vector of operators
	// values, types and vectors of children
	class op_vec {
	public:
		op_vec() {};
		~op_vec() {};
		
		//get & set value
		dagv &operator()(const unsigned int index) { return vals[index]; };
		dagv &val(const unsigned int index) { return vals[index]; };
		
		//get type	
		NODE_TYPE &type(const unsigned int index) { return types[index]; };
		
		//get vector of children
		std::vector<dagc> &children(const unsigned int index) { return ops[index]; };
		std::vector<dagc> &operator[](const unsigned int index) { return ops[index]; };
		
		//push type and vector
		void push_back(const NODE_TYPE t, const std::vector<dagc> &elem) {
			types.push_back(t); ops.push_back(elem); };
		
		//resize
		void resize(unsigned int n) { vals.resize(n); types.resize(n); ops.resize(n); };
		
		//size
		const unsigned int size() const { return ops.size(); };
		
		void init_vals() { vals.assign(size(), dagv()); };
		
	private:
		std::vector<dagv> vals;
		std::vector<NODE_TYPE> types;
		std::vector<std::vector<dagc>> ops;
	};
	
	// element of a bunch
	class bunch_elem {
	public:
		std::vector<dagv> ineq_vals;
		std::vector<dagv> vbool_vals;
	
		bunch_elem(ineq_vec ineq_list, var_vec vbool_list) {
			for (unsigned int i = 0; i < ineq_list.size(); i++)
				ineq_vals.push_back(ineq_list(i));
			for (unsigned int i = 0; i < vbool_list.size(); i++)
				vbool_vals.push_back(vbool_list(i));
		};
		~bunch_elem() {};

	private:
	};
	

//attributes
	
	//list of asserts
	std::vector<dagc> assert_list;
	
	//list of operators
	op_vec bop_list;
	op_vec nop_list;
	
	//list of boolean and numeric variables
	var_vec vbool_list;
	var_vec vnum_list;
	
	//list of inequalities
	ineq_vec ineq_list;
	
	//z3
	std::vector<z3::expr>	op_expr;
	std::vector<z3::expr> 	vbool_expr;
	std::vector<z3::expr> 	vnum_expr;
	std::vector<z3::expr> 	ineq_expr;
	z3::context 			z3context;
	z3::solver 				z3solver;
	
	//bunch pool
	std::vector<bunch_elem> bunch_list;
	
	//volume & latte routine
	//basic parameters
	bool		enable_bunch;
	bool 		enable_fact;
	int 		wordlength;
	
	std::string resultdir;	// dir of temp results
	std::string tooldir;	// dir of tools, e.g., vinci, latte
	
	//more statistics
	unsigned int stats_vol_calls;
	unsigned int stats_vol_reuses;
	unsigned int stats_total_dims;
	unsigned int stats_fact_bunches;
	
//methods
	solver(const std::string result_dir, const std::string tool_dir) : 
				z3solver(z3context), 
				enable_bunch(true),
				enable_fact(true),
				wordlength(0),
				resultdir(result_dir),
				tooldir(tool_dir), 
				stats_vol_calls(0),
				stats_vol_reuses(0),
				stats_total_dims(0),
				stats_fact_bunches(0),
				logic(UNKNOWN_LOGIC),
				solving_initialized(false) {};
	solver(const std::string result_dir, const std::string tool_dir, std::string filename) : 
				solver(result_dir, tool_dir) 
	{
		parse_smtlib2_file(filename);
	};
	~solver() {
		op_expr.clear();
		vbool_expr.clear();
		vnum_expr.clear();
		ineq_expr.clear();
	};

	//mk
	const dagc		mk_true() const { 
						if (solving_initialized) err_solving_initialized();
						return dagc(NT_CONST_BOOL, 0, true); };
	const dagc		mk_false() const { 
						if (solving_initialized) err_solving_initialized();
						return dagc(NT_CONST_BOOL, 0, false); };
	const dagc		mk_const(const std::string &s) const { return mk_const(stod(s)); };
	const dagc		mk_const(const double v) const { return dagc(NT_CONST_NUM, 0, v); };
	const dagc		mk_bool_decl(const std::string &name);
	const dagc		mk_var_decl(const std::string &name);
	const dagc		mk_key_bind(const std::string &key, const dagc expr);
	
	const dagc		mk_and(const std::vector<dagc> &params);
	const dagc		mk_or(const std::vector<dagc> &params);
	const dagc		mk_not(const dagc param) const {
						if (solving_initialized) err_solving_initialized();
						//no returned type check
						return dagc(param.t, param.id, param.v, -param.m); };
	const dagc		mk_imply(const std::vector<dagc> &params);
	const dagc		mk_xor(const std::vector<dagc> &params);
	const dagc		mk_eq(const dagc l, const dagc r);
	const dagc		mk_eq(const std::vector<dagc> &params);
	const dagc		mk_distinct(const dagc l, const dagc r);
	const dagc		mk_distinct(const std::vector<dagc> &params);
	const dagc		mk_ite(const dagc c, const dagc l, const dagc r);
	
	const dagc		mk_add(const std::vector<dagc> &params);
	const dagc		mk_neg(const dagc param) const {
						if (solving_initialized) err_solving_initialized();
						//no returned type check
						return dagc(param.t, param.id, param.v, -param.m); };
	const dagc		mk_minus(const std::vector<dagc> &params);
	const dagc		mk_mul(const std::vector<dagc> &params);
	const dagc		mk_div(const dagc l, const dagc r);
	const dagc		mk_le(const dagc l, const dagc r) {
						if (solving_initialized) err_solving_initialized();
						return mk_ineq(false, l, r); };
	const dagc		mk_lt(const dagc l, const dagc r) { 
						if (solving_initialized) err_solving_initialized();
						return mk_not(mk_ineq(false, r, l)); };
	const dagc		mk_ge(const dagc l, const dagc r) {
						if (solving_initialized) err_solving_initialized();
						return mk_ineq(false, r, l); };
	const dagc		mk_gt(const dagc l, const dagc r) {
						if (solving_initialized) err_solving_initialized();
						return mk_not(mk_ineq(false, l, r)); };

	//LOGIC
	const bool		islia() const { return logic == QF_LIA; };
	const bool		islra() const { return logic == QF_LRA; };

	//parse smtlib2 file
	void 			parse_smtlib2_file(const std::string filename);
	
	//SMT solving
	void			z3_init();
	void			z3_init_bounds(const int wordlength);
	const bool		solve();
	
	//volume & latte routine
	const double 	call_latte() { 
						if (!islia()) err_logic_latte();
						return lattice_counting(bsols.back()) * multiplier.back(); };
	const double 	call_latte(const unsigned int index) { 
						if (!islia()) err_logic_latte();
						return lattice_counting(bsols[index]) * multiplier[index]; };
	const double 	call_vinci() { 
						if (!islra()) err_logic_vinci();
						return volume_computation(bsols.back()) * multiplier.back(); };
	const double 	call_vinci(const unsigned int index) {
						if (!islra()) err_logic_vinci();
						return volume_computation(bsols[index]) * multiplier[index]; };
	const double 	call_polyvest(const double coef) {
						if (!islra()) err_logic_polyvest();
						return volume_estimation(bsols.back(), coef) * multiplier.back(); };
	const double 	call_polyvest(const unsigned int index, const double coef) {
						if (!islra()) err_logic_polyvest();
						return volume_estimation(bsols[index], coef) * multiplier[index]; };
	
	//print
	void 			print_ast(const dagc root);
	void 			print_model();
	void			print_ineq(unsigned int index);

private:

//attributes
	ENUM_LOGIC logic;
	
	//vars for parser
	char 			*buffer;
	unsigned long	buflen;
	char			*bufptr;
	unsigned int 	line_number;
	SCAN_MODE 		scan_mode;
	
	boost::unordered_map<std::string, dagc>	key_map;
	
	//var for ineqs
	std::vector<dagc> cond_stack; 
	
	//block any changes of constraints after z3 or vol initialized
	bool 			solving_initialized;
	
	//volume & latte routine
	//data structure
	unsigned int 	nVars;		// = vnum_list.size()
	unsigned int	nFormulas;	// = ineq_list.size()
	arma::mat		bigA;
	arma::vec		bigb;
	int 			*bigop;
	std::vector<long> 	multiplier;	//multiplier of bool sol
	std::vector<int*> 	bsols;		//bunches (assignments)
	
	//map for vol reusing
	std::map<std::vector<int>, double> vol_map;
	
//methods

	//mk
	const dagc		mk_err(const ERROR_TYPE t) const {return dagc(NT_ERROR, t); };
	const dagc		mk_oper(const NODE_TYPE t, std::vector<dagc> &p, double val = 0, double mult = 1) {
						bop_list.push_back(t, p); return dagc(t, bop_list.size() - 1, val, mult); };
	const dagc		mk_num_oper(const NODE_TYPE t, std::vector<dagc> &p, double val = 0, double mult = 1) {
						nop_list.push_back(t, p); return dagc(t, nop_list.size() - 1, val, mult); };

	const dagc		mk_eq_bool(const dagc l, const dagc r);
	// now we transform "l = r" into "(l <= r) and (l >= r)"
	//const dagc	mk_eq_num(const dagc l, const dagc r) { return mk_ineq(true, l, r); };
	
	const dagc		mk_mul(const dagc p, const double m);
	
	const dagc		mk_ineq(const bool iseq, const dagc l, const dagc r);
	const dagc		mk_ineq_core(const bool iseq, const dagc lhs, const dagc rhs);
	const dagc		remk_ineq(const dagc root);
	
	//parse smtlib2 file
	std::string		get_symbol();
	void 			scan_to_next_symbol();
	void			parse_lpar();
	void 			parse_rpar();
	void			skip_to_rpar();
	
	CMD_TYPE		parse_command();
	dagc			parse_expr();
	dagc			parse_and(const unsigned int ln);
	dagc			parse_or(const unsigned int ln);
	dagc			parse_not(const unsigned int ln);
	dagc			parse_imply(const unsigned int ln);
	dagc			parse_xor(const unsigned int ln);
	dagc			parse_eq(const unsigned int ln);
	dagc			parse_distinct(const unsigned int ln);
	dagc			parse_ite(const unsigned int ln);
	dagc			parse_add(const unsigned int ln);
	dagc			parse_neg(const unsigned int ln);
	dagc 			parse_mul(const unsigned int ln);
	dagc			parse_div(const unsigned int ln);
	dagc			parse_le(const unsigned int ln);
	dagc			parse_lt(const unsigned int ln);
	dagc			parse_ge(const unsigned int ln);
	dagc			parse_gt(const unsigned int ln);
	dagc			parse_let(const unsigned int ln);

	//errors & warnings
	void 			err_all(const ERROR_TYPE e, const std::string s = "", const unsigned int ln = 0) const;
	void 			err_all(const dagc e, const std::string s = "", const unsigned int ln = 0) const;
	
	void 			err_unexp_eof() const;
	void 			err_sym_mis(const std::string mis, const unsigned int ln) const;
	void 			err_sym_mis(const std::string mis, const std::string nm, const unsigned int ln) const;
	void 			err_unkwn_sym(const std::string nm, const unsigned int ln) const;
	void 			err_param_mis(const std::string nm, const unsigned int ln) const;
	void 			err_param_nbool(const std::string nm, const unsigned int ln) const;
	void 			err_param_nnum(const std::string nm, const unsigned int ln) const;
	void 			err_param_nsame(const std::string nm, const unsigned int ln) const;
	void 			err_logic(const std::string nm, const unsigned int ln) const;
	void 			err_mul_decl(const std::string nm, const unsigned int ln) const;
	void 			err_mul_def(const std::string nm, const unsigned int ln) const;
	void 			err_nlinear(const std::string nm, const unsigned int ln) const;
	void 			err_zero_divisor(const unsigned int ln) const;
	
	void			err_open_file(const std::string) const;
	void			err_solving_initialized() const;
	void 			err_unbounded_polytope() const;
	void 			err_logic_latte() const;
	void			err_logic_vinci() const;
	void 			err_logic_polyvest() const;
	
	void 			warn_cmd_nsup(const std::string nm, const unsigned int ln) const;

	//SMT solving
	z3::expr		z3_mk_ineq(const unsigned int index);
	z3::expr		z3_mk_term(const term t);
	z3::expr		z3_mk_nconst(const double val);
	z3::expr		z3_mk_op(const unsigned int index);
	z3::expr		z3_mk_expr(const dagc node);
	
	const bool 		get_result();
	const dagv		eval_op(const unsigned int index);
	const dagv		get_val(const dagc node);
	
	// method for bunch strategy
	// input: ineq_list and vbool_list
	// output: ineqs and vbools which should be flipped
	void			get_flip_list(std::vector<bool> &ineq_flip, std::vector<bool> &vbool_flip);
	
	// volume & lattice routine
	void			vol_init();
	const unsigned int get_decided_vars(int *bools, std::vector<int> &vars);
	const bool 		merge_sols(int *source, int *target);
	const unsigned int factorize_bsol(int *bools, std::vector<int*> &pbools);
	const bool		bound_checking(int *bools, unsigned int nRows, std::vector<int> vars);
	const double 	volume_estimation_basic(int *bools, unsigned int nRows, std::vector<int> vars, double coef);
	const double 	volume_estimation(int *boolsol, double coef);
	const double 	volume_computation_light(int *bools, int var);
	const double 	volume_computation_basic(int *bools, unsigned int nRows, std::vector<int> vars);
	const double 	volume_computation(int *boolsol);
	const double 	lattice_counting_light(int *bools, int var);
	const double 	lattice_counting_basic(int *bools, unsigned int nRows, std::vector<int> vars);
	const bool 		lattice_feasible(int *bools, unsigned int nRows, std::vector<int> vars);
	const double 	lattice_counting(int *boolsol);
	
};


}

#endif
