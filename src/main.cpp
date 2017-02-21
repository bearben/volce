#include <solver.h>

#define MAX_DIRSTR_SIZE 255

using namespace std;

void printUsage(char* exec_name){
	cout << endl;
	cout << "USAGE: " << exec_name << " [OPTION]... <INPUT-FILE> [OUTPUT-FILE]" << endl << endl;
	cout << "OPTIONS:" << endl;
	cout << "  -P,-p       \t Enable PolyVest for volume approximation. The input of" << endl;
    cout << "              \t linear inequalities are reals. By default, VolCE calls" << endl;
    cout << "              \t PolyVest. And it requires that all the numeric variables" << endl;
    cout << "              \t are reals (QF_LRA logic)." << endl;
    cout << "  -V,-v       \t Enable Vinci for volume computation. The input of linear" << endl;
    cout << "              \t inequalities are reals (QF_LRA logic)." << endl;
    cout << "  -L,-l       \t Enable LattE to count integer solutions. The input of" << endl;
    cout << "              \t linear inequalities should be integers (QF_LIA logic)." << endl;
    cout << "  -w={0,1,...}\t Specify the word length of numeric variables in bit-wise." << endl;
    cout << "              \t Set word length to 0 will disable this feature. By default," << endl;
    cout << "              \t the word length is 0." << endl;
    cout << "  -maxc={real}\t Set maximum sampling coefficient of PolyVest, which is an" << endl;
    cout << "              \t upper bound. Generally, the larger coefficient, the more" << endl;
    cout << "              \t accurate, however, the slower. The default value is 1.0." << endl;
    cout << "  -minc={real}\t Like option -maxc, this option set the minimum sampling" << endl;
	cout << "              \t coefficient of PolyVest. The default value is 0.01." << endl;
	cout << "  -bunch={0,1}\t Enable (1) or disable (0) the bunch strategy. By default," << endl;
	cout << "              \t this strategy is enabled." << endl;
	cout << "  -fact={0,1}\t Enable (1) or disable (0) the factorization strategy." << endl;
	cout << "              \t It can be very efficient for problems whose variables " << endl;
	cout << "              \t are less connective. By default, this strategy is enabled." << endl;
	cout << "  -verb={0,1} \t The verbosity of output. Positive value will enable pretty" << endl;
	cout << "              \t print. Otherwise, only print the final result. The default" << endl;
	cout << "              \t value is 1." << endl;
    cout << "INPUT-FILE:" << endl;
    cout << "  .smt2       \t SMT-LIBv2 language input." << endl;
    //cout << "  .vs or other\t Recognize as VolCE style input." << endl;
	cout << endl;
}

// calculate the coefficient for two-round strategy
double cal_coef(double vol, double mvol, double minc, double maxc){
	double t = 2 * maxc * vol / mvol;
	t = (t <= minc) ? minc : (t > maxc) ? maxc : t;
	return t;
}

int main(int argc, char **argv) {

	bool 	polyvest 	= false;
	bool 	vinci 		= false;
	bool 	latte 		= false;
	int 	wordlength 	= 0;
	double 	maxc 		= 1;
	double 	minc 		= 0.01;
	bool	bunch		= true;
	bool 	fact 		= true;
	int 	verbosity 	= 1;

	//auxiliary variables
	//clock_t c_start, c_end;
	string input_file = "";
	string output_file = "";

    if (argc == 1){
    	cout << "error: lack input file." << endl;
    	cout << "Use '-h' or '--help' for help." << endl;
        exit(0);
    }
    
 	////////////////////////////////////////////////////////////////////// 
    //parsing arguments
    for (int i = 1; i < argc; i++){
    	string argument = argv[i];
    	int offset = argument.find('=');
    	string key = argument.substr(0, offset);
    	string value = argument.substr(offset + 1);

		if (key == "-P" || key == "-p"){
			//PolyVest
			polyvest = true;
		}else if (key == "-V" || key == "-v"){
			//Vinci
			vinci = true;
		}else if (key == "-L" || key == "-l"){
			//LattE
			latte = true;
		}else if (key == "-W" || key == "-w"){
			//wordlength
			try {
				wordlength = stoi(value);
			}catch (const invalid_argument&){
				cout << "error: Invalid value \"" << value << "\" for argument \"" << key << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}else if (key == "-maxc"){
			//max coefficient
			try {
				maxc = stod(value);				
			}catch (const invalid_argument&){
				cout << "error: Invalid value \"" << value << "\" for argument \"" << key << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}else if (key == "-minc"){
			//min coefficient
			try {
				minc = stod(value);
			}catch (const invalid_argument&){
				cout << "error: Invalid value \"" << value << "\" for argument \"" << key << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}else if (key == "-bunch"){
			// disable bunch strategy
			try {
				bunch = stoi(value);
			}catch (const invalid_argument&){
				cout << "error: Invalid value \"" << value << "\" for argument \"" << key << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}else if (key == "-fact"){
			//disable factorization
			try {
				fact = stoi(value);
			}catch (const invalid_argument&){
				cout << "error: Invalid value \"" << value << "\" for argument \"" << key << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}else if (key == "-verb"){
			//set verbosity
			try{
				verbosity = stoi(value);
				if (verbosity != 0) verbosity = 1;
			}catch (const invalid_argument&){
				cout << "error: Invalid value \"" << value << "\" for argument \"" << key << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}else if (key == "-h" || key == "--help"){
			//help
			printUsage(argv[0]);
			return 0;
		}else{
			if (input_file == "")
				input_file = argv[i];
			else if (output_file == "")
				output_file = argv[i];
			else{
				cout << "error: Multiple inputs or outputs designated \"" << input_file << "\", \"" 
						<< output_file << "\" and \"" << argv[i] << "\"." << endl;
				cout << "Use '-h' or '--help' for help." << endl;
				exit(0);
			}
		}
    }
    if (!polyvest && !vinci && !latte) {
    	//enable polyvest in default
    	polyvest = true;
    }
    
    //////////////////////////////////////////////////////////////////////
    //print parameters
    
    if (verbosity > 0) {
	    cout << endl << "====================================" << endl;
  		cout << "============ Parameters ============" << endl;
  		cout << "====================================" << endl << endl;
  	}  	  	
    
    if (polyvest) {
    	cout << "-P\t\tEnable PolyVest." << endl;

		if (maxc > 0) 
			cout << "-maxc=" << maxc << "\t\tSet max coefficient to " << maxc << "." << endl;
		else {
			cout << "error: The value of max coefficient should be positive." << endl;
			cout << "Use '-h' or '--help' for help." << endl;
			exit(0);
		}

		if (minc > 0) 
			cout << "-minc=" << minc << "\tSet min coefficient to " << minc << "." << endl;
		else {
			cout << "error: The value of min coefficient should be positive." << endl;
			cout << "Use '-h' or '--help' for help." << endl;
			exit(0);
		}

    }
    
    if (vinci) {
    	cout << "-V\t\tEnable Vinci." << endl;
    }
    
    if (latte) {
    	cout << "-L\t\tEnable LattE." << endl;
    }
    
	if (wordlength == 0)
		cout << "-w=0\t\tDisabled default bounds." << endl;
	else if (wordlength > 0)
		cout << "-w=" << wordlength <<	"\t\tSet word length to " << wordlength << "." << endl;
	else {
		cout << "error: The value of wordlength should be 0 or positive." << endl;
		cout << "Use '-h' or '--help' for help." << endl;
		exit(0);
	}
	
	if (!bunch) {
		cout << "-bunch=0\tBunch strategy turned off." << endl;
	}else{
		cout << "-bunch=1\tBunch strategy turned on." << endl;
	}
   
	if (!fact) {
		cout << "-fact=0\t\tConstraints factorization turned off." << endl;
	}else{
		cout << "-fact=1\t\tConstraints factorization turned on." << endl;
	}
	
	if (!verbosity) {
		cout << "-verb=0\t\tPretty print turned off." << endl;
	} else {
		cout << "-verb=1\t\tPretty print turned on." << endl;
	}
	
	cout << endl;
	
  	cout << "Input File: \"" << input_file << "\"" << endl;
  	if (output_file != "") cout << "Output File: \"" << output_file << "\"" << endl;
  	cout << endl;
    
 	//////////////////////////////////////////////////////////////////////
 
 	//check path  
   	char current_absolute_path[MAX_DIRSTR_SIZE];
   	char execution_path[MAX_DIRSTR_SIZE];
	//obtain absolute path
	int cnt = readlink("/proc/self/exe", current_absolute_path, MAX_DIRSTR_SIZE);
	if (cnt < 0 || cnt >= MAX_DIRSTR_SIZE)
	{
	    cout << "error: Failed to get absolute path." << endl;
	    exit(0);
	}
	for (int i = cnt; i >= 0; --i)
	    if (current_absolute_path[i] == '/'){
	        current_absolute_path[i] = '\0';
	        break;
	    }
	if (verbosity > 0) cout << "VolCE Directory: " << current_absolute_path << endl;
  
    //obtain working directory    
	string execdir = getcwd(execution_path, sizeof(execution_path));
	if (verbosity > 0) cout << "Working Directory: " << execdir;

	//compare working directory and absolute path
	if (strcmp(current_absolute_path, execdir.c_str()) != 0){
		if (verbosity > 0) cout << endl << endl << 
			"warning: Working directory is not the absolute path of VolCE." << endl;
	}else{
		if (verbosity > 0) cout << "   ... OK" << endl;
	}
	
	string bindir = execdir + "/bin";

 	//////////////////////////////////////////////////////////////////////

	//initialize solver
	volce::solver s(execdir, bindir, input_file);
	s.enable_bunch = bunch;
	s.enable_fact = fact;
	s.wordlength = wordlength;

	if (verbosity > 0) {
  		cout << endl << "====================================" << endl;
  		cout << "========== Problem Scale ===========" << endl;
  		cout << "====================================" << endl << endl;
 		cout << "Number of Boolean variables: " << s.vbool_list.size() << endl;
		cout << "Number of inequalities: " << s.ineq_list.size() << endl;
		cout << "Number of numeric variables: " << s.vnum_list.size() << endl;
		//cout << "Number of Boolean operators: " << s.bop_list.size() << endl;
		//cout << "Number of numeric operators: " << s.nop_list.size() << endl;
		cout << "Number of assertions: " << s.assert_list.size() << endl;
	}

	s.z3_init();
	
	// obtain all bunches
  	if (verbosity > 0){
  		cout << endl << "====================================" << endl;
  		cout << "=========== SMT Solving ============" << endl;
  		cout << "====================================" << endl << endl;
 	}
 	
	unsigned int count = 0;
	printf("#Bunches: %d\n", count);
	while (s.solve())
		printf("\033[1A\r#Bunches: %d\n", ++count);
	
	//cout << "#Bunches: " << s.bunch_list.size() << endl;
	
	if (count == 0) {
		cout << endl << "The problem is unsat." << endl;
  		
  		//print to output
  		ofstream fout(output_file);
		fout << "unsat" << endl;
	  	fout.close();	

		return 1;
	}
	
 	//////////////////////////////////////////////////////////////////////
	
	double total_latte = 0;
	double total_vinci = 0;
	double total_polyvest = 0;
	
	// lattice counting routine
	if (latte) {

  		if (verbosity > 0){
  			cout << endl << "====================================" << endl;
  			cout << "============== LattE ===============" << endl;
  			cout << "====================================" << endl << endl;
   			cout << "Index\tCount" << endl;
  		}
		
		for (unsigned int i = 0; i < s.bunch_list.size(); i++) {
		
			double res = s.call_latte(i);
			
			if (verbosity > 0) {
				cout << i + 1 << "\t" << res << endl;
			}
			
			total_latte += res;
			
		}
	}
	
	// volume computation routine
	if (vinci) {

  		if (verbosity > 0){
  			cout << endl << "====================================" << endl;
  			cout << "============== Vinci ===============" << endl;
  			cout << "====================================" << endl << endl;
   			cout << "Index\tVolume" << endl;
  		}
		
		for (unsigned int i = 0; i < s.bunch_list.size(); i++) {
		
			double res = s.call_vinci(i);
			
			if (verbosity > 0) {
				cout << i + 1 << "\t" << res << endl;
			}
			
			total_vinci += res;
			
		}	

	}
	
	// volume estimation routine
	if (polyvest) {
	
	  	if (verbosity > 0){
  			cout << endl << "====================================" << endl;
  			cout << "============ PolyVest ==============" << endl;
  			cout << "====================================" << endl << endl;
  		}
  		
  		double *vol = new double[s.bunch_list.size()];
  		double maxvol = 0;

  		//first round
   		if (verbosity > 0){
   			cout << "FIRST ROUND" << endl;
   			cout << "Index\tVolume" << endl;
   		}
   		for (unsigned int i = 0; i < s.bunch_list.size(); i++){

			vol[i] = s.call_polyvest(i, minc);
  			
  			if (verbosity > 0) {
  				cout << i + 1 << "\t" << vol[i] << endl;
  			}
  				
  			if (maxvol < vol[i]) maxvol = vol[i];
  		}

  		//second round
  		if (verbosity > 0){
  			cout << endl << "SEC & LAST ROUND" << endl;
  			cout << "Index\tCoef\tVolume" << endl;
  		}
  		
 	  	for (unsigned int i = 0; i < s.bunch_list.size(); i++){
 	  		
 	  		if (vol[i] == 0) continue;
 	  		
 	  		double coef = cal_coef(vol[i], maxvol, minc, maxc);
 	  		
 	  		if (coef > minc){
 	  		
	  			vol[i] = s.call_polyvest(i, coef);
	  			
  				if (verbosity > 0) { 
  					cout << i + 1 << "\t" << coef << "\t" << vol[i] << endl;
  				}
  			}
			total_polyvest += vol[i];
			
  		}

	}

	if (verbosity > 0) {	
  		cout << endl << "====================================" << endl;
  		cout << "=========== Statistics =============" << endl;
  		cout << "====================================" << endl << endl;
  		cout << "The number of bunches: " << s.bunch_list.size() << endl;
  		cout << "The number of bunches (factorized):" << s.stats_fact_bunches << endl;
  		cout << "The number of calls (vol): " << s.stats_vol_calls << endl;
  		cout << "The number of vol reuses: " << s.stats_vol_reuses << endl;
  		cout << "The average dims for each call (vol): " << (double)s.stats_total_dims / s.stats_vol_calls << endl;
	}
	
  	cout << endl << "====================================" << endl << endl;
  	if (latte) cout << "The total count (LattE): " << total_latte << endl;
  	if (vinci) cout << "The total volume (Vinci): " << total_vinci << endl;
  	if (polyvest) cout << "The total volume (PolyVest): " << total_polyvest << endl;
  	cout << endl << "====================================" << endl << endl;
  	
  	//print to output
  	ofstream fout(output_file);
   	if (latte) fout << total_latte << endl;
  	if (vinci) fout << total_vinci << endl;
  	if (polyvest) fout << total_polyvest << endl;
  	fout.close();	
	
	//////////////////////////////////////////////////////////////////////


	return 1;

}
