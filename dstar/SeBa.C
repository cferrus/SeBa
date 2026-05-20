//// SeBa:        Binary evolution program SeBa.
////              computes the evolution of a binary given any 
////              initial conditions (M, m, a, e).
////             
//// Output:      in the form of the following files:
////              -init.dat           contains selected initial conditons.
////              -SeBa.data           contains binary evolution histories
////              -binev.data          contains remnant formation information
////
////              in the form of standard output (cerr).
////
////              Initialized parameters include:
////              - mass of the most massice component (primary star)
////              - mass of its binary companion (secondary star)
////              - semi-major axis of the binary system
////              - orbital eccentricity
////
////              routines incuded can be found in double_star.h.
////              The mass function routines are adoped from mkmass.C
////              and are defined in starbase.h
//// 
////              externally visible routines are:
////              -get_random_mass_ratio
////              -get_random_semi_major_axis
////              -get_random_eccentricity
////              The two utilities for the various parameter are:
////              -extract_...._distribution_type_string(....)
////              and
////              -type_string(char*)
////
////              The executable takes initial conditions (see Options)
////              and returns randomized binary parameters.                 
////                 
//// Options:   -M    upper primary mass limit [100[Msun]]
////            -m    lower limit to primary mass [0.1[Msun]]
////            -x    mass function exponent in case of power law [-2.35]
////            -F/f  mass function option: 0) Equal mass
////                                        1) Power-law [default]
////                                        2) Miller & Scalo
////                                        3) Scalo
////                                        4) Kroupa
////            Option -F requires one of the following strings:
////                      (mf_Power_Law, Miller_Scalo, Scalo, Kroupa)
////                   -f requires the appropriate integer (see mkmass.C)
////             -A   maximum semi-major axis limit [1000000[Rsun]]   
////             -a   minimum semi-major axis limit [0] 
////             -y   exponent for a power-law distribution [0] (flat in log)
////             -G/g Semi major axis option: 0) Equal_sma
////                                          1) Power Law [default]
////                                          2) Duquennoy & Mayor (1987)
////                                          3) Raghavan (2010)
////                                          4) Eggleton (1999)
////            Option -G requires one of the following strings:
////                      (Equal_sma, sma_Power_Law, Duquennoy_Mayor, Raghavan, Eggleton)
////                   -g requires appropriate integer (see starbase.h)
////             -E   maximum eccentricity [1] 
////             -e   minimum eccentricity [0] 
////             -v   exponent for a power-law distribution 
////             -U/u eccentricity option: 0) Equal eccentricity
////                                       1) Power Law 
////                                       2) Thermal distribution [default]
////            Option -U requires one of the following strings:
////                      (Equal_ecc, ecc_Power_Law, Thermal_Distribution)
////                   -u requires appropriate integer (see starbase.h)
////             -Q   maximum mass ratio [1]
////             -q   minimum mass ratio [0]
////                    extra option: q_min<0 : q_min=0.1/selected primary mass
////             -w   exponent for a power-law distribution  
////             -P/p mass ratio option: 0) constant mass ratio
////                                       1) Flat_q
////                                       2) Power Law 
////                                       3) Hogeveen (1992)
////            Option -P requires one of the following strings:
////                      (Equal_q, Flat_q, qf_Power_Law, Hogeveen)
////                   -p requires appropriate integer (see starbase.h)
////
////            -I select input file for reading initial conditions.
////               -uses: double_star::dump as input format.  [no default]
////            -R select random initial conditions    [false]
////               with parameters as discribed above.   
////            -n number of binaries to be simulated.  [1]
////               Options: -I all binaries in input file are computed.
////                        -R the number of binaries indicated.
////                        oterwise one binary is simulated with
////                        -M, -m, -a, -e as initial conditions.
////            -N initial ID number of first simulated binary 
////            -T or -t  binary end time. [13500] Myr
////            -s Random seed
////            -z select metallicity of binaries to be simulated. [0.02] Solar
////            -C Initial stellar type primary star [default is main_sequence]
////            -c Initial stellar type secondary star [default is main_sequence]
//   Note:  libnode.a is referenced for the routines which produce the 
//          mass function
//
//	version 1.0	Simon Portegies Zwart, Utrecht, 1992
//                      -First version with class structure
//	version 2.0	Simon Portegies Zwart, Utrecht, 1994
//                      -Coupling to starlab
//	version 3.0	Simon Portegies Zwart, Amsterdam, June 1997
//	version 3.3	Simon Portegies Zwart, Cambridge, March 1999
//      version ...     Simon Portegies Zwart, lost track....
//	version 4.0	Simon Portegies Zwart, Amsterdam, February 2003
//

#include "dyn.h"
#include "double_star.h"
#include "main_sequence.h"
//#include "dstar_to_dyn.h"
//#include "seba.h"
#include <sstream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif


#ifdef TOOLBOX

local bool read_binary_params(ifstream& in, real &m_prim, 
			      real &m_sec, real &sma, real &ecc, real &z) {

    m_prim = 1000;
    m_sec = 1000;
    sma = 0;
    ecc = 0;
    z = 1;
    while (m_prim>100.0 || m_sec>100.0 || ecc<0 || ecc>1 || z < 0.0001 || z > 0.03){   
      if(in.eof())
        return false;
    
      // reading from input file
    //  int id, tpp, tps; 
    //  real time, Rlp, Rls;
    
      in >> sma >> ecc >> m_prim >> m_sec >> z;
      //  in >>id >> time >> sma >> ecc 
      //     >> tpp >> m_prim >> Rlp >> tps >> m_sec >> Rls;
    }
    
    //PRC(m_prim);PRC(m_sec);PRC(sma);PRC(ecc);PRL(z);
  return true;
}

//const char * SeBa_Filename() {return "SeBa.data";)

/*-----------------------------------------------------------------------------
 *  binev  --
 *-----------------------------------------------------------------------------
 */
local bool  evolve_binary(dyn * bi,
                          real start_time, real end_time,
			  bool stop_at_merger_or_disruption,
			  bool stop_at_remnant_formation,
			  ostream& outstream) {


  double_star* ds = dynamic_cast(double_star*,
				 bi->get_starbase());

  //		Setup star from input data.
  real dt, time=start_time;
  ds->evolve_element(time);

  ds->dump(outstream, true);

  if (!bi->is_root() &&
      bi->get_parent()->is_root())

    do {

      //dt = ds->get_evolve_timestep() + cnsts.safety(minimum_timestep);
      dt =
      Starlab::max(ds->get_evolve_timestep(),cnsts.safety(minimum_timestep));

      time = Starlab::min(time+dt, end_time);

      ds->evolve_element(time);

      if (stop_at_merger_or_disruption &&
	  (ds->get_bin_type() == Merged ||
	   ds->get_bin_type() == Disrupted))
	return false;
      if (stop_at_remnant_formation &&
	 (ds->get_primary()->remnant() || ds->get_secondary()->remnant()))
	return false;

    }
    while (time<end_time);

  ds->dump(outstream, true);
  ds->set_star_story(NULL);

  rmtree(bi, false);
  return true;

}

int main(int argc, char ** argv) {

    bool e_flag = false;
    bool R_flag = false;
    bool F_flag = false;
    bool I_flag = false;
    bool O_flag = false;//doesn't work fully as "SeBa.data" is written in the code multiple times
    bool P_flag = false;
    bool U_flag = false;
    bool G_flag = false;

    bool verbose = false;
    bool stop_at_merger_or_disruption = false;
    bool stop_at_remnant_formation = false;
    bool random_initialization = false;
    stellar_type primary_type = Main_Sequence;
    stellar_type secondary_type = Main_Sequence;
    char * star_type_string = new char[64];
    binary_type bin_type = Detached;
    real binary_fraction = 1.0;

    int n_init = 0;
    int n = 1;

    real  m_prim;
    real  m_sec;
    real  sma;
    real  ecc = 0;
    real  z = 0;

    real m_tot = 1;
    real r_hm = 1;
    real t_hc = 1;
    real metal = cnsts.parameters(Zsun);


    char *mfc = new char[64];
    mass_function mf = mf_Power_Law;
    real m_min = 0.1;
    real m_max = 100;
    real m_exp = -2.35;    
    char *qfc = new char[64];
    mass_ratio_distribution qf = Flat_q;
    real q_min = 0;
    real q_max = 1;
    real q_exp = 0;
    char *afc = new char[64];
    sma_distribution af = sma_Power_Law;
    real a_min = 0;
    real a_max = 1.e+6; 
    real a_exp = -1;                           
    char *efc = new char[64];
    ecc_distribution ef = Thermal_Distribution;
    real e_min = 0;    // allow detection of constant eccentricity
    real e_max = 1;
    real e_exp;

    real start_time = 0;
    real end_time   = 13500;//35;

    char* input_filename;
    char* output_filename;
    output_filename = "SeBa.data";
    //char output_filename = new char "SeBa.data";

    int input_seed=0; 
    char seedlog[64];
    char paramlog[120];

    //check_help();

//    seba_counters* new_seba_counters = new seba_counters;

    extern char *poptarg;
    int c;
    const char *param_string = "hn:N:VRDSM:m:x:F:f:A:a:y:G:g:E:e:v:U:u:Q:q:T:t:I:O:w:P:p:n:s:z:C:c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case 'h':
                cerr <<
"Usage: " << argv[0] << " [options]\n"
"\n"
"MODES\n"
"  Single binary (default):\n"
"    " << argv[0] << " -M <M1> -m <M2> -a <sma> -e <ecc> -z <Z> -T <Myr>\n"
"\n"
"  From input file (one binary per line: sma ecc M1 M2 Z):\n"
"    " << argv[0] << " -I <file> -T <Myr>\n"
"\n"
"  Random population:\n"
"    " << argv[0] << " -R -n <N> -T <Myr> [distribution flags]\n"
"\n"
"GENERAL\n"
"  -h          Show this help message\n"
"  -V          Verbose: print diagnostics to stderr and run log to stdout\n"
"  -T/-t Myr   End time in Myr                         [default: 13500]\n"
"  -n N        Number of binaries to evolve             [default: 1]\n"
"  -N N        Starting binary identity number          [default: 0]\n"
"  -s seed     Fix random seed for reproducibility      [default: 0 (random)]\n"
"  -z Z        Metallicity (0.0001-0.03)                [default: 0.02 = solar]\n"
"  -I file     Read initial conditions from file (absolute path recommended)\n"
"  -O file     Output filename                          [default: SeBa.data]\n"
"  -R          Generate random initial conditions\n"
"  -D          Stop evolution at merger or disruption\n"
"  -S          Stop evolution at first remnant formation\n"
"\n"
"PRIMARY MASS  (random mode)\n"
"  -M Msun     Upper primary mass limit                 [default: 100]\n"
"  -m Msun     Lower primary mass limit                 [default: 0.1]\n"
"  -x exp      Power-law exponent                       [default: -2.35 (Salpeter)]\n"
"  -f N        Mass function (integer):\n"
"                0 = Equal mass\n"
"                1 = Power-law (default)\n"
"                2 = Miller & Scalo\n"
"                3 = Scalo\n"
"                4 = Kroupa\n"
"  -F name     Mass function by name (mf_Power_Law, Miller_Scalo, Scalo, Kroupa)\n"
"\n"
"SEMI-MAJOR AXIS  (random mode)\n"
"  -A Rsun     Upper SMA limit                          [default: 1e6]\n"
"  -a Rsun     Lower SMA limit                          [default: 0]\n"
"  -y exp      Power-law exponent (0 = flat in log)     [default: 0]\n"
"  -g N        SMA distribution (integer):\n"
"                0 = Equal (fixed)\n"
"                1 = Power-law (default)\n"
"                2 = Duquennoy & Mayor (1987)\n"
"                3 = Raghavan (2010)\n"
"                4 = Eggleton (1999)\n"
"  -G name     SMA distribution by name (Equal_sma, sma_Power_Law,\n"
"              Duquennoy_Mayor, Raghavan, Eggleton)\n"
"\n"
"ECCENTRICITY  (random mode)\n"
"  -E ecc      Upper eccentricity limit                 [default: 1]\n"
"  -e ecc      Lower eccentricity limit                 [default: 0]\n"
"  -v exp      Power-law exponent\n"
"  -u N        Eccentricity distribution (integer):\n"
"                0 = Equal (fixed)\n"
"                1 = Power-law\n"
"                2 = Thermal distribution (default)\n"
"  -U name     Eccentricity distribution by name (Equal_ecc, ecc_Power_Law,\n"
"              Thermal_Distribution)\n"
"\n"
"MASS RATIO  (random mode)\n"
"  -Q q        Upper mass ratio limit                   [default: 1]\n"
"  -q q        Lower mass ratio limit                   [default: 0]\n"
"                (negative value: q_min = 0.1 / primary_mass)\n"
"  -w exp      Power-law exponent                       [default: 0]\n"
"  -p N        Mass ratio distribution (integer):\n"
"                0 = Equal (fixed)\n"
"                1 = Flat q (default)\n"
"                2 = Power-law\n"
"                3 = Hogeveen (1992)\n"
"  -P name     Mass ratio distribution by name (Equal_q, Flat_q,\n"
"              qf_Power_Law, Hogeveen)\n"
"\n"
"STELLAR TYPES  (single binary mode)\n"
"  -C type     Initial primary stellar type   [default: Main_Sequence]\n"
"  -c type     Initial secondary stellar type [default: Main_Sequence]\n"
"\n"
"OUTPUT FORMAT\n"
"  9-column space-separated: id  bin_type  time  sma  ecc  s1_type  m1  s2_type  m2\n"
"  Two rows per binary: T=0 (initial) and T=end (final)\n"
"  Masses in Msun, SMA in Rsun, time in Myr\n"
"\n"
"EXAMPLES\n"
"  ./SeBa -M 2 -m 1 -a 200 -e 0.2 -z 0.001 -T 13500\n"
"  ./SeBa -I /path/to/input.txt -T 12550\n"
"  ./SeBa -R -n 100000 -T 13500 -f 4 -m 0.95 -M 10 -s 42\n";
                exit(0);
            case 'V': verbose = true;
		      break;
            case 'R': random_initialization = true;
		      break;
            case 'D': stop_at_merger_or_disruption = true;
		      break;
            case 'S': stop_at_remnant_formation = true;
		      break;
            case 'M': m_max = atof(poptarg);
		      break;
            case 'm': m_min = atof(poptarg);
		      break;
            case 'x': m_exp = atof(poptarg);
		      break;
	        case 'F': F_flag = true;
		      strcpy(mfc, poptarg);
	              break;
	        case 'f': mf = (mass_function)atoi(poptarg);
	              break;
            case 'A': a_max = atof(poptarg);
		      break;
            case 'a': a_min = atof(poptarg);
		      break;
            case 'y': a_exp = atof(poptarg);
		      break;
	        case 'G': G_flag = true;
		      strcpy(afc, poptarg);
	              break;
	       case 'g': af = (sma_distribution)atoi(poptarg);
	              break;
            case 'E': e_max = atof(poptarg);
		      break;
            case 'e': e_min = atof(poptarg);
		      break;
            case 'v': e_exp = atof(poptarg);
		      break;
	       case 'U': U_flag = true;
		      strcpy(efc, poptarg);
	              break;
	       case 'u': ef = (ecc_distribution)atoi(poptarg);
	              break;
            case 'Q': q_max = atof(poptarg);
		      break;
            case 'q': q_min = atof(poptarg);
		      break;
            case 't': 
            case 'T': end_time = atof(poptarg);
		      break;
            case 'I': I_flag = true;
		      input_filename = poptarg;
		      break;
            case 'O': O_flag = true;
		      output_filename = poptarg;
		      break;
            case 'w': q_exp = atof(poptarg);
		      break;
	       case 'P': P_flag = true;
		      strcpy(qfc, poptarg);
	              break;
	       case 'p': qf = (mass_ratio_distribution)atoi(poptarg);
	              break;
	       case 'n': n = atoi(poptarg);
	              break;
	       case 'N': n_init = atoi(poptarg);
	              break;
	       case 's': input_seed = atoi(poptarg);
		      break;
           case 'z': metal = atof(poptarg);
                break;
           case 'C': strcpy(star_type_string, poptarg);
	           primary_type = extract_stellar_type_string(star_type_string);
                break;
           case 'c': strcpy(star_type_string, poptarg);
	           secondary_type = extract_stellar_type_string(star_type_string);
                break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    int actual_seed = srandinter(input_seed);
    if (verbose) cerr << "random number generator seed = " << actual_seed << endl;
    sprintf(paramlog,
	   "   alpha  = %3.1f\n   lambda = %3.1f\n   beta   = %3.1f\n   gamma  = %4.2f\n   CE_method = %d\n   Jloss_method = %d \n",
	    cnsts.parameters(common_envelope_efficiency),
	    cnsts.parameters(envelope_binding_energy),
	    cnsts.parameters(specific_angular_momentum_loss),
	    cnsts.parameters(dynamic_mass_transfer_gamma),
            cnsts.use_common_envelope_method(),
            cnsts.use_jloss_method()
          );

    if (n <= 0) err_exit("mknodes: N > 0 required!");

    if(F_flag)
	mf = extract_mass_function_type_string(mfc);
    delete mfc;
    if(G_flag)
	af = extract_semimajor_axis_distribution_type_string(afc);
    delete afc;
    if(U_flag)
	ef = extract_eccentricity_distribution_type_string(efc);
    delete efc;
    if(P_flag)
	qf = extract_mass_ratio_distribution_type_string(qfc);
    delete qfc;

    actual_seed = srandinter(input_seed);
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);

    double_star::set_suppress_output(true);
    single_star::set_suppress_output(!verbose);

    static ofstream null_stream("/dev/null");
    if (!verbose) cerr.rdbuf(null_stream.rdbuf());

    ifstream infile(input_filename, ios::in);
    if(I_flag) {
      if (!infile) cerr << "error: couldn't read file "
	                << input_filename <<endl;
      if (verbose) cerr << "Reading input from file "<< input_filename <<endl;
    }

    ofstream outfile(output_filename, ios::out|ios::trunc);
    if (!outfile) {
      cerr << "error: couldn't open output file " << output_filename << endl;
      return 1;
    }

    // Create a root for run-level logging (separate from per-binary trees).
    dyn *log_root = mkdyn(1);
    log_root->log_history(argc, argv);
    log_root->log_comment(seedlog);
    log_root->log_comment(paramlog);
    if (verbose) log_root->print_log_story(cerr);

    if (verbose) print_initial_binary_distributions(m_min, m_max, mf, m_exp,
				       q_min, q_max, qf, q_exp,
				       a_min, a_max, af, a_exp,
				       e_min, e_max, ef, e_exp);

    // Pre-read all binary parameters from input file so the evolution
    // loop can be parallelised — each binary is fully independent.
    struct BinaryInput { real m_prim, m_sec, sma, ecc, z; };
    std::vector<BinaryInput> binary_inputs;
    if (I_flag) {
        BinaryInput bp;
        while (read_binary_params(infile, bp.m_prim, bp.m_sec, bp.sma, bp.ecc, bp.z))
            binary_inputs.push_back(bp);
        n = (int)binary_inputs.size();
    } else if (random_initialization) {
        if (metal < 0.0001 || metal > 0.03) {
            cerr << "Parameters are not within valid range" << endl;
            cerr << "0.0001 <= z <= 0.03" << endl;
            return 0;
        }
    } else {
        if (!(m_max<=100.0 && m_min<=100.0 && e_min>=0 && e_min<=1
              && metal>=0.0001 && metal<=0.03)) {
            cerr << "Parameters are not within valid range" << endl;
            cerr << "0.1 <= M <= 100 " << endl;
            cerr << "0.0001 <= z <= 0.03" << endl;
            cerr << " 0 <= e <= 1" << endl;
            return 0;
        }
        n = 1;
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n; i++) {

      real m_prim_i, m_sec_i, sma_i, ecc_i, z_i;

      if (I_flag) {
          m_prim_i = binary_inputs[i].m_prim;
          m_sec_i  = binary_inputs[i].m_sec;
          sma_i    = binary_inputs[i].sma;
          ecc_i    = binary_inputs[i].ecc;
          z_i      = binary_inputs[i].z;
      } else if (random_initialization) {
          // RNG is global state: serialise generation, parallelise evolution.
          z_i = metal;
#pragma omp critical(rng)
          {
              mkrandom_binary(m_min, m_max, mf, m_exp,
                              q_min, q_max, qf, q_exp,
                              a_min, a_max, af, a_exp,
                              e_min, e_max, ef, e_exp,
                              m_prim_i, m_sec_i, sma_i, ecc_i, z_i);
              while (m_prim_i>100.0 || m_sec_i>100.0 || ecc_i<0 || ecc_i>1)
                  mkrandom_binary(m_min, m_max, mf, m_exp,
                                  q_min, q_max, qf, q_exp,
                                  a_min, a_max, af, a_exp,
                                  e_min, e_max, ef, e_exp,
                                  m_prim_i, m_sec_i, sma_i, ecc_i, z_i);
          }
      } else {
          m_prim_i = m_max; m_sec_i = m_min;
          sma_i = a_min;    ecc_i = e_min; z_i = metal;
      }

      // All tree objects are thread-local — no sharing between iterations.
      dyn *root_i = mkdyn(1);
      root_i->set_mass(1);
      root_i->get_starbase()->set_stellar_evolution_scaling(m_prim_i, r_hm, t_hc);
      dyn *the_binary_i = root_i->get_oldest_daughter();
      add_secondary(the_binary_i, m_sec_i/m_prim_i);
      addstar(root_i, start_time, primary_type, z_i, 0, false, secondary_type);

      double_star *ds_i = new_double_star(the_binary_i, sma_i, ecc_i, start_time,
                                          i + n_init, bin_type);
      ds_i->set_use_hdyn(false);
      ds_i->get_primary()->set_identity(0);
      ds_i->get_secondary()->set_identity(1);

      // Evolve into a thread-local buffer; flush to shared file atomically.
      std::ostringstream oss;
      bool reached_end = evolve_binary(the_binary_i, start_time, end_time,
          stop_at_merger_or_disruption, stop_at_remnant_formation, oss);

#pragma omp critical(output)
      {
          outfile << oss.str();
      }

      if (!reached_end) {
          the_binary_i->get_starbase()->set_star_story(NULL);
          rmtree(the_binary_i, false);
      }

      delete the_binary_i;
      delete root_i;
    }

    log_root->log_history(argc, argv);
    log_root->log_comment(seedlog);
    log_root->log_comment(paramlog);
    if (verbose) log_root->print_log_story(cout);
    return 0;
}

#endif
