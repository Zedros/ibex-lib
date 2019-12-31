//============================================================================
//                                  I B E X
// File        : optimizer04.cpp
// Author      : Gilles Chabert  Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 12, 2012
// Last Update : Jul 12, 2012
//============================================================================


#include "ibex.h"
#include "../main/args.hxx"


#ifndef _IBEX_WITH_OPTIM_
#error "You need the plugin Optim to run this example."
#endif

const double default_relax_ratio = 0.2;

using namespace std;
using namespace ibex;
int main(int argc, char** argv){


	// ------------------------------------------------
	// Parameterized Optimizer (with a system loaded from a file, and choice of contractor, linearization  and bisector)
        // Load a problem to optimize
	// --------------------------
	try {

	if (argc<1) {
		cerr << "usage: optimizer04 filename filtering linear_relaxation bisection strategy prec timelimit "  << endl;
		exit(1);
	}

	args::ArgumentParser parser("********* IbexOpt (defaultoptimizer) *********.", "Solve a Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<std::string> _filtering(parser, "string", "the filtering method", {'f', "filt"});
	args::ValueFlag<std::string> _linear_relax(parser, "string", "the linear relaxation method", {"linear-relax"});
	args::ValueFlag<std::string> _bisector(parser, "string", "the bisection method", {'b', "bis"});
	args::ValueFlag<std::string> _strategy(parser, "string", "the search strategy", {'s', "search"});
	args::ValueFlag<double> _eps(parser, "float", "eps (the precision of the pareto front)", {"eps"});
	args::ValueFlag<double> _epsx(parser, "float", "eps_x (the precision of the x boxes)", {"eps_x"});
	args::ValueFlag<double> _timelimit(parser, "float", "timelimit", {'t',"time"});
	args::Flag _plot(parser, "plot", "Save a python plot.", {"plot"});
	args::Flag _nobisecty(parser, "nobisecty", "Do not bisect y variables.", {"no-bisecty"});
	args::Flag _cy_contract(parser, "cy-contract", "Contract using the box y+cy, w_ub=+inf.", {"cy-contract"});
	args::Flag _cy_contract_full(parser, "cy-contract", "Contract using the box y+cy.", {"cy-contract-full"});
	args::Flag _eps_contract(parser, "eps-contract", "Contract using eps.", {"eps-contract"});
	args::ValueFlag<int> _nb_ub_sols(parser, "int", "Max number of solutions added by the inner-simplex", {"nb_ub_sols"});
	args::ValueFlag<double> _weight2(parser, "float", "Min distance between two non dominated points to be considered (default: 0.01)", {"w2","weight2"});
	args::ValueFlag<double> _min_ub_dist(parser, "float", "Min distance between two non dominated points to be considered (default: eps/10)", {"min_ub_dist"});
	args::Flag _hv(parser, "hv", "Compute the hypervolume", {"hv"});
	args::Flag _trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	//restricciones del sistema original + goal=NULL
	System ext_sys(filename.Get().c_str());

	SystemFactory fac2;

	Variable w;
	Variable a;

	fac2.add_var(w);
	fac2.add_var(a);

	fac2.add_var(ext_sys.args[ext_sys.nb_var-2]);
	fac2.add_var(ext_sys.args[ext_sys.nb_var-1]);
	fac2.add_ctr(ext_sys.args[ext_sys.nb_var-2] + a * ext_sys.args[ext_sys.nb_var-1] - w = 0);


	OptimizerMOP::cy_contract_var= _cy_contract || _cy_contract_full;
	OptimizerMOP::_cy_upper= _cy_contract_full;
  OptimizerMOP::_hv=_hv;

	System *_ext_sys;
	if(OptimizerMOP::cy_contract_var)	_ext_sys=new System(ext_sys, System(fac2));
	else _ext_sys =new System(ext_sys);


	cout << *_ext_sys << endl;

	string filtering = (_filtering)? _filtering.Get() : "acidhc4";
	string linearrelaxation= (_linear_relax)? _linear_relax.Get() : "compo";
	string bisection= (_bisector)? _bisector.Get() : "largestfirst";
	string strategy= (_strategy)? _strategy.Get() : "NDSdist";
	double eps= (_eps)? _eps.Get() : 0.01 ;
	double eps_x= (_epsx)? _epsx.Get() : 1e-8 ;
	double timelimit = (_timelimit)? _timelimit.Get() : 1e4 ;
	double eqeps= 1.e-8;

	OptimizerMOP::_plot = _plot;

	OptimizerMOP::_nb_ub_sols = (_nb_ub_sols)? _nb_ub_sols.Get() : 1 ;
	OptimizerMOP::_min_ub_dist = (_min_ub_dist)? _min_ub_dist.Get() : 0.1;
	LoupFinderMOP::_weight2 = (_weight2)? _weight2.Get() : 0.01 ;
	bool no_bisect_y  = _nobisecty;
	OptimizerMOP::_eps_contract = _eps_contract;

	if(bisection=="largestfirst_noy"){
		bisection="largestfirst";
		no_bisect_y=true;
	}

	RNG::srand(0);

	cout << "Instance: " << argv[1] << endl;
	cout << "Filtering: " << filtering << endl;
	cout << "Linear Relax: " << linearrelaxation << endl;
	cout << "Bisector: " << bisection << endl;
	cout << "Strategy: " << strategy << endl;
	cout << "eps: " << eps << endl;
	cout << "eps_x: " << eps_x << endl;
	cout << "nb_ub_sols: " << OptimizerMOP::_nb_ub_sols << endl;
	cout << "min_ub_dist: " << OptimizerMOP::_min_ub_dist << endl;
	cout << "plot: " <<  ((OptimizerMOP::_plot)? "yes":"no") << endl;
	cout << "weight f2: " << LoupFinderMOP::_weight2 << endl;
	cout << "bisect y?: " << ((no_bisect_y)? "no":"yes") << endl;
	cout << "cy_contract?: " << ((OptimizerMOP::cy_contract_var)? "yes":"no") << endl;

	// the extended system
	// restricciones del sistema original + variables objetivo y restricciones

	//ExtendedSystem ext_sys(sys, eqeps);

	SystemFactory fac;

	for(int i=0; i<ext_sys.nb_var-2; i++ )
		fac.add_var(ext_sys.args[i]);


	for(int j=2; j<ext_sys.nb_ctr; j++ )
		fac.add_ctr(ext_sys.ctrs[j]);


	System sys(fac);
	for(int i=0; i<sys.nb_var; i++ )
		sys.box[i] = ext_sys.box[i];


	cout << sys << endl;

	IntervalVector box = ext_sys.box.mid();
	box[sys.nb_var]=0;
	box[sys.nb_var+1]=0;

	cout << ext_sys.ctrs[0].f.eval(box) << endl;
	cout << ext_sys.ctrs[1].f.eval(box) << endl;

	LoupFinderMOP finder(sys, ext_sys.ctrs[0].f, ext_sys.ctrs[1].f);

	//NormalizedSystem norm_sys(sys,eqeps);
	//LoupFinderDefault loupfinder (norm_sys,true);
	//LoupFinderDefault loupfinder (norm_sys,false);

	CellBufferOptim* buffer;
	if(strategy=="OC1")
	  buffer = new CellSet<OC1>;
	else if(strategy=="OC2")
	  buffer = new CellSet<OC2>;
	else if(strategy=="OC3")
	  buffer = new CellSet<OC3>;
	else if(strategy=="OC4")
	  buffer = new CellSet<OC4>;
	else if(strategy=="weighted_sum")
	  buffer = new CellSet<weighted_sum>;
	else if(strategy=="NDSdist")
	  buffer = new DistanceSortedCellBufferMOP;
	else if(strategy=="NDSsize")
	  buffer = new CellFeasibleDiving<maxsize>(*new CellNSSet);
	else if(strategy=="diving-minlb")
      buffer = new CellFeasibleDiving<minLB>(*new CellSet<minLB>);
	else if(strategy=="diving-weighted_sum")
	  buffer = new CellFeasibleDiving<weighted_sum>(*new CellSet<weighted_sum>);
	else if(strategy=="diving-NDSdist")
	  buffer = new CellFeasibleDiving<max_distance>(*new DistanceSortedCellBufferMOP);
	else if(strategy=="diving-NDSsize")
	  buffer = new CellFeasibleDiving<maxsize>(*new CellNSSet);

	/*else
		buffer = new CellDoubleHeap  (ext_sys);*/

	//        cout << "file " << argv[1] << endl;

	// Build the bisection heuristic
	// --------------------------

	Bsc * bs;

	Vector p(ext_sys.nb_var, eps_x);
	if(no_bisect_y){
		p[ext_sys.nb_var-2]=POS_INFINITY;
		p[ext_sys.nb_var-1]=POS_INFINITY;
	}

	if (bisection=="roundrobin")
	  bs = new RoundRobin (0);
	else if (bisection== "largestfirst"){
          bs= new LargestFirst(p);
	}else if (bisection=="smearsum")
	  bs = new SmearSum(ext_sys,p);
	else if (bisection=="smearmax")
	  bs = new SmearMax(ext_sys,p);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(ext_sys,p);
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(ext_sys,p);
	else if (bisection=="lsmear")
	  bs = new LSmear(ext_sys,p);
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// The contractors

	// the first contractor called
	CtcHC4 hc4(_ext_sys->ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(_ext_sys->ctrs,0.1,true);
	// hc4 inside xnewton loop
	CtcHC4 hc44xn (_ext_sys->ctrs,0.01,false);

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4"
	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(*_ext_sys,hc44cid,true);
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
	CtcCompo hc4acidhc4 (hc4, acidhc4);



	Ctc* ctc;
	if (filtering == "hc4")
	  ctc= &hc4;
	else if
	  (filtering =="acidhc4")
	  ctc= &hc4acidhc4;
	else if
	  (filtering =="3bcidhc4")
	  ctc= &hc43bcidhc4;
	else {cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}

	Linearizer* lr;
	if (linearrelaxation=="art")
	  lr= new LinearizerCombo(*_ext_sys,LinearizerCombo::ART);
	else if  (linearrelaxation=="compo")
	  lr= new LinearizerCombo(*_ext_sys,LinearizerCombo::COMPO);
	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (*_ext_sys, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM_OPP);
	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
	CtcFixPoint* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          {
		cxn_poly = new CtcPolytopeHull(*lr);
		//BitSet bset=BitSet(_ext_sys->nb_var);
		//bset.add(_ext_sys.nb_var-2); //w

		//cxn_poly->set_contracted_vars(bset);
		cxn_compo =new CtcCompo(*cxn_poly, hc44xn);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	  }
	//  the actual contractor  ctc + linear relaxation
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          ctcxn= new CtcCompo  (*ctc, *cxn);
	else
	  ctcxn = ctc;

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	OptimizerMOP o(sys.nb_var,sys.ctrs,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,eps);
	max_distance::UB= &o.get_UB();


	//	cout << " sys.box " << sys.box << endl;

	// the trace
	o.trace=(_trace)? _trace.Get() : false;

	// the allowed time for search
	o.timeout=timelimit;

	// the search itself
	o.optimize(ext_sys.box);

	// printing the results
	o.report(false);
  //cout << o.get_time() << "  " << o.get_nb_cells() << "  " << o.get_nb_sols() << endl;

	//	if (filtering == "acidhc4"  )
	//cout    << " nbcidvar " <<  acidhc4.nbvar_stat() << endl;

	delete bs;
	delete buffer;
	delete _ext_sys;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn") {
		delete lr;
	    delete ctcxn;
	    delete cxn;
	    delete cxn_poly;
	    delete cxn_compo;
	}



	return 0;

	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
