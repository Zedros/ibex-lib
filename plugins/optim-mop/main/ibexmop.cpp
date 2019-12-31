//============================================================================
//                                  I B E X
// File        : ibexmop.cpp
// Author      : Ignacio Araya, Matias Campusano, Damir Aliquintui
// Copyright   : PUCV (Chile)
// License     : See the LICENSE file
// Created     : Jan 01, 2017
// Last Update : Jan 01, 2017
//============================================================================


#include "ibex.h"
#include "args.hxx"


#ifndef _IBEX_WITH_OPTIM_MOP_
#error "You need the plugin Optim MOP to run this example."
#endif

const double default_relax_ratio = 0.2;

using namespace std;
using namespace ibex;
int main(int argc, char** argv){


	try {

	if (argc<1) {
		cerr << "usage: optimizer04 filename filtering linear_relaxation bisection strategy prec timelimit "  << endl;
		exit(1);
	}

	args::ArgumentParser parser("********* IbexMop *********.", "Solve a NLBOOP (Minibex file).");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<std::string> _filtering(parser, "string", "the filtering method (default: acidhc4)", {'f', "filt"});
	args::ValueFlag<std::string> _linear_relax(parser, "string", "the linear relaxation method (default: compo)", {"linear-relax"});
	args::ValueFlag<std::string> _bisector(parser, "string", "the bisection method (default: largestfirst)", {'b', "bis"});
	args::ValueFlag<std::string> _strategy(parser, "string", "the search strategy (default: NDSdist)", {'s', "search"});
	args::ValueFlag<double> _eps(parser, "float", "the desired precision (default: 0.01)", {"eps"});
	args::ValueFlag<double> _eps_rpm(parser, "float", "precision of the relative point method (default: 1e-5)", {"eps_rpm"});
	args::ValueFlag<double> _eps_r(parser, "float", "the desired relative precision (default: 0.01)", {"eps_r"});
	args::ValueFlag<double> _timelimit(parser, "float", "timelimit (default: 100)", {'t',"timelimit"});
	args::Flag _cy_contract(parser, "cy-contract", "Contract using the box y+cy, w_ub=+inf.", {"cy-contract"});
	args::Flag _cy_contract_full(parser, "cy-contract", "Contract using the additional constraint cy.", {"cy-contract-full"});
	args::Flag _eps_contract(parser, "eps-contract", "Contract using eps.", {"eps-contract"});
	//args::ValueFlag<int> _nb_ub_sols(parser, "int", "Max number of solutions added by the inner-polytope algorithm (default: 50)", {"N"});
	//args::ValueFlag<double> _weight2(parser, "float", "Weight of the second objective (default: 0.01)", {"w2","weight2"});
	//args::ValueFlag<double> _min_ub_dist(parser, "float", "Min distance between two non dominated points to be considered (default: eps/10)", {"min_ub_dist"});
	args::Flag _nobisecty(parser, "nobisecty", "Do not bisect y variables.", {"no-bisecty"});
	args::ValueFlag<std::string> _upperbounding(parser, "string", "Upper bounding strategy (no|ub1|ub2) (default: no).", {"ub"});
	args::ValueFlag<double> _rh(parser, "float", "Termination criteria for the ub2 algorithm (dist < rh*ini_dist)", {"rh"});
	args::Flag _server_mode(parser, "server", "Server Mode (some options are discativated).",{"server_mode"});
	args::ValueFlag<std::string> _output_file(parser, "string", "Server Output File ", {"server_out"});
	args::ValueFlag<std::string> _instructions_file(parser, "string", "Server Instructions File", {"server_in"});
	args::ValueFlag<std::string> _input_file(parser, "string", "Loading file", {"input_file"});


	args::Flag verbose(parser, "verbose", "Verbose output. Shows the dominance-free set of solutions obtained by the solver.",{'v',"verbose"});
	args::Flag _trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});
	args::Flag _plot(parser, "plot", "Save a file to be plotted by plot.py.", {"plot"});
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

  //original system: 2 objective functions and constraints
	System ext_sys(filename.Get().c_str());

  //for accing the cy envelope constraint
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

	System *_ext_sys;
	if(OptimizerMOP::cy_contract_var)	_ext_sys=new System(ext_sys, System(fac2));
	else _ext_sys =new System(ext_sys);


	//cout << *_ext_sys << endl;

	string filtering = (_filtering)? _filtering.Get() : "acidhc4";
	string linearrelaxation= (_linear_relax)? _linear_relax.Get() : "compo";
	string bisection= (_bisector)? _bisector.Get() : "largestfirst";
	string strategy= (_strategy)? _strategy.Get() : "NDSdist";
	double eps= (_eps)? _eps.Get() : 0.01 ;
	double eps_rpm= (_eps_rpm)? _eps_rpm.Get() : 0.1 ;
	double rel_eps= (_server_mode && !_eps_r)? 0.01: ((_eps_r)? _eps_r.Get() : 0.0 );
	double eps_x= 1e-8 ;
	double timelimit = (_timelimit)? _timelimit.Get() : 100 ;
	double eqeps= 1.e-8;
	double rh=(_rh)? _rh.Get():0.1;
	bool _segments=(_upperbounding && _upperbounding.Get()=="ub1");
	bool _hamburger=(_upperbounding && _upperbounding.Get()=="ub2");

	OptimizerMOP::_plot = _plot;

	int nb_ub_sols = 50 ;
	OptimizerMOP::_min_ub_dist = 0.1;
	LoupFinderMOP::_weight2 = 0.01 ;
	bool no_bisect_y  = _nobisecty;
	OptimizerMOP::_eps_contract = _eps_contract;
	if(_server_mode) OptimizerMOP::_eps_contract = false;

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
	cout << "rel_eps: " << rel_eps << endl;
	cout << "eps_x: " << eps_x << endl;
	cout << "nb_ub_sols: " << nb_ub_sols << endl;
	cout << "min_ub_dist: " << OptimizerMOP::_min_ub_dist << endl;
	cout << "plot: " <<  ((OptimizerMOP::_plot)? "yes":"no") << endl;
	cout << "weight f2: " << LoupFinderMOP::_weight2 << endl;
	cout << "bisect y?: " << ((no_bisect_y)? "no":"yes") << endl;
	cout << "cy_contract?: " << ((OptimizerMOP::cy_contract_var)? "yes":"no") << endl;
	cout << "eps_contract?: " << ((OptimizerMOP::_eps_contract)? "yes":"no") << endl;
	cout << "segments?: " << ((_segments)? "yes":"no") << endl;
	cout << "hamburger?: " << ((_hamburger)? "yes":"no") << endl;


	SystemFactory fac;

	Array<const ExprNode> symbs;
	Array<const ExprSymbol> _x1x2;
	for(int i=0; i<ext_sys.args.size()-2; i++ ){
		const ExprSymbol& x=ExprSymbol::new_(ext_sys.args[i].name);
		fac.add_var(x);
		symbs.add(x);

	}

	for(int j=2; j<ext_sys.nb_ctr; j++ ){
		const ExprNode& e = ext_sys.ctrs[j].f.expr();
		Array<const ExprSymbol> _x1x2;
		for(int i=0;i<ext_sys.args.size()-2;i++) _x1x2.add(ext_sys.ctrs[j].f.args()[i]);

		const ExprNode& new_e = ExprCopy().copy(_x1x2, symbs, e);
		ExprCtr cc(new_e, ext_sys.ctrs[j].op);
		fac.add_ctr(cc);
		//delete &new_e;

	}


	System sys(fac);
	for(int i=0; i<sys.nb_var; i++ )
		sys.box[i] = ext_sys.box[i];


	//cout << sys << endl;

	IntervalVector box = ext_sys.box.mid();
	box[sys.nb_var]=0;
	box[sys.nb_var+1]=0;

	LoupFinderMOP finder(sys, ext_sys.ctrs[0].f, ext_sys.ctrs[1].f, 1e-8, nb_ub_sols);

	CellBufferOptim* buffer;
	if(strategy=="OC3")
	  buffer = new CellSet<OC3>;
	else if(strategy=="OC4")
	  buffer = new CellSet<OC4>;
	else if(strategy=="weighted_sum")
	  buffer = new CellSet<weighted_sum>;
	else if(strategy=="NDSdist")
	  buffer = new DistanceSortedCellBufferMOP;


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
	//else if (bisection=="lsmear")
	//  bs = new LSmear(ext_sys,p);
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

	CtcFixPoint* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          {
		cxn_poly = new CtcPolytopeHull(*lr);

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
	OptimizerMOP* o;
	if(!_server_mode){
		o = new OptimizerMOP(sys.nb_var,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,
					(_hamburger)?  OptimizerMOP::HAMBURGER: (_segments)? OptimizerMOP::SEGMENTS:OptimizerMOP::POINTS,
					OptimizerMOP::MIDPOINT,	eps, rel_eps);
	}else{
		o = new OptimizerMOP_S(sys.nb_var,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,
							(_hamburger)?  OptimizerMOP::HAMBURGER: (_segments)? OptimizerMOP::SEGMENTS:OptimizerMOP::POINTS,
							OptimizerMOP::MIDPOINT,	eps, rel_eps, eps_rpm);
		OptimizerMOP_S::_rh=rh;
		OptimizerMOP_S::output_file= (_output_file)? _output_file.Get():"output2.txt";
		OptimizerMOP_S::instructions_file=(_instructions_file)? _instructions_file.Get():"instructions.txt";
	}


	if(strategy=="NDSdist"){
		dynamic_cast<DistanceSortedCellBufferMOP*>(buffer)->set(o->ndsH);
	}

	// the trace
	o->trace=(_trace)? _trace.Get() : false;

	// the allowed time for search
	o->timeout=timelimit;

	// the search itself
	if(!_input_file)
		o->optimize(ext_sys.box);
	else
	  o->optimize(ext_sys.box, _input_file.Get());


	// printing the results
	o->report(verbose);



/*
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

*/

	return 0;

	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
