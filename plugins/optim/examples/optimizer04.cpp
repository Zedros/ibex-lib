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

	if (argc<8) {
		cerr << "usage: optimizer04 filename filtering linear_relaxation bisection strategy prec goal_prec timelimit "  << endl;
		exit(1);
	}

	System sys(argv[1]);
	string filtering = argv[2];
	string linearrelaxation= argv[3];
	string bisection= argv[4];
	string strategy= argv[5];
	int nbinput=5;
	int beamsize;
	//if (strategy=="bs" || strategy== "beamsearch") {beamsize=atoi(argv[6]); nbinput++;}

	double prec= atof(argv[nbinput+1]);
	double goalprec= atof (argv[nbinput+2]);
	double timelimit = atof(argv[nbinput+3]);
	double eqeps= 1.e-8;

	RNG::srand(atoi(argv[nbinput+4]));

	// the extended system 
	ExtendedSystem ext_sys(sys,eqeps);
	NormalizedSystem norm_sys(sys,eqeps);
	LoupFinderDefault loupfinder (norm_sys,true);
	//LoupFinderDefault loupfinder (norm_sys,false);

	CellBufferOptim* buffer;
	if(strategy=="diving")
		buffer = new CellFeasibleDiving(ext_sys);
	else
		buffer = new CellDoubleHeap  (ext_sys);

	//        cout << "file " << argv[1] << endl;

	// Build the bisection heuristic
	// --------------------------

	Bsc * bs;

	if (bisection=="roundrobin")
	  bs = new RoundRobin (prec);
	else if (bisection== "largestfirst")
          bs= new LargestFirst(prec);
	else if (bisection=="smearsum")
	  bs = new SmearSum(ext_sys,prec);
	else if (bisection=="smearmax")
	  bs = new SmearMax(ext_sys,prec);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(ext_sys,prec);
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(ext_sys,prec);
	else if (bisection=="lsmear")
	  bs = new LSmear(ext_sys,prec);
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// The contractors

	// the first contractor called
	CtcHC4 hc4(ext_sys.ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(ext_sys.ctrs,0.1,true);
	// hc4 inside xnewton loop 
	CtcHC4 hc44xn (ext_sys.ctrs,0.01,false);

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4") 
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4" 
	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(ext_sys,hc44cid,true);
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
	  lr= new LinearizerCombo(ext_sys,LinearizerCombo::ART);
	else if  (linearrelaxation=="compo")
	  lr= new LinearizerCombo(ext_sys,LinearizerCombo::COMPO);
	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (ext_sys, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM_OPP);
	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
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
	Optimizer o(sys.nb_var,*ctcxn,*bs,loupfinder,*buffer,ext_sys.goal_var(),prec,goalprec,goalprec);

	//	cout << " sys.box " << sys.box << endl;

	// the trace 
	o.trace=0;

	// the allowed time for search
	o.timeout=timelimit;

	// the search itself 
	o.optimize(sys.box);

	// printing the results     
	//	o.report();
        cout << o.get_time() << "  " << o.get_nb_cells() << endl;

	//	if (filtering == "acidhc4"  )
	//cout    << " nbcidvar " <<  acidhc4.nbvar_stat() << endl;

	delete bs;
	delete buffer;
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
