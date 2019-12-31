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


	System ext_sys("../benchs/MOP/nbi.txt");

	SystemFactory fac2;

	Variable w;
	Variable a;

	fac2.add_var(w);
	fac2.add_var(a);

	fac2.add_var(ext_sys.args[ext_sys.nb_var-2]);
	fac2.add_var(ext_sys.args[ext_sys.nb_var-1]);
	fac2.add_ctr(ext_sys.args[ext_sys.nb_var-2] + a * ext_sys.args[ext_sys.nb_var-1] - w = 0);

	System _ext_sys(ext_sys, System(fac2));

	cout << _ext_sys << endl;

	double prec=  1.e-8;
	double timelimit = 10;
	double eqeps= 1.e-8;

	RNG::srand(1);

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

	CellBufferOptim* buffer = new DistanceSortedCellBufferMOP;

	Bsc * bs = new LSmear(ext_sys,1e-8);

	CtcHC4 hc4(_ext_sys.ctrs,0.01,true);

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	OptimizerMOP o(sys.nb_var,sys.ctrs,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, hc4,*bs,*buffer,finder,prec, 1e-8);
	max_distance::UB= &o.get_UB();

  o.trace=1;

  o.y1=Interval(0,10);
  o.y2=Interval(0,10);



	//test 1 (colinear intersecting segments)

  o.insert_lb_segment(point2(4,4), point2(6,2));
  o.plot(NULL); getchar();
  o.insert_lb_segment(point2(5,3), point2(7,1));
  o.plot(NULL); getchar();
	// expected:(0,10),(4,10),(4,4),(6,2),(7,1),(10,1),(10,0)
  o.get_LB().clear();

	//test 1.1

  o.insert_lb_segment(point2(4,4), point2(6,2));
  o.plot(NULL); getchar();
  o.insert_lb_segment(point2(5,3), point2(7,0));
  o.plot(NULL); getchar();

  o.get_LB().clear();

	//test 1.2

	o.insert_lb_segment(point2(4,4), point2(6,2));
	//o.plot(NULL); getchar();
	o.insert_lb_segment(point2(5,3), point2(7,2));
	o.plot(NULL); getchar();
  o.get_LB().clear();

	//test 1.3
	o.insert_lb_segment(point2(2,6), point2(2,2));
	o.insert_lb_segment(point2(2,4), point2(4,1));
	o.plot(NULL); getchar();
	o.get_LB().clear();

	//test 1.4
	o.y1=Interval(-274.1,-273.7);
	o.y2=Interval(75.6,76.2);
	o.insert_lb_segment(point2(-274.0001518221091,76.10519467725033), point2(-273.8918843804995,75.76993943437671));
	o.insert_lb_segment(point2(-273.8918843804995,75.76993943437671), point2(-273.7618407372951,75.76993943437671));
	o.insert_lb_segment(point2(-273.7618407372951,75.76993943437671), point2(-273.7618407372951,75.62017931612018));
	o.plot(NULL); getchar();
	o.insert_lb_segment(point2(-273.9068242045245,75.79869619288252), point2(-273.7604371938824,75.6979795945986));
	//o.insert_lb_segment(point2(-273.7618407372951,75.76993943437671), point2(-273.7618407372951,75.62017931612018));


	o.plot(NULL); getchar();
	o.get_LB().clear();


	  o.y1=Interval(0,10);
	  o.y2=Interval(0,10);
  //test2
	o.insert_lb_segment(point2(4,4), point2(6,2));
  //o.plot(NULL); getchar();
  o.insert_lb_segment(point2(6,2), point2(8,0));
  o.plot(NULL); getchar();
  o.get_LB().clear();
	//expected: [(0,10),(4,10),(4,4),(6,2),(8,0),(10,0),]


	//test3

	o.insert_lb_segment(point2(2,8), point2(4,8));
  //o.plot(NULL); getchar();
  o.insert_lb_segment(point2(4,6), point2(6,6));
  //o.plot(NULL); getchar();
	o.insert_lb_segment(point2(6,4), point2(8,4));
  //o.plot(NULL); getchar();
	o.insert_lb_segment(point2(8,2), point2(10,2));
	//o.plot(NULL); getchar();
	//o.insert_lb_segment(point2(1,10), point2(10,1));
	//o.plot(NULL); getchar();
	o.insert_lb_segment(point2(4,9), point2(10,0));
	//o.plot(NULL); getchar();
	o.insert_lb_segment(point2(0,10), point2(9,4));
	o.plot(NULL); getchar();
	o.get_LB().clear();

  //test4
	o.insert_lb_segment(point2(2,8), point2(4,4));
	//o.plot(NULL); getchar();
	o.insert_lb_segment(point2(3,6), point2(8,2));
	o.plot(NULL); getchar();
	o.get_LB().clear();
	delete bs;
	delete buffer;

	return 0;

	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
