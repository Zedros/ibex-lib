/* ============================================================================
 * I B E X - Optimizer Tests
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Mar 2, 2012
 * Last update : Oct 17, 2017
 * ---------------------------------------------------------------------------- */

#include "TestOptimizerMOP.h"
#include "ibex.h"


using namespace std;
using namespace ibex;

namespace ibex {

void TestOptimizerMOP::vec_problem01() {
	char* filename = "../benchs/binh.txt";
	System ext_sys(filename);
	SystemFactory factory;
	Variable w;
	Variable a;

	factory.add_var(w);
	factory.add_var(a);

	factory.add_var(ext_sys.args[ext_sys.nb_var-2]);
	factory.add_var(ext_sys.args[ext_sys.nb_var-1]);
	factory.add_ctr(ext_sys.args[ext_sys.nb_var-2] + a * ext_sys.args[ext_sys.nb_var-1] - w = 0);
	OptimizerMop::_server_mode = false;
		//con el flag cy-contract o cy-contract-full
	OptimizerMOP::cy_contract_var= true;

		//con el flag cy-contract
	OptimizerMOP::_cy_upper= true;

	System *_ext_sys = new System(ext_sys, System(fac2));;

	string filtering = "acidhc4";
	string linearrelaxation = "compo";
	string bisection = "largestfirst";
	string strategy = "NDSdist";
	double eps = 0.0001;
	double rel_eps = 0.0;
	double eps_x= 1e-8;
	double timelimit = 3600;
	double eqeps = 1.e-8;
	double rh = 0.1;
	//--ub=ub1
	bool _segments = false;
	//--ub-ub2
	bool _hamburger = true;

	OptimizerMOP::_plot = false;
	int nb_ub_sols = 50;
	OptimizerMOP::_min_ub_dist = 0.1;
	LoupFinderMOP::_weight2 = 0.01;
	bool no_bisect_y  = false;
	OptimizerMOP::_eps_contract = _eps_contract;
	RNG::srand(0);

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

	CellBufferOptim* buffer = new DistanceSortedCellBufferMOP;

	Bsc * bs = new LargestFirst(p);

	Vector p(ext_sys.nb_var, eps_x);

	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);

	Ctc* ctc = &hc43bcidhc4;
	Linearizer* lr = new LinearizerCombo(*_ext_sys,LinearizerCombo::COMPO);
	CtcFixPoint* cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	CtcPolytopeHull* cxn_poly = new CtcPolytopeHull(*lr);
	CtcCompo* cxn_compo = new CtcCompo(*cxn_poly, hc44xn);
	Ctc* ctcxn = new CtcCompo  (*ctc, *cxn);

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	OptimizerMOP o(sys.nb_var,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,
			OptimizerMOP::HAMBURGER,
			OptimizerMOP::MIDPOINT,	eps, rel_eps);

	OptimizerMOP::_rh = rh;
	OptimizerMOP::_server_mode =_server_mode;
	OptimizerMOP::output_file = "output2.txt";
	OptimizerMOP::instructions_file = "instructions.txt";

	o.trace = false;
	o.timeout = timelimit;
	Status stat = o.optimize(ext_sys.box);

	CPPUNIT_ASSERT(stat == SUCCESS);
}

void TestOptimizerMOP::vec_problem02() {
	
	CPPUNIT_ASSERT(true);
}



void TestOptimizerMOP::issue50_1() {
	CPPUNIT_ASSERT(true);
}

void TestOptimizerMOP::issue50_2() {
	CPPUNIT_ASSERT(true);
}

void TestOptimizerMOP::issue50_3() {
	CPPUNIT_ASSERT(true);
}

void TestOptimizerMOP::issue50_4() {
	CPPUNIT_ASSERT(true);
}

void TestOptimizerMOP::unconstrained() {
	CPPUNIT_ASSERT(true);
}

} // end namespace
