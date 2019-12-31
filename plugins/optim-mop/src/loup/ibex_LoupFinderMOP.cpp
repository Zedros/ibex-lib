/*
 * ibex_LoupFinderMOP.cpp
 *
 *  Created on: 5 oct. 2017
 *      Author: iaraya
 */

#include "ibex_LoupFinderMOP.h"
#include "ibex_PdcHansenFeasibility.h"
#include "ibex_FncActivation.h"


namespace ibex {

 double LoupFinderMOP::_weight2=0.0;

//TODO: remove this recipe for the argument of the max number of iterations of the LP solver
LoupFinderMOP::LoupFinderMOP(const System& sys, const Function& goal1, const Function& goal2, double eqeps, int nb_sol) :
		sys(sys), norm_sys(sys,eqeps), lr(norm_sys,LinearizerXTaylor::RESTRICT, LinearizerXTaylor::INF),
		lp_solver(sys.nb_var, std::max((sys.nb_var)*3,100/*LPSolver::default_max_iter*/)),
		goal1(goal1), goal2(goal2), has_equality(false), nb_sol(nb_sol), phase(0), vec1(norm_sys.nb_var), vec2(norm_sys.nb_var) {

	if (sys.nb_ctr>0)
		// ==== check if the system contains equalities ====
		for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
			if (sys.ops[i]==EQ) {
				(bool&) has_equality = true;
				break;
			}
		}

//	nb_simplex=0;
//	diam_simplex=0;
}

bool LoupFinderMOP::ub_correction(Vector p, IntervalVector& res){
    //if(!norm_sys.is_inner(p)) return false;

	if (!has_equality && norm_sys.is_inner(p)){
		res=IntervalVector(p);
		return true;
	}

	//sys is the original system (with equations)

	FncActivation af(sys,p,NormalizedSystem::default_eps_h);

	if (af.image_dim()==0) {
		res=IntervalVector(p);
		return true;
	}
	//cout << "gere" << endl;

	IntervalVector epsbox(p);

	//epsbox.inflate(NormalizedSystem::default_eps_h);
	//PdcHansenFeasibility pdc(af, false);
	// ====================================================
	// solution #2: we call Hansen test in inflating mode.
	PdcHansenFeasibility pdc(af, true);
	// ====================================================
	//cout << epsbox << endl;
	//cout << af.eval(0,epsbox) << endl;
	if (pdc.test(epsbox)==YES) {
		//note: don't call is_inner because it would check again all equalities (which is useless
		// and perhaps wrong as the solution box may violate the relaxed inequality (still, very unlikely))
		bool satisfy_inequalities=true;
		for (int j=0; j<sys.f_ctrs.image_dim(); j++) {
			if (!af.activated()[j]){

				if (((sys.ops[j]==LEQ || sys.ops[j]==LT)
					  && sys.f_ctrs.eval(j,pdc.solution()).ub()>0)
						||
					((sys.ops[j]==GEQ || sys.ops[j]==GT)
					  && sys.f_ctrs.eval(j,pdc.solution()).lb()<0)) {

					/* TODO: && !entailed->original(j)*/
						satisfy_inequalities=false;
						cout << "not corrected" << endl;
						break;
					}
			}
		}
		if (satisfy_inequalities) {
			res = pdc.solution();
			//cout << "corrected" << endl;
			return true;
		}
	}
	//===========================================================
	return false;
}

std::pair<IntervalVector, double> LoupFinderMOP::find(const IntervalVector& box, const IntervalVector& lp, double l) {

	int n=norm_sys.nb_var;

	//phase 0 or 1: call to simplex
    if(phase < nb_sol && phase<=1 && lp_solver.default_max_bound > box.max_diam() ){


		lp_solver.clean_ctrs();
		lp_solver.set_bounds(box);

		IntervalVector box2(box);
		box2.resize(n+2);
		box2[n]=0.0; box2[n+1]=0.0; 
		IntervalVector ig= (phase==0)?
				(goal1.gradient(box2.mid())+ _weight2*goal2.gradient(box2.mid())) :
				(goal2.gradient(box2.mid())+ _weight2*goal1.gradient(box2.mid()));

    if(nb_sol==1)
       ig = goal1.gradient(box2.mid()) + goal2.gradient(box2.mid()) ;

		if (ig.is_empty()){ // unfortunately, at the midpoint the function is not differentiable
			phase = 0;
			throw NotFound(); // not a big deal: wait for another box...
		}

		ig.resize(n);

		Vector g=ig.mid();

		// set the objective coefficient
		// TODO: replace with lp_solver.set_obj(g) when implemented
		for (int j=0; j<n; j++)
			//lp_solver.set_obj(g);
			lp_solver.set_obj_var(j,g[j]);

		int count = lr.linearize(box,lp_solver);

		if (count==-1) {
			lp_solver.clean_ctrs();
			phase = 0;
			throw NotFound();
		}

		LPSolver::Status_Sol stat = lp_solver.solve();

		if (stat == LPSolver::OPTIMAL) {
			//the linear solution is mapped to intervals
			Vector loup_point(lp_solver.get_primal_sol());

			//correct the point
			for(int i=0;i<box.size();i++){
				if(loup_point[i] < box[i].lb())  loup_point[i] = box[i].lb();
				if(loup_point[i] > box[i].ub())  loup_point[i] = box[i].ub();
			}

			if(phase==0) vec1=loup_point;
			if(phase==1) {
				if(vec1==vec2){
				    phase=0;
				    throw NotFound();
				 }
				vec2=loup_point;
			}

			phase ++;
			return make_pair(loup_point, 0.0);

		}else{
			phase = 0;
			throw NotFound();
		}
    }else if(phase>=2 && phase < nb_sol){
    	Vector vec3 = vec1 + (((double)phase-1.0)/((double)nb_sol))*(vec2-vec1);
		for(int i=0;i<box.size();i++){
		    if(vec3[i] < box[i].lb())  vec3[i] = box[i].lb();
			if(vec3[i] > box[i].ub())  vec3[i] = box[i].ub();
		}
		phase ++;
		return make_pair(vec3, 0.0);
    }

    phase=0;
    throw NotFound();

}

} /* namespace ibex */
