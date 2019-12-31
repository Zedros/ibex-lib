/*
 * ibex_LoupFinderMOP.h
 *
 *  Created on: 5 oct. 2017
 *      Author: iaraya
 */

#ifndef OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_
#define OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_

#include "ibex_Vector.h"
#include "ibex_System.h"
#include "ibex_LoupFinder.h"
#include "ibex_LinearizerXTaylor.h"
#include "ibex_NormalizedSystem.h"

#include <list>

using namespace std;

namespace ibex {

/**
 * \ingroup optim
 *
 * \brief Upper-bounding algorithm based on XTaylor restriction.
 *
 * The algorithm builds an inner (feasible) polytope inside the
 * current box (see #LinearizerXTaylor) and then minimizes the TWO
 * linear approximation of the goal functions on this polytope via
 * a LP solver. The resulting points are verified a posteriori to
 * be feasible (wrt nonlinear constraint).
 *
 * \note Only works with inequality constraints.
 */
class LoupFinderMOP : public LoupFinder{

public:

	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * The system is an inequality system of constraints.
	 * Goal functions have the form: f1 - z1  and f2 - z2.
	 *
	 * \param sys         - The NLP problem.
	 * \param goal1
	 * \param goal2
	 */
	LoupFinderMOP(const System& sys, const Function& goal1, const Function& goal2, double eqeps=NormalizedSystem::default_eps_h, int nb_sol=50);


	/**
	 * \brief Correct the solution p by using a Hansen feasibility test with eps-inflation
	 */
	bool ub_correction(Vector p, IntervalVector& res);

	/**
	 * \brief Find a candidate solution for the non-dominated set.
	 *
	 * \param box        - the box where the solution is searched
	 * \param loup_point - not used
	 * \param loup       - not used
	 * \return             a candidate solution <x,f(x)> (may be not feasible)
	 * \throws             NotFound in case of failure.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup=POS_INFINITY);


	/**
	 * \brief True if equalities are accepted.
	 *
	 * This function is abstract and may be overriden in the subclass.
	 *
	 * By default: return true.
	 */
	virtual bool rigorous() const {
		return false;
	}

	virtual void clear() {
		lp_solver.clean_ctrs();
		phase=0;
	}

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderMOP() {}


	/**
	 * \brief The relaxed NLP problem for finding feasible points
	 */
	const NormalizedSystem norm_sys;

	/**
	 * Weight of the secondary objective function for the linear program
	 */
	double static _weight2;

protected:

	/**
	 * \brief The real NLP problem.
	 */
	const System& sys;



	/**
	 * \brief Objective function f1
	 * Functions have the form: f1 - z1  and f2 - z2. Thus, in order to
	 * evaluate them we have to set z1 and z2 to [0,0].
	 */
	const Function& goal1;

	/**
	 * \brief Objective function f2
	 */
	const Function& goal2;

	/**
	 * \brief True iff there is an equality.
	 */
	const bool has_equality;

	/** Linearization technique. */
	LinearizerXTaylor lr;

	/** linear solver */
	LPSolver lp_solver;

private:

	//number of solutions to find
	int nb_sol;

	//0: means the first solution of the polytope (or the midpoint)
	//1: means the second solution of the polytope
	//>=2: means the solutions in the line
	int phase;
	Vector vec1; //the first solution of the poltytope
	Vector vec2; //the second solution of the polytope



};

} /* namespace ibex */

#endif /* OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_ */
