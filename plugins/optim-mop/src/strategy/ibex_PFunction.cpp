/*
 * ibex_PFunction.cpp
 *
 *  Created on: 23 may. 2018
 *      Author: iaraya
 */

#include "ibex_PFunction.h"
#include "ibex_OptimizerMOP.h"

namespace ibex {

bool PFunction::MIN=true;
bool PFunction::MAX=false;

double PFunction::_min_newton_step=0.1; //25% of the current diameter
double PFunction::_min_diam=0.1; //15% of the initial diameter
double PFunction::_eps_opt=1e-7;

PFunction::PFunction(const Function& f1, const Function& f2, const IntervalVector& xa, const IntervalVector& xb):
		f1(f1),f2(f2), xa(xa), xb(xb), last_f(F1) { }

/**
 * convert pf.t to t in inter
 * t = inter.lb() + pf.t*(inter.ub() - inter.lb());
 */
void PFunction::contract_curve(const Interval& t) {
	if(t.is_empty() || t.lb() < 0 || t. ub() > 1) return;

	IntervalVector a = xa+t.lb()*(xb-xa);
	IntervalVector b = xa+t.ub()*(xb-xa);

	xa = a;
	xb = b;
}

Interval PFunction::eval(const Interval& t, const Interval& m, bool minimize, function f) const{
	/*if(last_f==f && evals.find(t.mid()) != evals.end()){
		return evals.at(t.mid());
	}
	last_f=f;*/

	IntervalVector xt = xa+t*(xb-xa);
	Interval result;
	if(f==F1)
		result = OptimizerMOP::eval_goal(f1,xt,  xt.size());
	else if(f==F2)
		result = OptimizerMOP::eval_goal(f2,xt,  xt.size());
	else
		result = OptimizerMOP::eval_goal(f2,xt, xt.size()) - m*OptimizerMOP::eval_goal(f1, xt,  xt.size());

	if( minimize ){
		evals[t.mid()]=-result;
		return -result;
	}	else {
		evals[t.mid()]=result;
		return result;
	}
}

Interval PFunction::deriv(const Interval& t, const Interval& m, bool minimize, function f) const{
	IntervalVector xt = xa+t*(xb-xa);
	Interval result;
	if(f==F1){
		IntervalVector g1 = OptimizerMOP::deriv_goal(f1, xt, xt.size());
		result = g1*(xb-xa);
	}else if(f==F2){
		IntervalVector g2 = OptimizerMOP::deriv_goal(f2, xt, xt.size());
		result = g2*(xb-xa);
	}else{
		IntervalVector g1 = OptimizerMOP::deriv_goal(f1, xt, xt.size());
		IntervalVector g2 = OptimizerMOP::deriv_goal(f2, xt, xt.size());
		result= (g2-m*g1)*(xb-xa);
	}

	if( minimize )	return -result;
	else return result;
}

IntervalVector PFunction::get_point(const Interval& t) const{
	IntervalVector y(2);
	IntervalVector xt = xa+t*(xb-xa);
	y[0]=OptimizerMOP::eval_goal(f1,xt,  xt.size());
	y[1]=OptimizerMOP::eval_goal(f2,xt,  xt.size());
	return y;
}



void PFunction::get_curve_y(std::vector< pair <double, double> >& curve_y ){
	Interval a;
	IntervalVector interVector = IntervalVector(2);
	int max_iterations = 100;
	for(int i=0;i <= max_iterations; i++) {
		a = Interval((double) i/max_iterations);
		//cout << a << endl;
		interVector = get_point(a);
		//cout << interVector[0].ub() << "," <<  interVector[1].ub() << endl;
		curve_y.push_back(make_pair(interVector[0].ub(), interVector[1].ub()));
	}
}

bool PFunction::newton_rcontract(const Interval& m, bool minimize, function f, Interval& inter, const Interval& derivate, double lb){
  if(derivate.lb()>=0) return false;

	double point_t = inter.ub();
	Interval point_c = eval(point_t, m, minimize, f);
	double t_before = NEG_INFINITY;
	while(t_before - point_t > _min_newton_step*inter.diam() && point_t > inter.lb()) {

		t_before = point_t;

		point_t = (t_before - (lb - point_c)/derivate.lb()).ub();
		point_c = eval(point_t, m, minimize, f).lb();
	}

	// TODO: no deber��a pasar! problema con derivada? no conservativo?
	if(point_t > inter.lb() and point_c.lb() > lb ){
		cout << "r error: point_c > lb+epsilon " << point_c << "<" << lb << endl;
		return true;
	}

	// Se elimina el intervalo ya que no contiene una solucion mejor a lb+epsilon
	if(point_t < inter.lb())
		return false;
	 else
		inter = Interval(inter.lb(), std::min(point_t,inter.ub()));


	return true;
}

bool PFunction::newton_lcontract(const Interval& m, bool minimize, function f, Interval& inter, const Interval& derivate, double lb){
  if(derivate.ub()<=0) return false;

	double point_t = inter.lb();
	Interval point_c = eval(point_t, m, minimize, f);
	double t_before = NEG_INFINITY;

/*  cout << "derivate:" << derivate << endl;
	cout << "inter:" << inter << endl;
	cout << "lb:" << lb << endl;
	cout << "t:" << point_t << endl;
	cout << "c:" << point_c.lb() << endl;*/

	while(point_t - t_before > _min_newton_step*inter.diam() && point_t < inter.ub()) {
		t_before = point_t;
		point_t = ((lb - point_c)/derivate.ub() + t_before).lb();
		point_c = eval(point_t, m, minimize, f);
		//cout << "t:" << point_t << endl;
		//cout << "c:" << point_c.lb() << endl;
	}

		// TODO: no deber��a pasar! problema con derivada? no conservativo?
	if(point_t < inter.ub() and point_c.lb() > lb ){
		cout << "l error: point_c > lb " << lb-point_c.lb() << endl;
    return true;
		//exit(0);
	}

	// Se elimina el intervalo ya que no contiene una solucion mejor a lb+epsilon
	if(point_t > inter.ub())
		return false;
	 else {
		//se contracta el intervalo si point_t > 0
		inter = Interval(std::max(point_t,inter.lb()), inter.ub());
	}

	return true;
}


Vector* PFunction::find_feasible(Vector& y, double eps){
	stack<Interval> T;
	evals.clear();
	T.push(Interval(0,1));

	while(!T.empty()) {
		Interval t = T.top(); T.pop();
		double tmid = t.mid();

    if( eval(tmid, Interval(), false, F1).ub() < y[0] + eps && eval(tmid, Interval(), false, F2).ub() < y[1] + eps)
			return new Vector(get_point(tmid).mid());
		else{
			Interval a = eval(t, Interval(), false, F1);
			Interval b = eval(t, Interval(), false, F2);
			if(a.lb() > y[0] + eps || b.lb() > y[1] + eps) continue;
			if(t.diam()>eps){
				T.push(Interval(t.lb(), tmid));
				T.push(Interval(tmid, t.ub()));
			}
		}
	}
	return NULL;
}

pair<double, double> PFunction::optimize(const Interval& m, bool minimize, function f, double max_c, Interval init){
	evals.clear();

	double diam0=init.diam();
	//if(init.is_empty()) init=Interval(0,1);

	Interval derivate;

	double t_final;
	double lb = NEG_INFINITY;
	stack<Interval> pila;
	pila.push(init);
	Interval inter, left, right;
	double point_t, point_c, t_before;
	Interval y_r, y_c, y_l;
	double lb_interval;


	// pila
	double t_temp;

	int iter=0;
	while(!pila.empty() and lb < max_c) {
		// cout << pila.size() <<endl;

		inter = pila.top();
		pila.pop();



		/************ lowerbounding ***************/
		//Optimized with a map for keeping the evaluations
		y_r=eval(inter.lb(), m, minimize, f);
		y_c=eval(inter.mid(), m, minimize, f);
		y_l= eval(inter.ub(), m, minimize, f);

		lb_interval = y_r.ub();
		t_temp = inter.lb();
		if(lb_interval < y_c.ub()) {
			lb_interval = y_c.ub();
			t_temp = inter.mid();
		}
		if(lb_interval < y_l.ub()) {
			lb_interval = y_l.ub();
			t_temp = inter.ub();
		}

		if(lb_interval + _eps_opt > lb) {
			t_final = t_temp;
			lb = lb_interval+_eps_opt;
		}

		/*
		if(fabs(lb_interval) < 1 && lb_interval + _eps_opt > lb) {
			t_final = t_temp;
			lb = lb_interval+_eps_opt;
		}
		else if(fabs(lb_interval) >= 1 && lb_interval + fabs(lb_interval)*_eps_opt > lb) {
			t_final = t_temp;
			lb = lb_interval + fabs(lb_interval)*_eps_opt;
		}*/
		/***************************************/

		// Contract if derivate is not empty

		derivate = deriv(inter, m, minimize, f);
		if(!derivate.is_empty()) {

			// contract Newton from left
			if(derivate.ub()!=POS_INFINITY &&  !newton_lcontract(m,  minimize, f, inter, derivate, lb))
				continue;

			if(derivate.lb()!=NEG_INFINITY && !newton_rcontract(m,  minimize, f, inter, derivate, lb))
				continue;

		}


		// bisect interval and push in stack
		if(inter.is_bisectable() and inter.diam() > _min_diam*diam0) {
			pair<Interval,Interval> bsc = inter.bisect(0.5);
			pila.push(bsc.first);
			pila.push(bsc.second);
		}
		iter++;
	}

	if(minimize)
		return make_pair(-lb, t_final);
	else
		return make_pair(lb, t_final);

	//if( (!minimize && m.ub() > 0) || (minimize && m.ub() <= 0) || (minimize && m.is_empty()) ) return make_pair(-lb, t_final);
	//else return make_pair(lb, t_final);
}

} /* namespace ibex */
