/*
 * ibex_OptimizerMOPserver.cpp
 *
 *  Created on: Sep 25, 2019
 *      Author: iaraya
 */

#include "ibex_Timer.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include <float.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <set>

#include "ibex_OptimizerMOPserver.h"

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif

namespace ibex {

string OptimizerMOP_S::instructions_file="";
string OptimizerMOP_S::output_file="";

Vector rpm_compare::ref(2);

OptimizerMOP_S::OptimizerMOP_S(int n, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
		Mode nds_mode, Mode split_mode, double eps, double rel_eps, double eps_rpm) :
		OptimizerMOP(n, f1, f2, ctc, bsc, buffer, finder, nds_mode, split_mode, eps, rel_eps),
		eps_rpm(eps_rpm), sstatus(STAND_BY_SEARCH), rp(2), x_rpm(n), y_rpm(2){

}


void OptimizerMOP_S::zoom(string instruction, IntervalVector& focus, ifstream& myfile){

	//focus[0]=BxpMOPData::y1_init;
	//focus[1]=BxpMOPData::y2_init;


	double y1_lb,y1_ub,y2_lb,y2_ub;

  myfile >> y1_lb >> y1_ub;
	myfile >> y2_lb >> y2_ub;
	focus[0] = Interval(y1_lb,y1_ub);
	focus[1] = Interval(y2_lb,y2_ub);

	cout << focus << endl;

	if(sstatus == STAND_BY_RPM){
		if(!focus.contains(rp)){
		    rpm_stop();
		}
	}


	if(instruction == "zoom_out"){
		focus[0]=BxpMOPData::y1_init;
		focus[1]=BxpMOPData::y2_init;
		update_focus(focus);
	}

	for(auto cc:paused_cells){
		buffer.push(cc);
		cells.insert(cc);
	}
	paused_cells.clear();

}

void OptimizerMOP_S::get_solution(ifstream& myfile){

	string output_file;
	double y1,y2;
	string y1_str, y2_str;
	myfile >> output_file;
	myfile >> y1_str >> y2_str;
	cout << y1_str << "," << y2_str << endl;
	y1 = stod(y1_str);
	y2 = stod(y2_str);


	Vector y(2); y[0]=y1; y[1]=y2;

	pair<Vector, NDS_data> data = ndsH.get(y);
	//cout << "y:" << data.first << endl;
	if(data.second.n>0) cout << data.second.x1 << endl;
	if(data.second.n>1) cout << data.second.x2 << endl;

	Vector* v=NULL;
	Vector realy(2);
	if(data.second.n==0) {cout << "data.second.n=0, error?" << endl; exit(0); }

	if(data.second.n==1 || data.second.x1 == data.second.x2){
		realy[0]=eval_goal(goal1, data.second.x1, data.second.x1.size()).ub();
		realy[1]=eval_goal(goal2, data.second.x1, data.second.x1.size()).ub();
		if( realy[0] < y[0] + eps && realy[1] < y[1] + eps)
		  v=new Vector(data.second.x1);
	}else{
		PFunction pf(goal1, goal2, data.second.x1, data.second.x2);
		v=pf.find_feasible(y, 1e-8);
	}

	ofstream output, output_tmp;
	output.open(output_file,ios_base::app);
	output_tmp.open("output.tmp");
	if(v){
		realy[0]=eval_goal(goal1, *v, v->size()).ub();
		realy[1]=eval_goal(goal2, *v, v->size()).ub();
    output << y1_str << " " << y2_str << endl;
		output << *v << endl;
		output << realy << endl;
		output_tmp << y1_str << " " << y2_str << endl;
		output_tmp << *v << endl;
		output_tmp << realy << endl;
		cout << y1_str << " " << y2_str << endl;
		cout << *v << endl;
		cout << realy << endl;
		delete v;
	}else {
		output << y1_str << " " << y2_str << endl;
		output << "not found" << endl;
		output_tmp << "not found" << endl;
	}
	output.close();
	output_tmp.close();


}


void OptimizerMOP_S::read_instructions(IntervalVector& focus){
	//se lee el archivo de instrucciones y se elimina
	string line; ifstream myfile;
	myfile.open(instructions_file);
	if (myfile.is_open()){
		string instruction;
		while(myfile >> instruction){
			cout << instruction << endl;
			if(instruction=="zoom_in" || instruction=="zoom_out"){
				zoom(instruction, focus, myfile);
				if(sstatus==STAND_BY_SEARCH) sstatus=SEARCH;
				if(sstatus==STAND_BY_RPM) sstatus=RPM;
			}else if(instruction=="upper_envelope"){
				get_solution(myfile);
			}else if(instruction == "rpm"){
				sstatus = RPM;
				rpm_init(myfile);
			}else if(instruction == "rpm_stop"){
				rpm_stop();
			}else if(instruction=="save"){
				string filename;
				myfile >> filename;
				save_state_in_file(filename);
			}else if(instruction=="pause"){
				 if(sstatus==SEARCH) sstatus=STAND_BY_SEARCH;
				 if(sstatus==RPM) sstatus=STAND_BY_RPM;
			}else if(instruction=="continue"){
				 if(sstatus==STAND_BY_SEARCH) sstatus=SEARCH;
				 if(sstatus==STAND_BY_RPM) sstatus=RPM;
			}else if(instruction=="finish"){
				 sstatus=FINISHED;
			}
		}

		myfile.close();
		rename(instructions_file.c_str(), (instructions_file+".old").c_str());
	}

}

void OptimizerMOP_S::write_envelope(IntervalVector& focus){
	//escritura de archivos
	//dormir 1 segundo y lectura de instrucciones
	cout << "escritura de archivo" << endl;

	NDS_seg LBaux;
	NDS_seg UBaux=ndsH;

	update_focus(focus);

	for(auto cc:cells)	LBaux.add_lb(*cc);
	for(auto cc:paused_cells) LBaux.add_lb(*cc);
	for(auto cc:rpm_cells) LBaux.add_lb(*cc);

	//se escribe el archivo de salida
	IntervalVector focus2(2);
	focus2[0]=BxpMOPData::y1_init;
	focus2[1]=BxpMOPData::y2_init;
	update_focus(focus2);
	cout << 3 << endl;
	py_Plotter::offline_plot(UBaux.NDS2,  &LBaux.NDS2, output_file.c_str(), &focus2);
}

void OptimizerMOP_S::write_status(double rel_prec){
	ofstream output;
	output.open( (output_file+".state").c_str());
	switch(sstatus){
		case STAND_BY_SEARCH: cout << "STAND_BY" << endl;  output << "STAND_BY" ; break;
		case STAND_BY_RPM: cout << "STAND_BY_RPM" << endl; output << "STAND_BY_RPM" ; break;
		case REACHED_PRECISION: cout << "REACHED_PRECISION" << endl; output << "REACHED_PRECISION" ; break;
		case SEARCH: cout << "SEARCH" << endl; output << "SEARCH" ; break;
		case RPM: cout << "RPM" << endl; output << "RPM" ; break;
		case FINISHED: cout << "FINISHED" << endl; output << "FINISHED" ; break;
	}
	output << "," << rel_prec << endl;
	output.close();
}

void OptimizerMOP_S::rpm_stop(){
	cout << "rpm_stop" << endl;
	for(auto cc:rpm_cells){
			buffer.push(cc);
			cells.insert(cc);
	}
	rpm_cells.clear();

	for(auto cc:paused_cells){
			buffer.push(cc);
			cells.insert(cc);
	}
	paused_cells.clear();

	if(sstatus == RPM) sstatus=SEARCH;
	else if(sstatus == STAND_BY_RPM) sstatus=STAND_BY_SEARCH;
	cout << "rpm_stop_end" << endl;
}

double OptimizerMOP_S::rpm_init(ifstream& myfile){
	myfile >> rpm_file;

	myfile >> rp[0];
	myfile >> rp[1];
	cout << rp << endl;

	rpm_compare::ref = rp;
	cout << rp << endl;
	ub_rpm = POS_INFINITY;
	// we compute the best point in Y'
	list< Vector > inner_segments = ndsH.non_dominated_points(rp, false);

	for (auto point:inner_segments){
		 double d = rpm_compare::distance(point,rp);

		 if (d<ub_rpm) ub_rpm = d;
		 //x_rpm = point.second.x1;
	}
	cout << ub_rpm << endl;

	// Boxes in cells + cells_pause are revised
	// boxes with lb < ub_rpm-eps_rpm are put into cells_rpm
	rpm_cells.clear();
	for(auto cc:paused_cells){
		if(rpm_compare::distance(get_boxy(cc->box,n).lb(),rp) < ub_rpm-eps_rpm){
			rpm_cells.insert(cc);
		  paused_cells.erase(cc);
		}
	}

	while(!buffer.empty()){
		Cell* cc=buffer.pop();
		if(rpm_compare::distance(get_boxy(cc->box,n).lb(),rp) < ub_rpm-eps_rpm)
			rpm_cells.insert(cc);
		else paused_cells.insert(cc);
	}
	cells.clear();
}

void OptimizerMOP_S::update_focus(IntervalVector& focus){

	IntervalVector new_focus(2);
	new_focus.set_empty();

	for(auto cc:cells){
		IntervalVector boxy=get_boxy(cc->box,n);
		if(new_focus.is_empty())
			new_focus=boxy;
		else new_focus|=boxy;
	}

	for(auto cc:paused_cells){
		IntervalVector boxy=get_boxy(cc->box,n);
		if(new_focus.is_empty())
			new_focus=boxy;
		else new_focus|=boxy;
	}

	focus&=new_focus;

}

OptimizerMOP_S::Status OptimizerMOP_S::optimize(const IntervalVector& init_box, string filename) {

	nb_cells=0;
	buffer.flush();
	ndsH.clear();

	y1_ub.first=POS_INFINITY;
	y2_ub.second=POS_INFINITY;
	time=0;

	cells.clear();

	load_state_from_file(filename, init_box);


	for(auto c:cells){
		buffer.push(c);
	}

	IntervalVector focus(2);
	focus[0]=BxpMOPData::y1_init;
	focus[1]=BxpMOPData::y2_init;

	set<Cell*> paused_cells;
	update_focus(focus);
	cout << focus << endl;

  sstatus=SEARCH;
	return _optimize(init_box, focus);

}

OptimizerMOP_S::Status OptimizerMOP_S::optimize(const IntervalVector& init_box) {
	Cell* root=new Cell(IntervalVector(n+2));
	pre_optimize(init_box, root);
	cells.clear();
	cells.insert(root);
  sstatus=SEARCH;

	IntervalVector focus(2);
	focus[0]=BxpMOPData::y1_init;
	focus[1]=BxpMOPData::y2_init;

	return _optimize(init_box, focus);
}

bool OptimizerMOP_S::upper_bounding_rpm(const IntervalVector& box, Vector& rp, double& ub_rpm) {

	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
	IntervalVector x(n);
	finder.clear();

	list< pair <double, double> > points;
	list< pair< pair< double, double> , pair< double, double> > > segments;

	Vector mid=box2.mid();
	if (finder.norm_sys.is_inner(mid)){
		Vector v(2); v[0]=eval_goal(goal1,mid,n).ub(); v[1]=eval_goal(goal2,mid,n).ub();
		ndsH.addPoint(v, NDS_data(mid));
		double rp_d = rpm_compare::distance(v,rp);
		if(rp_d < ub_rpm) {
			ub_rpm=rp_d;
			x_rpm=mid;
			y_rpm=v;
		}
	}


	try{
			x = finder.find(box2,box2,POS_INFINITY).first;
			Vector v(2); v[0]=eval_goal(goal1,x,n).ub(); v[1]=eval_goal(goal2,x,n).ub();
			ndsH.addPoint(v, NDS_data(x.mid()));
			double rp_d = rpm_compare::distance(v,rp);
			if(rp_d < ub_rpm) {
				ub_rpm=rp_d;
				x_rpm=mid;
				y_rpm=v;
			}
	}catch (LoupFinder::NotFound& ) {
		return true;
	}

	return true;

}

OptimizerMOP_S::Status OptimizerMOP_S::_optimize(const IntervalVector& init_box, IntervalVector& focus) {

	Timer timer;
	timer.start();

	Timer timer_stand_by;
	paused_cells.clear();
	rpm_cells.clear();

	double current_precision = POS_INFINITY;

	int iter = 1;
	try {
		bool server_pause=false;
		while (!buffer.empty() || !paused_cells.empty() || !rpm_cells.empty()) {
			if(iter%5==0) server_pause=true;
			while( (buffer.empty() && rpm_cells.empty()) || sstatus==REACHED_PRECISION || sstatus==STAND_BY_SEARCH || server_pause){
        if (sstatus == SEARCH) timer_stand_by.restart();
				if(server_pause) {
			    	cout << "buffer size:" << buffer.size() << endl;
			    	cout << "eps:" << eps << endl;
						write_envelope(focus);
				}
				sleep(2);
				read_instructions(focus);

        if(sstatus == RPM || sstatus==STAND_BY_RPM){
					cout << "initialized: reference point method" << endl;
					cout << rpm_cells.size() << "," << cells.size() << "," << paused_cells.size() << endl;
					if(rpm_cells.empty() && !paused_cells.empty()){
						rpm_stop();
						ofstream output;
						output.open(rpm_file,ios_base::app);
						output << y_rpm[0] << " " << y_rpm[1] << " ";
						output << x_rpm << endl;
						cout << "distance:" << ub_rpm << endl;
						cout << "solution:" << y_rpm[0] << " " << y_rpm[1] << " ;" << x_rpm << endl;
					}
				}

				write_status(current_precision);

				if(sstatus == FINISHED || (buffer.empty() && paused_cells.empty() && rpm_cells.empty() ) || timer_stand_by.get_time()>1000 ){
					 sstatus = FINISHED;
					 write_status(current_precision);
					 exit(0);
				}

				if(buffer.empty() && (sstatus == SEARCH || sstatus == STAND_BY_SEARCH)) sstatus = REACHED_PRECISION;
				server_pause=false; iter++;

			}

      //Iteration of the solver
			Cell* c=NULL;

			if(sstatus == RPM){ // reference point method
			  c = *rpm_cells.begin();
				rpm_cells.erase(rpm_cells.begin());
				cout << "rpm_distance:" << rpm_compare::distance(get_boxy(c->box,n).lb(),rp) << "(" << ub_rpm << ")" << endl;
				if(cdata->ub_distance <= eps) {
					paused_cells.insert(c);
					continue;
				}
			}else{ //normal strategy
				c = buffer.top();
				buffer.pop();
				cells.erase(c);
				current_precision=cdata->ub_distance;

				//update focus & epsilon
				list< Vector > inner_segments = ndsH.non_dominated_points(focus.lb());
				dominance_peeler2(focus,inner_segments);
				if(rel_eps>0.0)	eps=std::max(focus[0].diam(),focus[1].diam())*rel_eps;

				IntervalVector boxy(2); boxy[0]=c->box[n]; boxy[1]=c->box[n+1];
				if(focus[0].ub()<boxy[0].lb() || focus[1].ub()<boxy[1].lb() ){
					paused_cells.insert(c);
					continue;
				}

				cout << cdata->ub_distance << endl;
				if(current_precision <= eps){
					 while(!buffer.empty())
						 paused_cells.insert(buffer.pop());

					 cells.clear();
					 continue;
				}
			}

			nb_cells++;
			iter++;

      //Contraction
			contract_and_bound(*c, init_box);
			if (c->box.is_empty()) {
				delete c;
				continue;
			}

      if(sstatus != RPM){
	      //upper_bounding
				upper_bounding(c->box);
	      //Discarding by using distance and epsilon
				double dist=ndsH.distance(c);
				if(dist < eps){
					if(dist>=0.0) paused_cells.insert(c);
					else delete c;
					continue;
				}
			}else{ //Discarding by reference point distance
			   upper_bounding_rpm(c->box,rp,ub_rpm);
				 double rp_d = ub_rpm-rpm_compare::distance(get_boxy(c->box,n).lb(),rp);
				 if(rp_d < eps_rpm){
					 paused_cells.insert(c);
				   continue;
				 }

				 if(ub_rpm<0.0){
					  paused_cells.insert(c);
					  rpm_stop();
						ofstream output;
						output.open(rpm_file,ios_base::app);
            output << y_rpm[0] << " " << y_rpm[1] << " ";
						output << x_rpm << endl;
						cout << ub_rpm << endl;
						cout << "solution:" << y_rpm[0] << " " << y_rpm[1] << " ;" << x_rpm << endl;
						continue;
				 }
			}

      //Bisection
			pair<Cell*,Cell*> new_cells;
			try {
				new_cells=pair<Cell*,Cell*>(bsc.bisect(*c));
				delete c; // deletes the cell.
			}
			catch (NoBisectableVariableException& ) {
				throw NoBisectableVariableException();
			}


			if(sstatus == RPM){
				rpm_cells.insert(new_cells.first);
				rpm_cells.insert(new_cells.second);
			}else{
				buffer.push(new_cells.first);
				cells.insert(new_cells.first);

				buffer.push(new_cells.second);
				cells.insert(new_cells.second);
			}


			if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
			time = timer.get_time();

		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;

		cout << "timeout" << endl;
	}

	timer.stop();
	time = timer.get_time();

	return status;
}


} /* namespace ibex */
