/*
 * ibex_OptimizerMOPserver.h
 *
 *  Created on: Sep 25, 2019
 *      Author: iaraya
 */

#ifndef __IBEX_OPTIMIZERMOP_S_H__
#define __IBEX_OPTIMIZERMOP_S_H__

#include "ibex_OptimizerMOP.h"

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif

namespace ibex {

/**
 * Criteria for bi-objective problems
 */
struct rpm_compare {

  static double distance(const Vector& v, const Vector& rp){
     double d1 = (v[0]-rp[0]);// / BxpMOPData::y1_init.diam();
  	 double d2 = (v[1]-rp[1]);// / BxpMOPData::y2_init.diam();

  	 double dist = std::max(d1, d2);
  	 return dist;
  }

	bool operator() (const Cell* c1, const Cell* c2){
     int n = c1->box.size();
     Vector y1 = OptimizerMOP::get_boxy(c1->box,n).lb();
     Vector y2 = OptimizerMOP::get_boxy(c2->box,n).lb();

     double d1 = distance(y1,ref);
     double d2 = distance(y2,ref);

     return d1 < d2;
	}

  static Vector ref;
};


class OptimizerMOP_S : public OptimizerMOP {
public:


    typedef enum {STAND_BY_SEARCH, STAND_BY_RPM, REACHED_PRECISION, SEARCH, RPM, FINISHED} ServerStatus;

	OptimizerMOP_S(int n, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
			Mode nds_mode=POINTS, Mode split_mode=MIDPOINT, double eps=default_eps, double rel_eps=0.0,
    double eps_rpm=0.0);

	virtual ~OptimizerMOP_S() { }


	/**
	 * \brief Run the optimization.
	 *
	 * \param init_box             The initial box
	 *
	 * \return SUCCESS             if the global minimum (with respect to the precision required) has been found.
	 *                             In particular, at least one feasible point has been found, less than obj_init_bound, and in the
	 *                             time limit.
	 **
	 *         TIMEOUT             if time is out.
	 *
	 */
	virtual Status optimize(const IntervalVector& init_box);

  virtual Status optimize(const IntervalVector& init_box, string filename);

  Status _optimize(const IntervalVector& init_box, IntervalVector& focus) ;

  bool upper_bounding_rpm(const IntervalVector& box, Vector& rp, double& ub_rpm);

    void save_state_in_file(string filename){
     //write object into the file
     // Object to write in file
     ofstream file_obj;

     // Opening file in append mode
     file_obj.open(filename, ios::out);
     file_obj.write((char*) &BxpMOPData::y1_init, sizeof(BxpMOPData::y1_init));
     file_obj.write((char*) &BxpMOPData::y2_init, sizeof(BxpMOPData::y2_init));

     int len = cells.size();
     file_obj.write((char*) &len, sizeof(int));
     for (auto c:cells){
       file_obj.write((char*) &c->box[0], sizeof(c->box[0])*c->box.size());
       file_obj.write((char*) &cdata->a, sizeof(cdata->a));
       file_obj.write((char*) &cdata->w_lb, sizeof(cdata->w_lb));
       file_obj.write((char*) &cdata->ub_distance, sizeof(cdata->ub_distance));
     }

     len = paused_cells.size();
     file_obj.write((char*) &len, sizeof(int));
     for (auto c:paused_cells){
       file_obj.write((char*) &c->box[0], sizeof(c->box[0])*c->box.size());
       file_obj.write((char*) &cdata->a, sizeof(cdata->a));
       file_obj.write((char*) &cdata->w_lb, sizeof(cdata->w_lb));
       file_obj.write((char*) &cdata->ub_distance, sizeof(cdata->ub_distance));
     }

     len = ndsH.size();

     file_obj.write((char*) &len, sizeof(int));
     for (auto elem:ndsH.NDS2){
       file_obj.write((char*) &elem.first[0], sizeof(elem.first[0])*elem.first.size());
       file_obj.write((char*) &elem.second.n, sizeof(elem.second.n));
       if(elem.second.n >= 1) file_obj.write((char*) &elem.second.x1[0], sizeof(elem.second.x1[0])*elem.second.x1.size());
       if(elem.second.n == 2) file_obj.write((char*) &elem.second.x2[0], sizeof(elem.second.x2[0])*elem.second.x2.size());
     }

     for (auto p : ndsH.NDS2) cout << p.first << endl;

     file_obj.close();
  }

  void load_state_from_file(string filename, const IntervalVector& init_box){
    // Object to read from file
    ifstream file_obj;

    // Opening file in input mode
    file_obj.open(filename, ios::in);
    file_obj.read((char*)&BxpMOPData::y1_init, sizeof(BxpMOPData::y1_init));
    file_obj.read((char*)&BxpMOPData::y2_init, sizeof(BxpMOPData::y2_init));

    int len;
    file_obj.read((char*) &len, sizeof(int));
    for (int i=0; i<len; i++){
      Cell* c=new Cell(init_box);
      file_obj.read((char*)&c->box[0], sizeof(c->box[0])*c->box.size());
      c->prop.add(new BxpMOPData());
    	buffer.add_property(c->box, c->prop);
    	bsc.add_property(c->box, c->prop);
    	ctc.add_property(c->box, c->prop);

      file_obj.read((char*) &cdata->a, sizeof(cdata->a));
      file_obj.read((char*) &cdata->w_lb, sizeof(cdata->w_lb));
      file_obj.read((char*) &cdata->ub_distance, sizeof(cdata->ub_distance));
      cells.insert(c);
    }

    file_obj.read((char*) &len, sizeof(int));
    for (int i=0; i<len; i++){
      Cell* c=new Cell(init_box);
      file_obj.read((char*)&c->box[0], sizeof(c->box[0])*c->box.size());
      c->prop.add(new BxpMOPData());
      buffer.add_property(c->box, c->prop);
    	bsc.add_property(c->box, c->prop);
    	ctc.add_property(c->box, c->prop);

      file_obj.read((char*) c->prop[BxpMOPData::id], sizeof(BxpMOPData));
      cells.insert(c);
    }

    file_obj.read((char*) &len, sizeof(int));
    ndsH.NDS2.clear ();

    for (int i=0; i<len; i++){
      Vector y(2);
      int nn;

      file_obj.read((char*) &y[0], sizeof(y[0])*2);
      file_obj.read((char*) &nn, sizeof(int));
      Vector x1(1);
      Vector x2(1); //or n+2?

      if(nn>=1) {x1.resize(n); file_obj.read((char*) &x1[0], sizeof(x1[0])*(n));}
      if(nn==2) {x2.resize(n); file_obj.read((char*) &x2[0], sizeof(x2[0])*(n));}

      ndsH.NDS2[y]=NDS_data(x1,x2);
      ndsH.NDS2[y].n=nn;
    }
    cout.precision(17);
    for (auto p : ndsH.NDS2) cout << p.first << endl;


    file_obj.close();
  }

	void write_status(double rel_prec);

	/*
    * \brief Print the region of solutions
    *
    *
    *
    * Inputs:
    *    \param cells 				   a
    *    \param paused_cells		   a
    *    \param focus 				   a
    */

 	void write_envelope(IntervalVector& focus);


	/*
    * \brief Read the instruccions of work
    *
    * after reading a instruction, this delete it from the file
    *
    * Inputs:
    *    \param cells 				   a
    *    \param paused_cells		   a
    *    \param focus 				   a
    */

	void read_instructions(IntervalVector& focus);


	/*
    * \brief Update the focus of solution
    *
    * This take in count the hull of the region and the found solutions
    *
    * Inputs:
    *    \param cells 				   a
    *    \param paused_cells		   a
    *    \param focus 				   a
    */

	void update_focus(IntervalVector& focus);

	void zoom(string instruction, IntervalVector& focus, ifstream& myfile);

	void get_solution(ifstream& myfile);

  double rpm_init(ifstream& myfile);

  void rpm_stop();

  set<Cell*> cells;
  set<Cell*> paused_cells;
  multiset<Cell*, rpm_compare> rpm_cells;

  //server variables
	static string instructions_file;
	static string output_file;

	ServerStatus sstatus;

  // Precision of the reference point method
  double eps_rpm;
  // Reference point
  Vector rp;
  Vector x_rpm;
  Vector y_rpm;
  double ub_rpm;
  string rpm_file;

};

} /* namespace ibex */

#endif /* PLUGINS_OPTIM_MOP_SRC_STRATEGY_IBEX_OPTIMIZERMOPSERVER_H_ */
