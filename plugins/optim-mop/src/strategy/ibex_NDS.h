/*
 * ibex_NDS.h
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_IntervalVector.h"
#include "ibex_pyPlotter.h"
#include <map>
#include <list>
#include "ibex_BxpMOPData.h"
#include "ibex_OptimizerMOP.h"

#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_

using namespace std;

namespace ibex {


/**
 * comparation function for sorting NDS2 by increasing x and decreasing by y
 */
struct sorty2{
	bool operator()(const Vector& y1, const Vector& y2){
		if(y1[0] != y2[0])
			return y1[0]<y2[0];
		return y1[1]>y2[1];

	}
};

class NDS_data{
public:
	NDS_data() : x1(1), x2(1){ n=0; }

	NDS_data(const Vector& x1) : x1(x1),  x2(1) { n=1;}

	NDS_data(const Vector& x1, const Vector& x2) : x1(x1),  x2(x2) { n=2; }

/*
	NDS_data(const NDS_data& other) {
		x1=NULL; x2=NULL;
		if(other.x1) x1=new Vector(*other.x1);
		if(other.x2) x2=new Vector(*other.x2);
	}
*/

	/*NDS_data& operator=(NDS_data& other){
		if(this != &other){
			//if(x1) delete x1;
			//if(x2) delete x2;
			x1=NULL; x2=NULL;
			if(other.x1) x1=new Vector(*other.x1);
			if(other.x2) x2=new Vector(*other.x2);
		}
		return *this;
	}*/

	~NDS_data(){
		//if(x1) delete x1;
		//if(x2) delete x2;
	}

  int n; //number of valid vectors
	Vector x1;
	Vector x2;
};

/**
 * \brief Segment based non-dominated set
 */
class NDS_seg {
public:


	NDS_seg() {
		clear();
	};

	virtual ~NDS_seg() { };


	void clear(){
		NDS2.clear();
		//the first point
		Vector v(2); v[0]=NEG_INFINITY; v[1]=POS_INFINITY;
		NDS2.insert(make_pair(v, NDS_data()));
		//the middle point
		v[0]=POS_INFINITY;
		NDS2.insert(make_pair(v, NDS_data()));
		//the last point
		v[1]=NEG_INFINITY;
		NDS2.insert(make_pair(v, NDS_data()));
	}

	Vector lb(){
		if(NDS2.size()>3){
			auto it=NDS2.begin();
			it++;
			double min1=it->first[0];
			it=NDS2.end();
			it--;
			it--;
			double min2=it->first[1];
			Vector v(2); v[0]=min1; v[1]=min2;
			return v;
		}
		else{
			Vector v(2); v[0]=POS_INFINITY; v[1]=POS_INFINITY;
			return v;
		}
	}

	Vector nadir(){
		if(NDS2.size()>3){
			auto it=NDS2.begin();
			it++;it++;
			double min2=it->first[1];
			it=NDS2.end();
			it--;
			it--;it--;
			double min1=it->first[0];
			Vector v(2); v[0]=min1; v[1]=min2;
			return v;
		}
		else{
			Vector v(2); v[0]=NEG_INFINITY; v[1]=NEG_INFINITY;
			return v;
		}
	}

	Interval hypervolume(const Interval& y1, const Interval& y2) const{
		Interval hv=0.0;
		double prev1=y1.lb();
		double prev2=y2.ub();
		for(auto ndp:NDS2){
			if(ndp.first[0] == NEG_INFINITY) continue;
			if(ndp.first[1] == NEG_INFINITY) continue;

			double next1=std::min(ndp.first[0],y1.ub());
			double next2=std::min(y2.ub(),ndp.first[1]);

			if(next1 > prev1){
				Interval hvv=(Interval(next1)-prev1)*( (y2.ub()	-Interval(prev2))  +  (Interval(prev2)-next2)/2.0);
				hv+=hvv;

			}

			if(next1>=prev1 && next2<=prev2){
				prev1=next1;
				prev2=next2;
			}


		}
		return hv;
	}

	int size() const{
		return NDS2.size();
	}


	static IntervalVector get_box_y(const Cell* c){
		IntervalVector boxy(2);
		int n=c->box.size()-2;
		boxy[0]=c->box[n];
		boxy[1]=c->box[n+1];
		return boxy;
	}

	/**
	 * \brief return true if new_p is dominated by the NDS
	 */
	bool is_dominated(const Vector& new_p);

	void addPoint(IntervalVector& y, const NDS_data& data=NDS_data()){
		addPoint(y.lb(), data);
	};


	/**
	* Add a point in the NDS structure
	*/
	void addPoint(const Vector& y, const NDS_data& data=NDS_data());

	/**
	* Add a segment in the NDS structure. Return 0 if the segment did not modify the NDS
	*/
	bool addSegment(const pair< Vector, Vector>& y1y2, const NDS_data& data=NDS_data());



	// Agrega el lowerbound de una caja al NDS
	pair <Vector, Vector> add_lb(Cell& c);

    struct NoIntersectionException : public exception {
       const char * what () const throw () {
          return "NoIntersectionException";
       }
    };

	//returns a list of points non-dominated by lb
	list< Vector > non_dominated_points(const Vector& lb, bool cutting_points=true){
		list < Vector > inpoints;


		//first potential dominated point
		Vector v(2); v[0]=lb[0]; v[1]=NEG_INFINITY;
		map< Vector, NDS_data >:: iterator ent1=NDS2.upper_bound(v);
		Vector v11 = ent1->first; ent1--;
		//last point before x=lbx
		Vector v10 = ent1->first;

		Vector v1(2); v1[0]=lb[0]; v1[1]=v11[1];
		Vector v2(2); v2[0]=lb[0]; v2[1]=v10[1];

		//x-point cutting lbx
		if(cutting_points){
			  Vector firstp(2);
		    firstp = pointIntersection( v10, v11, v1, v2);
		   	if(firstp[1] <= lb[1] ){
		   		return inpoints; //empty list
		   	}
		   inpoints.push_back(firstp);
	  }

		ent1++;
		//points dominated by lb
		while(ent1->first[1] > lb[1]){
			inpoints.push_back(ent1->first);
			ent1++;
		}

		if(cutting_points){
	    Vector lastp(2);
			//last potential dominated point
			ent1--;
			Vector v21 = ent1->first;
			//first point after y=lby
			ent1++;
			Vector v20 = ent1->first;

			//x-point cutting lby
		    v1[0]=v20[0]; v1[1]=lb[1];
		    v2[0]=v21[0]; v2[1]=lb[1];
			lastp = pointIntersection( v20, v21, v1, v2);
			inpoints.push_back(lastp);
	 }

		return inpoints;
	}

  /**
	* \brief Returns the point intersecting two segments. Otherwise it throw
	* a NoIntersectionException
	* It is conservative, that is:
	* 1) if there are intersection it should return a point dominated by the real intersection
	* 2) if there are no intersection it may return an exception or a point dominated by one segment
	* TODO: revise with colinear generated examples
	*/
	static Vector pointIntersection(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3);



  static bool _trace;
	double distance(const Cell* c){
		IntervalVector box_y=get_box_y(c);
		double a = ((BxpMOPData*) c->prop[BxpMOPData::id])->a;
		double w_lb = ((BxpMOPData*) c->prop[BxpMOPData::id])->w_lb;


		double dist;
		IntervalVector points(4);
		if(a!=0){
			pair<Vector, Vector> points = get_segment(box_y.lb(),-1/a, w_lb/a);
			dist= distance(points.first, points.second, -1/a, w_lb/a);
		}else{
			pair<Vector, Vector> points = get_segment(box_y.lb());
			dist= distance(points.first, points.second);
		}

		return std::max(0.0,dist);
	}

	//returns the segment yA--yB of the line y_2=m*y_2+c dominated by lb
    pair <Vector, Vector> get_segment(const Vector& lb, double m=POS_INFINITY, double c=POS_INFINITY){
		double max_dist=NEG_INFINITY;

		Interval Ay=lb[1];
		Interval Bx=lb[0];

		if(m!=POS_INFINITY){
			Ay = Interval(m)*lb[0]+c;  // Ax=lbx
			Bx = (Interval(lb[1])-c)/m; // By=lby
			if(Ay.lb() < lb[1]){ Ay=lb[1]; }
			if(Bx.lb() < lb[0]){ Bx=lb[0]; }
		}

		Vector yA(2);
		Vector yB(2);


		yA[0]=lb[0];
		yA[1]=Ay.lb();
		yB[0]=Bx.lb();
		yB[1]=lb[1];
		return make_pair(yA,yB);
	}


	// m in [-oo, 0]
	double distance(Vector& yA, Vector& yB, double m=POS_INFINITY, double c=POS_INFINITY){
		Interval Ax=yA[0];
		Interval Ay=yA[1];
		Interval Bx=yB[0];
		Interval By=yB[1];

		double max_dist=NEG_INFINITY;

		Vector lb(2); lb[0]=yA[0], lb[1]=yB[1];
		list< Vector> inner_segments= non_dominated_points(lb);
		Vector* p0=NULL;

		bool Adist=false;
		bool Bdist=false;

		for(auto p : inner_segments){
			if(p[0]==POS_INFINITY && p[1]==POS_INFINITY) return POS_INFINITY;

			Interval dist;
			//up-left point
			if((p[0]-Ax).lb() < (p[1]-Ay).ub() || p[1]==POS_INFINITY){
				dist=p[0]-Interval(Ax);
			}
			//bottom-right point
			else if((p[1]-By).lb() < (p[0]-Bx).ub() || p[0]==POS_INFINITY){
				dist=p[1]-Interval(By);
				if(!Bdist && p0){
					Interval mm=NEG_INFINITY;
					if(p[0]-(*p0)[0] != 0)
						mm= (Interval(p[1])-(*p0)[1])/(Interval(p[0])-(*p0)[0]);

					if(mm.lb() > -1 && mm.lb() < 0.0 ){
						Interval cc= p[1] - mm*p[0];
						dist=std::max(dist.ub(), ((mm*Ax - Ay + cc)/(1.0-mm)).lb());
						Bdist=true;
					}
				}
			}
			//cy-45-degree zone
			else{
				dist= -(m*p[0] - p[1]+c)/(1.0-m);
				if(!Adist && p0){
					Interval mm=NEG_INFINITY;
					if(p[0]-(*p0)[0] != 0)
						mm= (Interval(p[1])-(*p0)[1])/(Interval(p[0])-(*p0)[0]);

					if(mm.lb() <-1 && mm.lb()>NEG_INFINITY){
						Interval cc= p[1] - mm*p[0];
						//cout << "dist-A:" << (mm*Bx - By + cc)/(1.0-mm) << endl;
						dist=std::max(dist.ub(), ((mm*Bx - By + cc)/(1.0-mm)).lb());
						Adist=true;
					}
				}
			}


			if(dist.ub()>max_dist){

				max_dist=dist.ub();
			}

			if(p0) delete p0;
			p0=new Vector(p);

		}
		if(p0) delete p0;
		return max_dist;
	}

	pair<Vector, NDS_data> get(const Vector& y){
		auto it=NDS2.lower_bound(y);
		it--;
		return *it;
	}


	// The current non-dominated set sorted by increasing y1
    // el valor (vector xl, vector xr) corresponde al final del segmento anterior
	// (resp. comienzo del segmento siguiente) en x.
	map< Vector, NDS_data, sorty2 > NDS2;
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_ */
