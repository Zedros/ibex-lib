/*
 * ibex_NDS.cpp
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_NDS.h"

namespace ibex {

	 //map< pair <double, double>, IntervalVector, sorty2 > NDS_seg::NDS2;
	 bool NDS_seg::_trace;

 	// Agrega el lowerbound de una caja al NDS
 	pair <Vector, Vector> NDS_seg::add_lb(Cell& c){
 		IntervalVector box_y=get_box_y(&c);
		//cout << box_y << endl;
		//cout << ((BxpMOPData*) c.prop[BxpMOPData::id])->a << "," <<
		//((BxpMOPData*) c.prop[BxpMOPData::id])->w_lb << endl;

 		if(OptimizerMOP::cy_contract_var){
 			pair <Vector, Vector> segment = get_segment(box_y.lb(),
 						-1/((BxpMOPData*) c.prop[BxpMOPData::id])->a,
 						((BxpMOPData*) c.prop[BxpMOPData::id])->w_lb/((BxpMOPData*) c.prop[BxpMOPData::id])->a);

 			addPoint(segment.first);
 			addPoint(segment.second);
 			addSegment(segment);
 			return segment;
 		}else{
 			addPoint(box_y.lb());
 			return make_pair(box_y.lb(), box_y.lb());
 		}

 	}

	bool NDS_seg::is_dominated(const Vector& new_p){
		if(new_p[0] == POS_INFINITY && new_p[1] == POS_INFINITY) return false;

		std::map<Vector, NDS_data >::iterator it1 = --NDS2.lower_bound(new_p);

		Vector point1 = it1->first;
		it1++;
		Vector point2 = it1->first;

		// Se comprueba que no sea dominado por el anterior al lower_bound (puede ocurrir)
		if(point1[0] <= new_p[0] && point1[1] <= new_p[1])
			return true;


		// Se comprueba que no sea dominado por el lower_bound
		if(point2[0] <= new_p[0] && point2[1] <= new_p[1] )
			return true;

		// comprobar que no este dominado por la recta que forma los dos puntos anteriores
		if(!(new_p[0] <= point1[0] && new_p[1] <= point1[1] ) &&
				!(new_p[0] <= point2[0] && new_p[1] <= point2[1])) {

			//pendiente de los dos puntos
			Interval m = (Interval(point2[1])-Interval(point1[1]))/
					(Interval(point2[0])-Interval(point1[0]));

			// se obtiene el c de la funcion de los puntos
			Interval c = Interval(point1[1]) - m*Interval(point1[0]);

			// se obtiene el c del nuevo punto
			Interval cEval = Interval(new_p[1]) - m*Interval(new_p[0]);

			if(cEval.lb() > c.ub())
				return true;
		}

		return false;
	}



	bool NDS_seg::addSegment(const pair<Vector, Vector>& y1y2, const NDS_data& data) {

		const Vector& y1=y1y2.first;
		const Vector& y2=y1y2.second;

		if(y1y2.first == y1y2.second){
			addPoint(y1y2.first, data);
			return true;
		}

		//se insertan en DS2 todos los puntos que se ubican entre y1 y y2 incluyendo el anterior a y1 y el siguiente a y2
		std::map<Vector, NDS_data >::iterator aux, it1 = --NDS2.lower_bound(y1);


		Interval m = (Interval(y1[1]) - Interval(y2[1]))/
								(Interval(y1[0]) - Interval(y2[0]));
		double c_ub = (Interval(y1[1]) - m*Interval(y1[0])).ub();


		//segmentos en el rango
		list< pair<Vector, NDS_data> > DS2;
		DS2.push_back(*it1);
		it1++;

		for(;it1 != NDS2.end();) {

			DS2.push_back(*it1);
			if(it1->first[1] < y1[1] && it1->first[1] < y2[1]){
				 break;
			 }

			//se elimina el punto si es dominado por el segmento
			if( c_ub < (Interval(it1->first[1]) - m*Interval(it1->first[0])).lb()
			 && y1[0] < it1->first[0] && y2[1] < it1->first[1]){
				aux = it1; ++aux;
				NDS2.erase(it1);
				it1 = aux;

			} else it1++;
		}

		//se intersecta el segmento con los segmentos de la NDS
		//se agregan las intersecciones en NDS
		int intersections=0;
		pair<Vector, NDS_data> prev = DS2.front();
		for(auto next:DS2){
			if(prev.first==next.first) continue; //esto ocurre la primera vez

			try{
				double m2=-1e100;
				if(prev.first[0]-next.first[0]!=0)
				   m2=(prev.first[1]-next.first[1])/(prev.first[0]-next.first[0]);

				Vector point = pointIntersection(prev.first, next.first, y1, y2);
				if(m2<m.mid())
				   NDS2.insert(make_pair(point,prev.second));
				else
				   NDS2.insert(make_pair(point,data));
				intersections++;
			}catch(NoIntersectionException& e) {  }

			prev.first = next.first;
			prev.second = next.second;
		}

		return (intersections>0);

	}

	void NDS_seg::addPoint(const Vector& new_y, const NDS_data& data) {
		if (is_dominated(new_y)) return;

		// Removes from NDS the points dominated by p
		// Then, adds the new point between the corresponding NDS points and adds new ones
		// intersecting the old segments

		std::map<Vector, NDS_data >::iterator it1 = NDS2.lower_bound(new_y); // Se llega al nodo izquierdo del nodo eval, deberia estar mas arriba
		std::map<Vector, NDS_data >::iterator aux;

		Vector first_dom=it1->first;
		it1--;
		Vector last_dom=it1->first;
		bool first=true;
		NDS_data prev_data;
		NDS_data next_data;
		for(;it1 != NDS2.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if(it1->first[1] < new_y[1]) break;

			// comprueba si esta dominado el punto
			if(new_y[0] <= it1->first[0] && new_y[1] <= it1->first[1]) {
				aux = it1;
				++aux;
				if(first) {
					prev_data=it1->second;
					first_dom=it1->first;
				}
				first=false;
				last_dom=it1->first;
				next_data = it1->second;
				NDS2.erase(it1);
				it1 = aux;
			} else ++it1;
		}

		std::map<Vector, NDS_data >::iterator it2 = --NDS2.lower_bound(new_y);
		it1 = it2;
		it2++;

		Vector aux_y(2);

		aux_y[0]=new_y[0]; aux_y[1]=POS_INFINITY;

		try{
			Vector intersection1 = pointIntersection(it1->first, first_dom, new_y, aux_y);
			aux_y[0]=POS_INFINITY; aux_y[1]=new_y[1];
			Vector intersection2 = pointIntersection(last_dom, it2->first, new_y, aux_y);

			// se agregan el punto y los dos obtenidos anteriormente
			NDS2.insert(make_pair(new_y, data));
			NDS2.insert(make_pair(intersection1, NDS_data()));
			NDS2.insert(make_pair(intersection2, next_data));
		}catch(NoIntersectionException& e) {  }

	}


	// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
	// intersect the intersection point may be stored in the floats i_x and i_y.
	Vector NDS_seg::pointIntersection(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3)
	{

		double p0_x=p0[0];
		double p0_y=p0[1];
		double p1_x=p1[0];
		double p1_y=p1[1];
		double p2_x=p2[0];
		double p2_y=p2[1];
		double p3_x=p3[0];
		double p3_y=p3[1];


		Interval i_x, i_y;
		pair<double,double> i ;


		//Cuando el segmento es vertical su pendiente es infinito
		if(p0_x==p1_x && p2_y==p3_y) {
			if( ((p2_x <= p0_x and p0_x <= p3_x) or
					(p3_x <= p0_x and p0_x <= p2_x))
			   && ((p0_y <= p2_y and p2_y <= p1_y) or
					(p1_y <= p2_y and p2_y <= p0_y)) )

				i = make_pair(p0_x,p2_y);

			else throw NoIntersectionException();
		}
		else if (p0_y==p1_y && p2_x==p3_x){

			if( ((p0_x <= p2_x and p2_x <= p1_x) or
					(p1_x <= p2_x and p2_x <= p1_x)) &&
				((p2_y <= p0_y and p0_y <= p3_y) or
					(p3_y <= p0_y and p0_y <= p2_y)))

					i = make_pair(p2_x,p0_y);

			else throw NoIntersectionException();
		}else{

      if(p0_y==POS_INFINITY) p0_y=std::max(p2_y,p3_y);
			if(p0_y==POS_INFINITY) p0_y=std::min(p2_y,p3_y);
			if(p1_y==POS_INFINITY) p1_y=std::max(p2_y,p3_y);
			if(p1_y==POS_INFINITY) p1_y=std::min(p2_y,p3_y);
			if(p0_x==POS_INFINITY) p0_x=std::max(p2_x,p3_x);
			if(p0_x==POS_INFINITY) p0_x=std::min(p2_x,p3_x);
			if(p1_x==POS_INFINITY) p1_x=std::max(p2_x,p3_x);
			if(p1_x==POS_INFINITY) p1_x=std::min(p2_x,p3_x);
			if(p2_y==POS_INFINITY) p2_y=std::max(p0_y,p1_y);
			if(p3_y==POS_INFINITY) p3_y=std::max(p0_y,p1_y);
			if(p2_x==POS_INFINITY) p2_x=std::max(p0_x,p1_x);
			if(p3_x==POS_INFINITY) p3_x=std::max(p0_x,p1_x);

		  Interval r_x, r_y, s_x, s_y, u, t, p_x=p0_x, p_y=p0_y, q_x=p2_x, q_y=p2_y;
	    r_x = p1_x - p_x;     r_y = p1_y - p_y;
	    s_x = p3_x - q_x;     s_y = p3_y - q_y;

			Interval rxs = -s_x * r_y + r_x * s_y;


			if(!rxs.contains(0)){
	    	//u = (-r_y * (p_x - q_x) + r_x * (p_y - q_y)) / rxs;  //(p-q) x r /rxs
 	    	t = ( s_x * (p_y - q_y) - s_y * (p_x - q_x)) / rxs;  //(p-q) x s /rxs
			}else if((-r_y * (p_x - q_x) + r_x * (p_y - q_y)).contains(0)) {
				//colinear
				Interval rxr = (r_x*r_x + r_y*r_y);
				Interval t0= (q_x-p_x)*r_x + (q_y-p_y)*r_y / (r_x*r_x + r_y*r_y); // (q ��� p) �� r / (r �� r)
				Interval t1= t0 + s_x*r_x + s_y*r_y / (r_x*r_x + r_y*r_y); // t0 + s �� r / (r �� r)

				t = Interval(0,1);
				t &= Interval(std::min(t0.lb(),t1.lb()),std::max(t0.ub(),t1.ub()));
			}else //parallel
				throw NoIntersectionException();


	    if (/*u.ub() >= 0 && u.lb() <= 1 &&*/ t.ub() >= 0 && t.lb() <= 1)
	    {
	        // Collision detected
	        i_x = p_x + (t * r_x);
	        i_y = p_y + (t * r_y);

					if (p0_x==p1_x) i_x=p0_x;
					if (p2_x==p3_x) i_x=p2_x;
					if (p0_y==p1_y) i_y=p0_y;
					if (p2_y==p3_y) i_y=p2_y;

					i = make_pair(i_x.ub(),i_y.ub());
					//cout << r_y << endl;
					//cout << "point:" << i.first << "," << i.second << endl;

	    }else{
				throw NoIntersectionException();
			}
		}

		Vector v(2); v[0]=i.first; v[1]=i.second;
		return v;

	}

} /* namespace ibex */
