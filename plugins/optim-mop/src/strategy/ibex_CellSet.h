//============================================================================
//                                  I B E X
// File        : ibex_CellSet.h
// Author      : Gilles Chabert, Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Sep 12, 2014
//============================================================================

#ifndef __IBEX_CELL_SET_H__
#define __IBEX_CELL_SET_H__

#include "ibex_Random.h"
#include <set>
#include <map>
#include <queue>
#include "ibex_BxpMOPData.h"
#include "ibex_CellBufferOptim.h"

using namespace std;

namespace ibex {

/**
 * \ingroup optim
 *
 * \brief cell Set buffer (for global optimization)
 *  This is a simple-set buffer where the heap criterion
 * is passed by template
 */

	template<class T>
	class CellSet : public CellBufferOptim {
	public:

	  CellSet();

	  void flush();


	  /** Return the size of the buffer. */
	  unsigned int size() const;

	  /** Return true if the buffer is empty. */
	  bool empty() const;

	  /** push a new cell on the stack. */
	  void push(Cell* cell);

	  /** Pop a cell from the stack and return it.*/
	  Cell* pop();

	  /** Return the next box (but does not pop it).*/
	  Cell* top() const;

	  virtual double minimum() const;

	  virtual void contract(double loup);

	private:
		/* Set of Cells */
		typename std::multiset<Cell*, T> cset;

	};

	/**
	 * This criterion corresponds to the SR1 criterion of the paper of Fernandez&Toth (2007) for bi-objective problems
	 * Authors says that in this way the curve is generated from left-top to right-bottom
	 */
	struct minLB {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-1].lb() != c2->box[n-1].lb()) return (c1->box[n-1].lb() < c2->box[n-1].lb());
		  return (c1->depth < c2->depth);

	  }
	};

	/**
	 * OC1 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * min y1.lb
	 */
	struct OC1 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-1].lb() != c2->box[n-1].lb()) return (c1->box[n-1].lb() < c2->box[n-1].lb());
		  if(c1->box[n-2].lb() != c2->box[n-2].lb()) return (c1->box[n-2].lb() < c2->box[n-2].lb());
		  return (c1->depth < c2->depth);
	  }
	};

	/**
	 * OC2 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * min y2.lb
	 */
	struct OC2 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-2].lb() != c2->box[n-2].lb()) return (c1->box[n-2].lb() < c2->box[n-2].lb());
		  if(c1->box[n-1].lb() != c2->box[n-1].lb()) return (c1->box[n-1].lb() < c2->box[n-1].lb());
		  return (c1->depth < c2->depth);

	  }
	};

	/**
	 * OC3 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * increasing value of y1.lb + y2.lb
	 */
	struct OC3 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-1].lb() + c1->box[n-2].lb() != c2->box[n-1].lb() + c2->box[n-2].lb())
			  return (c1->box[n-1].lb() + c1->box[n-2].lb() < c2->box[n-1].lb() + c2->box[n-2].lb());
		  return (c1->depth < c2->depth);
	  }
	};

	/**
	 * OC4 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * decreasing value of hypervolume of the point y
	 */
	struct OC4 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  double hyper1=(BxpMOPData::y1_init.ub()-c1->box[n-1].lb())*(BxpMOPData::y2_init.ub()-c1->box[n-2].lb());
		  double hyper2=(BxpMOPData::y1_init.ub()-c1->box[n-1].lb())*(BxpMOPData::y2_init.ub()-c1->box[n-2].lb());

		  if(hyper1 != hyper2) return (hyper1 < hyper2);
		  return (c1->depth < c2->depth);
	  }
	};

	/**
	 * Criteria for bi-objective problems used in the paper by Martin et al. (2016)
	 */
	struct weighted_sum {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();
		  double c1_ev= (c1->box[n-2].lb()-BxpMOPData::y1_init.lb())/BxpMOPData::y1_init.diam() +
				  (c1->box[n-1].lb()-BxpMOPData::y2_init.lb())/BxpMOPData::y2_init.diam();

		  double c2_ev= (c2->box[n-2].lb()-BxpMOPData::y1_init.lb())/BxpMOPData::y1_init.diam() +
				  (c2->box[n-1].lb()-BxpMOPData::y2_init.lb())/BxpMOPData::y2_init.diam();

		  return c1_ev < c2_ev;
	  }
	};



	template<class T>
	CellSet<T>::CellSet() {
	}

	template<class T>
	void CellSet<T>::flush() {
		while (!cset.empty()) {
			delete *cset.begin();
			cset.erase(cset.begin());
		}
	}

	template<class T>
	unsigned int CellSet<T>::size() const {
		return cset.size();
	}

	template<class T>
	bool CellSet<T>::empty() const {
		return cset.empty();
	}

	template<class T>
	void CellSet<T>::push(Cell* cell) {
		if (capacity>0 && size() == capacity) throw CellBufferOverflow();
		cset.insert(cell);
	}

	template<class T>
	Cell* CellSet<T>::pop() {
		Cell* c = *cset.begin();
		cset.erase(cset.begin());
		return c;
	}

	template<class T>
	Cell* CellSet<T>::top() const{
		return *cset.begin();
	}

	template<class T>
	double CellSet<T>::minimum() const{
	      if(size()==0)
	        return POS_INFINITY;
	      else
	    	return top()->box[top()->box.size()-1].lb();
	}

	template<class T>
	void CellSet<T>::contract(double loup){
		  typename std::set<Cell*, T>::iterator it= cset.begin();

		  while(it!= cset.end()){
			  if( (*it)->box[(*it)->box.size()-1].lb() > loup ){
				  typename std::set<Cell*, T>::iterator it2=it; it2++;
				  cset.erase(it);
				  it=it2;
			  }else it++;
		  }
	}


}
#endif // __IBEX_CELL_SET_H__
