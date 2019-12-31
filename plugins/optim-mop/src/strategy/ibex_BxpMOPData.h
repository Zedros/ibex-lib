//============================================================================
//                                  I B E X
// File        : ibex_BxpMOPData.h
// Author      : Ignacio Araya
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Oct 18, 2014
// Last Update : Jul 05, 2018
//============================================================================

#ifndef __IBEX_BXP_MOP_DATA_H__
#define __IBEX_BXP_MOP_DATA_H__

#include "ibex_Bxp.h"
#include "ibex_Interval.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_Map.h"

namespace ibex {

/**
 * \ingroup strategy
 *
 * \brief Data required for the Optimizer
 */
class BxpMOPData : public Bxp {
public:
	/**
	 * \brief Constructor for the root node.
	 */
	BxpMOPData();

	/**
	 * \brief Delete *this.
	 */
	~BxpMOPData();

	/**
	 * \brief Create a copy.
	 */
	virtual Bxp* copy(const IntervalVector& box, const BoxProperties& prop) const;

	/**
	 * \brief Update the property upon box modification.
	 *
	 */
	void update(const BoxEvent& event, const BoxProperties& prop);

	//id of the property
	static const long id;

	/**
	 * the current number of generated cells (for instantiating the id)
	 */
	static int nb_cells;

	/**
	 * The evaluation of the objective f1 with the initial box
	 */
	static Interval y1_init;

	/**
	 * The evaluation of the objective f2 with the initial box
	 */
	static Interval y2_init;

    /**unique identifier for comparisons*/
    int idd;

	/** Parameter of the constraint c_y. After filtering we know that y1+a*y2 > w_lb and we can
	 * use this information for filtering**/
	double a;

	/** Parameter of the constraint c_y. After filtering we know that y1+a*y2 > w_lb and we can
	 * use this information for filtering**/
	double w_lb;


	/** Distance of the box to the current non dominated set  */
	double ub_distance;

protected:

	/**
	 * \brief Constructor by copy.
	 */
	explicit BxpMOPData(const BxpMOPData& e);

	static Map<long,false>& ids();
};

/*================================== inline implementations ========================================*/

inline Bxp* BxpMOPData::copy(const IntervalVector& box, const BoxProperties& prop) const {
	return new BxpMOPData(*this);
}

inline void BxpMOPData::update(const BoxEvent& event, const BoxProperties& prop) {
	// Nothing for the moment
}


} // end namespace ibex

#endif // __IBEX_BXP_OPTIM_DATA_H__
