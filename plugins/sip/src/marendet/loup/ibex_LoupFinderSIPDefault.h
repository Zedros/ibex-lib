/* ============================================================================
 * I B E X - ibex_LoupFinderSIPDefault.h
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Antoine Marendet
 * Created     : Nov 12, 2018
 * ---------------------------------------------------------------------------- */
 
#ifndef __SIP_IBEX_LOUPFINDERDEFAULT_H__
#define __SIP_IBEX_LOUPFINDERDEFAULT_H__


#include "ibex_SIPSystem.h"

#include "ibex_Cell.h"
#include "ibex_IntervalVector.h"
#include "ibex_LoupFinderSIP.h"

#include <utility>

namespace ibex {
class LoupFinderSIPDefault: public LoupFinderSIP {
public:
	LoupFinderSIPDefault(const SIPSystem& system);
	virtual ~LoupFinderSIPDefault();
	std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);
	std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup, BoxProperties& prop);
};

} // end namespace ibex

#endif // __SIP_IBEX_LOUPFINDERDEFAULT_H__

