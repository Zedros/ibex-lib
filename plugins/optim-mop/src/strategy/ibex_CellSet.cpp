/*
 * ibex_CellSet.cpp
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */


#include "ibex_CellSet.h"

namespace ibex {

	template class CellSet<OC1>;
	template class CellSet<OC2>;
	template class CellSet<OC3>;
	template class CellSet<OC4>;
	template class CellSet<minLB>;
	template class CellSet<weighted_sum>;

}
