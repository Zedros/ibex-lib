/*
 * ibex_pyPlotter.cpp
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */

#include "ibex_pyPlotter.h"
#include "ibex_OptimizerMOP.h"

#include <iostream>

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif


#include "ibex_OptimizerMOP.h"

namespace ibex {


void py_Plotter::offline_plot(map< Vector, NDS_data, struct sorty2 >& NDS,
 map< Vector, NDS_data, struct sorty2 >* NDS2, const char* output_file, IntervalVector* focus){
	ofstream output;
	output.open(output_file);

	output << "[";

	map< Vector, NDS_data > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		if(!focus || (*focus).contains(ub->first)){
			//output << "(" << ub->first[0] << "," << ub->first[1] << "),";

      output  << "(" << ub->first[0] << " ; " << ub->first[1] << ")_" <<
                ((ub->second.n==1)? ub->second.x1:0.0) << ",";
    }
	}

  output << "]" << endl;

  if(NDS2){
		output << "[";
		ub=NDS2->begin();
		for(;ub!=NDS2->end();ub++){
			if(!focus || (*focus).contains(ub->first))
				output << "(" << ub->first[0] << "," << ub->first[1] << "),";
		}

	  output << "]" << endl;
  }else
		output << "[]" << endl;

	output.close();

}

} /* namespace ibex */
