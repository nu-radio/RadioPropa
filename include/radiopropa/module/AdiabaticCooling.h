#ifndef CRPROPA_ADIABATICCOOLING_H
#define CRPROPA_ADIABATICCOOLING_H

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "radiopropa/Vector3.h"
#include "radiopropa/Module.h"
#include "radiopropa/Units.h"
#include "radiopropa/advectionField/AdvectionField.h"
#include "kiss/logger.h"


namespace radiopropa {

/**
@class AdiabaticCooling
@brief Implements adiabatic cooling/heating due to advection.
*/

class AdiabaticCooling: public Module {
	private:
		ref_ptr<AdvectionField> advectionField;
		double limit;

	public:
	/** Constructor
	@param advectionField 	The advection field used for the adiabatic energy change
	@param limit 		Maximum relative energy change allowed	
*/
		AdiabaticCooling(ref_ptr<AdvectionField> advectionField);
		AdiabaticCooling(ref_ptr<AdvectionField> advectionField, double limit);		
		void process(Candidate *c) const;
	
		void setLimit(double l);

		double getLimit() const;

};	



}; // end namesspace radiopropa
#endif // CRPROPA_ADIABATICCOOLING_H
