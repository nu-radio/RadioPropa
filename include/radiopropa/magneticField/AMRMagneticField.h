#ifndef CRPROPA_AMRMAGNETICFIELD_H
#define CRPROPA_AMRMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_SAGA

#include <iostream> 
#include <string>
#include <cstdio>

#ifdef _OPENMP
    #include "omp.h"
#endif

#include "radiopropa/Units.h"
#include "radiopropa/magneticField/MagneticField.h"
#include "radiopropa/Vector3.h"

#include "saga/LocalProperties.h"
#include "saga/AMRgrid.h"
#include "saga/MagneticField.h"
#include "saga/Referenced.h"



namespace radiopropa {

/**
 @class AMRMagneticField
 @brief Wrapper for saga::MagneticField
 */
class AMRMagneticField: public MagneticField {

private:
	saga::ref_ptr<saga::MagneticField> field;
    double cfLength;
    double cfDensity;
    double cfMagneticField;

public:        
    AMRMagneticField(saga::ref_ptr<saga::MagneticField> field_, double convLength, double convDensity, double convMagneticField)            
    {
        field = field_;
        cfLength = convLength;
        cfDensity = convDensity;
        cfMagneticField = convMagneticField;
    }

    Vector3d getField(const Vector3d &position) const {

        double x = position.x / cfLength;
        double y = position.y / cfLength;
        double z = position.z / cfLength;

        std::vector<double> b = field->getField(x, y, z);
        Vector3d B;
        B.setXYZ(b[0], b[1], b[2]);
        B = B * cfMagneticField;

		return B;   
    }

};

} // namespace radiopropa

#endif // CRPROPA_HAVE_SAGA
#endif // CRPROPA_AMRMAGNETICFIELD_H
