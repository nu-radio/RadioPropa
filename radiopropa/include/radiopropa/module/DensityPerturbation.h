#ifndef DENSITYPERTURBATION_H
#define DENSITYPERTURBATION_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Module.h"
#include "radiopropa/Geometry.h"
#include "radiopropa/Referenced.h"
#include "radiopropa/Units.h"
#include "radiopropa/Vector3.h"
#include <iostream>

namespace radiopropa
{
/**
 @class PerturbationLayer
 @brief Perturbation Layer in the density and therefor also refractive index resulting in entrapping rays.
 */
class PerturbationLayer: public Module {
protected:
	ref_ptr<Surface> surface;
	double thickness, threshold, fraction;

public:
    PerturbationLayer(Surface *surface=NULL, double thickness=0, double threshold=0.01, double fraction=0.023); 
    ~PerturbationLayer();
    void process(Candidate* candidate) const;
    std::string getDescription() const;

    ref_ptr<Surface> getSurface() const;

    double getThickness() const;
    void setThickness(double new_thickness);

    double getThreshold() const;
    void setThreshold(double new_threshold);

    double getFraction() const;
    void setFraction(double new_fraction);

    void positionCorrection(Candidate* candidate, Vector3d new_direction) const;
    bool parallelToLayer(Vector3d position, Vector3d direction) const;
    bool inLayer(Vector3d position) const;
    bool createdInLayer(Candidate* candidate) const;
    PerturbationLayer* clone() const;
};

class PerturbationHorizontal: public PerturbationLayer {
protected:
	double z;

public:
	PerturbationHorizontal(double z=0, double thickness=0, double threshold=0.01, double fraction=0.023); 
	~PerturbationHorizontal();
	void positionCorrection(Candidate *candidate, Vector3d new_direction) const;
	PerturbationLayer* clone() const;
};


}

#endif // DENSITYPERTURBATION_H