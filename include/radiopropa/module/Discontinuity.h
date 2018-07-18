#ifndef DISCONTINUITY_H
#define DISCONTINUITY_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Module.h"
#include "radiopropa/Geometry.h"
#include "radiopropa/Referenced.h"

namespace radiopropa
{
/**
 @class Discontinuity
 @brief Dcsiontinuity in the refractive index at which secondary rays are created according to   Fresnell equations.
 */
class Discontinuity: public Module {
private:
	ref_ptr<Surface> surface;
	double n1, n2;

public:
    Discontinuity(Surface *_surface, double _n1, double _n2); 
    void process(Candidate *candidate) const;
};



}

#endif // DISCONTINUITY_H
