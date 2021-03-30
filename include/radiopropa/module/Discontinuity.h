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

class TransmissiveLayer: public Module {
private:
    ref_ptr<Surface> surface;
    double transmission;
public:
    TransmissiveLayer(Surface *surface, double transmission=1); 
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

class ReflectiveLayer: public Module {
private:
    ref_ptr<Surface> surface;
    double reflection;
    std::map<Candidate*, int> times_reflectedoff;
public:
    ReflectiveLayer(Surface *surface, double reflection=1); 
    void process(Candidate *candidate);
    std::string getDescription() const;

    int getTimesReflectedoff(Candidate *candidate);
};


}

#endif // DISCONTINUITY_H
