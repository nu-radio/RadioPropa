#ifndef DISCONTINUITY_H
#define DISCONTINUITY_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Module.h"
#include "radiopropa/Geometry.h"
#include "radiopropa/Referenced.h"
#include "radiopropa/Units.h"

namespace radiopropa
{
/**
 @class Discontinuity
 @brief Dcsiontinuity in the refractive index at which secondary rays are created according to   Fresnell equations.
 */
class Discontinuity: public Module {
private:
	ref_ptr<Surface> surface;
	double n1, n2, fraction;
    bool surfacemode;
    double tolerance = 0.01*meter;

public:
    Discontinuity(Surface *_surface, double _n1, double _n2, bool surfacemode=false, double fraction=0.023); 
    void process(Candidate *candidate) const;
    std::string getDescription() const;
    
    double getFraction() const;
    void setFraction(double new_fraction);
    
    bool parallelToSurface(Vector3d position, Vector3d direction) const;
    bool atSurface(Vector3d position) const;
    bool createdAtSurface(Candidate *candidate) const;
    void positionCorrection(Candidate* candidate, Vector3d new_direction) const;
    void setSurfacemode(bool mode);
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
