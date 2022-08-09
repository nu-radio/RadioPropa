#ifndef RandomScattering_H
#define RandomScattering_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Module.h"
#include "radiopropa/Geometry.h"
#include "radiopropa/Referenced.h"
#include "radiopropa/Units.h"
#include "radiopropa/Random.h"

namespace radiopropa
{
/**
 @class RandomScattering
 @brief Dcsiontinuity in the refractive index at which secondary rays are created according to   Fresnell equations.
 */
class RandomScattering: public Module {
private:
	ref_ptr<Surface> surface;
	//double fraction;
    //bool surfacemode;
    double tolerance = 0.1*meter;
    double delta;

public:
    RandomScattering(Surface *_surface, double _delta); 
    void process(Candidate *candidate) const;
    std::string getDescription() const;
    
    //double getFraction() const;
    //void setFraction(double new_fraction);
    
    bool parallelToSurface(Vector3d position, Vector3d direction) const;
    bool atSurface(Vector3d position) const;
    bool createdAtSurface(Candidate *candidate) const;
    void positionCorrection(Candidate* candidate, Vector3d new_direction) const;
   // void setSurfacemode(bool mode);
};

class RandomScatteringLayer: public Module {
private:
    ref_ptr<Surface> surface;
    double transmission;
public:
    RandomScatteringLayer(Surface *surface, double transmission=1); 
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};
/*
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

*/
}

#endif // RandomScattering_H
