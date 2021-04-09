#ifndef CRPROPA_SOURCE_H
#define CRPROPA_SOURCE_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Grid.h"
#include "radiopropa/EmissionMap.h"


#include <vector>

namespace radiopropa {

/**
 @class SourceFeature
 @brief Abstract base class cosmic ray source features
 */
class SourceFeature: public Referenced {
protected:
	std::string description;
public:
	virtual void prepareParticle(ParticleState& particle) const {};
	virtual void prepareCandidate(Candidate& candidate) const;
	std::string getDescription() const;
};


/**
 @class SourceInterface
 @brief Abstract base class for cosmic ray sources
 */
class SourceInterface : public Referenced {
public:
	virtual ref_ptr<Candidate> getCandidate() const = 0;
	virtual std::string getDescription() const = 0;
};

/**
 @class Source
 @brief General cosmic ray source

 This class is a container for source features.
 The source prepares a new candidate by passing it to all its source features
 to be modified accordingly.
 */
class Source: public SourceInterface {
	std::vector<ref_ptr<SourceFeature> > features;
public:
	void add(SourceFeature* feature);
	ref_ptr<Candidate> getCandidate() const;
	std::string getDescription() const;
};

/**
 @class SourceList
 @brief List of cosmic ray sources of individual lumosities.

 The SourceList is a source itself. It can be used if several sources are
 needed in one simulation.
 */
class SourceList: public SourceInterface {
	std::vector<ref_ptr<Source> > sources;
	std::vector<double> cdf;
public:
	void add(Source* source, double weight = 1);
	ref_ptr<Candidate> getCandidate() const;
	std::string getDescription() const;
};


/** @defgroup SourceFeature SourceFeatures
 *  Sourcefeatures are added to sources and manipulate the proeprties of the
 *  emitted candidate.
 *  @{
 */



/**
 @class SourceFrequency
 @brief Sets the initial frequency to a given value
 */
class SourceFrequency: public SourceFeature {
	double E;
public:
	SourceFrequency(double frequency);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceAmplitude
 @brief Sets the initial frequency to a given value
 */
class SourceAmplitude : public SourceFeature {
	double A;
public:
	SourceAmplitude(double amplitude);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};



/**
 @class SourcePosition
 @brief Position of a point source
 */
class SourcePosition: public SourceFeature {
	Vector3d position; /**< Source position */
public:
	SourcePosition(Vector3d position);
	SourcePosition(double d);
	void prepareParticle(ParticleState &state) const;
	void setDescription();
};

/**
 @class SourceMultiplePositions
 @brief Multiple point source positions with individual luminosities
 */
class SourceMultiplePositions: public SourceFeature {
	std::vector<Vector3d> positions;
	std::vector<double> cdf;
public:
	SourceMultiplePositions();
	void add(Vector3d position, double weight = 1);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformSphere
 @brief Uniform random source positions inside a sphere
 */
class SourceUniformSphere: public SourceFeature {
	Vector3d center;
	double radius;
public:
	SourceUniformSphere(Vector3d center, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformShell
 @brief Uniform random source positions on a sphere
 */
class SourceUniformShell: public SourceFeature {
	Vector3d center;
	double radius;
public:
	SourceUniformShell(Vector3d center, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformBox
 @brief Uniform random source positions inside a box
 */
class SourceUniformBox: public SourceFeature {
	Vector3d origin;
	Vector3d size;
public:
	/** Constructor
	 @param origin	lower box corner
	 @param size	upper box corner
	 */
	SourceUniformBox(Vector3d origin, Vector3d size);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceUniformCylinder
 @brief Uniform random source positions inside a Cylinder
 */

class SourceUniformCylinder: public SourceFeature {
	Vector3d origin;
	double height;
	double radius;
public:
	/** Constructor
	 @param origin	lower middle of cylinder
	 @param height	height of the cylinder
	 @param radius	radius of the cylinder
*/
	SourceUniformCylinder(Vector3d origin, double height, double radius);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


/**
 @class SourceUniform1D
 @brief 1D-Positions from a uniform source distribution in an expanding universe

 This source property sets random x-coordinates according to a uniform source
 distribution in a given comoving distance interval.
 This is done by drawing a light travel distance from a flat distribution and
 converting to a comoving distance.
 */
class SourceUniform1D: public SourceFeature {
	double minD; // minimum light travel distance
	double maxD; // maximum light travel distance
	bool withCosmology;
public:
	/** Constructor
	 @param minD	minimum comoving distance
	 @param maxD 	maximum comoving distance
	 @param withCosmology	specify if universe expanding
	 */
	SourceUniform1D(double minD, double maxD, bool withCosmology=true);
	void prepareParticle(ParticleState& particle) const;
	void setDescription();
};

/**
 @class SourceDensityGrid
 @brief Random source positions from a density grid
 */
class SourceDensityGrid: public SourceFeature {
	ref_ptr<ScalarGrid> grid;
public:
	SourceDensityGrid(ref_ptr<ScalarGrid> densityGrid);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceDensityGrid1D
 @brief Random source positions from a 1D density grid
 */
class SourceDensityGrid1D: public SourceFeature {
	ref_ptr<ScalarGrid> grid;
public:
	SourceDensityGrid1D(ref_ptr<ScalarGrid> densityGrid);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceIsotropicEmission
 @brief Isotropic emission from a source
 */
class SourceIsotropicEmission: public SourceFeature {
public:
	SourceIsotropicEmission();
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceDirection
 @brief Emission in a discrete direction
 */
class SourceDirection: public SourceFeature {
	Vector3d direction;
public:
	SourceDirection(Vector3d direction = Vector3d(-1, 0, 0));
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};

/**
 @class SourceEmissionMap
 @brief Deactivate Candidate if it has zero probability in provided EmissionMap
 */
class SourceEmissionMap: public SourceFeature {
	ref_ptr<EmissionMap> emissionMap;
public:
	SourceEmissionMap(EmissionMap *emissionMap);
	void prepareCandidate(Candidate &candidate) const;
	void setEmissionMap(EmissionMap *emissionMap);
	void setDescription();
};

/**
 @class SourceEmissionCone
 @brief Uniform random emission inside a cone
 */
class SourceEmissionCone: public SourceFeature {
	Vector3d direction;
	double aperture;
public:
	SourceEmissionCone(Vector3d direction, double aperture);
	void prepareParticle(ParticleState &particle) const;
	void setDescription();
};


}// namespace radiopropa

#endif // CRPROPA_SOURCE_H
