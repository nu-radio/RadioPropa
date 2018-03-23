#ifndef CRPROPA_PARTICLE_STATE_H
#define CRPROPA_PARTICLE_STATE_H

#include "radiopropa/Vector3.h"

namespace radiopropa {

/**
 @class ParticleState
 @brief State of the particle: ID, frequency, position, direction

 The ParticleState defines the state of an ultra-high frequency cosmic ray, which
 is assumed to be traveling at the exact speed of light.
 The cosmic ray state is defined by particle ID, frequency and position and
 direction vector.
 For faster lookup mass and charge of the particle are stored as members.
 */
class ParticleState {
private:
	int id; ///< particle ID (Particle Data Group numbering scheme)
	double frequency; ///< total frequency
	Vector3d position; ///< position vector in comoving coordinates
	Vector3d direction; ///< unit vector of velocity or momentum
	double pmass; ///< particle rest mass
	double charge; ///< particle charge

public:
	ParticleState(int id = 0, double frequency = 0,
			Vector3d position = Vector3d(0, 0, 0),
			Vector3d direction = Vector3d(-1, 0, 0));

	/// Set position in comoving coordinates
	void setPosition(const Vector3d &pos);
	/// Get position in comoving coordinates
	const Vector3d &getPosition() const;

	/// Set direction unit vector, non unit-vectors are normalized
	void setDirection(const Vector3d &dir);
	/// Get direction unit vector
	const Vector3d &getDirection() const;

	/// Set frequency in [J]
	void setFrequency(double newFrequency);
	/// Get frequency in [J]
	double getFrequency() const;
	/// Get rigidity defined as E/(Z*e) in [V]
	double getRigidity() const;

	/// Set particle ID
	void setId(int);
	/// Get particle ID
	int getId() const;

	std::string getDescription() const;

	// ======== Helper methods ========

	/// Electrical charge of the particle in [A]
	double getCharge() const;
	/// Mass of the particle in [kg]
	double getMass() const;

	/// Set Lorentz factor and modify the particle's frequency accordingly
	void setLorentzFactor(double gamma);
	/// Get Lorentz factor
	double getLorentzFactor() const;

	/// Velocity: direction times the speed of light in [m/s]
	Vector3d getVelocity() const;
	/// Momentum: direction times frequency divided by the speed of light [kg m/s]
	Vector3d getMomentum() const;
};

} // namespace radiopropa

#endif // CRPROPA_PARTICLE_STATE_H
