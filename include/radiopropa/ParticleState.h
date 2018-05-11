#ifndef CRPROPA_PARTICLE_STATE_H
#define CRPROPA_PARTICLE_STATE_H

#include "radiopropa/Vector3.h"

namespace radiopropa {

/**
 @class ParticleState
 @brief State of the particle: ID, frequency, position, direction

 */
class ParticleState {
private:
	int id; ///< particle ID (Particle Data Group numbering scheme)
	double frequency; ///< total frequency
	Vector3d position; ///< position vector in comoving coordinates
	Vector3d direction; ///< unit vector of velocity or momentum
  Vector3d amplitude;

public:
	ParticleState(int id = 0, double frequency = 0,
			Vector3d position = Vector3d(0, 0, 0),
			Vector3d direction = Vector3d(-1, 0, 0),
			Vector3d amplitude = Vector3d(0, 1, 1)
      );

	/// Set position in comoving coordinates
	void setPosition(const Vector3d &pos);
	/// Get position in comoving coordinates
	const Vector3d &getPosition() const;

	/// Set direction unit vector, non unit-vectors are normalized
	void setDirection(const Vector3d &dir);
	/// Get direction unit vector
	const Vector3d &getDirection() const;

	void setAmplitude(const Vector3d &newAmplitude);
	const Vector3d& getAmplitude() const;

	/// Set frequency in [J]
	void setFrequency(double newFrequency);
	/// Get frequency in [J]
	double getFrequency() const;

	/// Set particle ID
	void setId(int);
	/// Get particle ID
	int getId() const;

	std::string getDescription() const;

	/// Velocity: direction times the speed of light in [m/s]
	Vector3d getVelocity() const;
	/// Momentum: direction times frequency divided by the speed of light [kg m/s]
	Vector3d getMomentum() const;
};

} // namespace radiopropa

#endif // CRPROPA_PARTICLE_STATE_H
