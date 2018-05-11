#include "radiopropa/ParticleState.h"
#include "radiopropa/Units.h"
#include "radiopropa/Common.h"

#include <stdlib.h>
#include <sstream>

namespace radiopropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir, Vector3d amplitude) {
	setId(id);
	setFrequency(E);
	setPosition(pos);
	setDirection(dir);
	setAmplitude(amplitude);
}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	// We need to know the rotation of the direction in to update also
	// the amplitude
	Vector3d newDirection = dir.getUnitVector();
	double alpha = direction.getAngleTo(newDirection);
	if (alpha != 0)
	{
		Vector3d rotAxis = direction.cross(newDirection);
		Vector3d newAmplitude = amplitude.getRotated(rotAxis.getUnitVector(), alpha);
		amplitude = newAmplitude;
		direction = newDirection;
	}
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}


void ParticleState::setAmplitude(const Vector3d &newAmplitude) {
	amplitude = newAmplitude; 
}

const Vector3d & ParticleState::getAmplitude() const {
	return amplitude;
}


void ParticleState::setFrequency(double newFrequency) {
	frequency = std::max(0., newFrequency); // prevent negative energies
}

double ParticleState::getFrequency() const {
	return frequency;
}

void ParticleState::setId(int newId) {

}

int ParticleState::getId() const {
	return id;
}


Vector3d ParticleState::getVelocity() const {
	return direction * c_light;
}

Vector3d ParticleState::getMomentum() const {
	return direction * (frequency / c_light);
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << frequency / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction << ", ";
	ss << "A = " << amplitude;
	return ss.str();
}

} // namespace radiopropa
