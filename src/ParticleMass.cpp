#include "radiopropa/ParticleMass.h"
#include "radiopropa/ParticleID.h"
#include "radiopropa/Common.h"

#include <kiss/convert.h>
#include <kiss/logger.h>

#include <vector>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace radiopropa {

double nuclearMass(int id) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	return nuclearMass(A, Z);
}

double nuclearMass(int A, int Z) {
  KISS_LOG_WARNING << "radiopropa: does not support nuclear masses";
	return 0; 
}

} // namespace radiopropa
