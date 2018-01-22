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

struct NuclearMassTable {
	std::vector<double> table;

	NuclearMassTable() {
		std::string filename = getDataPath("nuclear_mass.txt");
		std::ifstream infile(filename.c_str());

		if (!infile.good())
      KISS_LOG_WARNING << "radiopropa: could not open file "<< "  filename";
			//throw std::runtime_error("radiopropa: could not open file " + filename);

		int Z, N;
		double mass;
		while (infile.good()) {
			if (infile.peek() != '#') {
				infile >> Z >> N >> mass;
				table.push_back(mass);
			}
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		infile.close();
	}
};

static NuclearMassTable nuclearMassTable;

double nuclearMass(int id) {
	int A = massNumber(id);
	int Z = chargeNumber(id);
	return nuclearMass(A, Z);
}

double nuclearMass(int A, int Z) {
	if ((A < 1) or (A > 56) or (Z < 0) or (Z > 26) or (Z > A))
		throw std::runtime_error("nuclearMass: nucleus not found " + kiss::str(A) + ", " + kiss::str(Z));
	int N = A - Z;
	return nuclearMassTable.table[Z * 31 + N];
}

} // namespace radiopropa
