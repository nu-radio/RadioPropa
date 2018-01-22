#ifndef CRPROPA_QUIMBYMAGNETICFIELD_H
#define CRPROPA_QUIMBYMAGNETICFIELD_H

#ifdef CRPROPA_HAVE_QUIMBY

#include "radiopropa/Units.h"
#include "radiopropa/magneticField/MagneticField.h"

#include "quimby/MagneticField.h"

#include <stdexcept>
#include <sstream>

namespace radiopropa {

/**
 @class QuimbyMagneticField
 @brief Wrapper for quimby::MagneticField
 */
class QuimbyMagneticField: public MagneticField {
	quimby::ref_ptr<quimby::MagneticField> field;
public:
	QuimbyMagneticField(quimby::ref_ptr<quimby::MagneticField> field) : field(field) {

	}
	QuimbyMagneticField(quimby::MagneticField *field) : field(field) {
	}
	Vector3d getField(const Vector3d &position) const {
		quimby::Vector3f b, r = quimby::Vector3f(position.x, position.y, position.z);
		bool isGood = field->getField(r / kpc, b);
		if (!isGood) {
			std::ostringstream str;
			str << "QuimbyMagneticField: invalid position at " << position;
			throw std::runtime_error(str.str());
		}
		return Vector3d(b.x, b.y, b.z) * gauss;
	}
};
#if 1
/**
 @class QuimbyMagneticFieldAdapter
 @brief Wrapper to use radiopropa::MagneticField in Quimby
 */
class QuimbyMagneticFieldAdapter: public quimby::MagneticField {
	radiopropa::ref_ptr<radiopropa::MagneticField> field;
public:
	QuimbyMagneticFieldAdapter(radiopropa::ref_ptr<radiopropa::MagneticField> field) : field(field) {

	}

	bool getField(const quimby::Vector3f &position, quimby::Vector3f &b) const {
		radiopropa::Vector3d r = radiopropa::Vector3d(position.x, position.y, position.z) * radiopropa::kpc;
		radiopropa::Vector3d B = field->getField(r);
		b = quimby::Vector3f(B.x, B.y, B.z) / gauss;
		return true;
	}
};
#endif

} // namespace radiopropa



#endif // CRPROPA_HAVE_QUIMBY
#endif // CRPROPA_QUIMBYMAGNETICFIELD_H
