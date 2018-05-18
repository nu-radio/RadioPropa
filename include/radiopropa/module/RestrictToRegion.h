#ifndef CRPROPA_RESTRICTTOREGION_H
#define CRPROPA_RESTRICTTOREGION_H

#include "radiopropa/Referenced.h"
#include "radiopropa/Candidate.h"
#include "radiopropa/Module.h"
#include "radiopropa/Geometry.h"

namespace radiopropa {
/**
 * \addtogroup Condition
 * @{
 */

/**
 @class RestrictToRegion
 @brief Limit Module to region in simulation.

 */
class RestrictToRegion: public Module {
private:
	ref_ptr<Surface> surface;
	ref_ptr<Module> module;
public:
	RestrictToRegion(Module* _module, Surface* _surface);
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/** @}*/

} // namespace radiopropa

#endif // CRPROPA_RESTRICTTOREGION_H
