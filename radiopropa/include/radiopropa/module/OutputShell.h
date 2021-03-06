#ifndef RADIOPROPA_OUTPUTSHELL_H
#define RADIOPROPA_OUTPUTSHELL_H

#include "radiopropa/Module.h"
#include "radiopropa/AssocVector.h"
#include "radiopropa/Variant.h"

namespace radiopropa {

/**
 @class ShellOutput
 @brief Show the trajectory in the shell.
 */
class ShellOutput: public Module {
public:
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellOutput1D
 @brief Show the trajectory in the shell.
 */
class ShellOutput1D: public Module {
public:
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellPropertyOutput
 @brief Show the candidate properties in the shell.
 */
class ShellPropertyOutput: public Module {
public:
	typedef Loki::AssocVector<std::string, Variant> PropertyMap;
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace cprpropa

#endif // RADIOPROPA_OUTPUTSHELL_H
