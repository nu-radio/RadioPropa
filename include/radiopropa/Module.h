#ifndef CRPROPA_MODULE_H
#define CRPROPA_MODULE_H

#include "radiopropa/Candidate.h"
#include "radiopropa/Referenced.h"
#include "radiopropa/Common.h"

#include <string>

namespace radiopropa {

class Candidate;

/**
 @class Module
 @brief Abstract base class for modules
 */
class Module: public Referenced {
	std::string description;
public:
	Module();
	virtual ~Module() {
	}
	virtual std::string getDescription() const;
	void setDescription(const std::string &description);
	virtual void process(Candidate *candidate) const = 0;
	inline void process(ref_ptr<Candidate> candidate) const {
		process(candidate.get());
	}
};


/**
 @class AbstractCondition
 @brief Abstract Module providing common features for conditional modules.
 */
class AbstractCondition: public Module {
protected:
	ref_ptr<Module> rejectAction, acceptAction;
	bool makeRejectedInactive, makeAcceptedInactive;
	std::string rejectFlagKey, rejectFlagValue;
	std::string acceptFlagKey, acceptFlagValue;

	void reject(Candidate *candidate) const;
	inline void reject(ref_ptr<Candidate> candidate) const {
		reject(candidate.get());
	}

	void accept(Candidate *candidate) const;
	inline void accept(ref_ptr<Candidate> candidate) const {
		accept(candidate.get());
	}

public:
	AbstractCondition();
	void onReject(Module *rejectAction);
	void onAccept(Module *acceptAction);
	void setMakeRejectedInactive(bool makeInactive);
	void setMakeAcceptedInactive(bool makeInactive);
	void setRejectFlag(std::string key, std::string value);
	void setAcceptFlag(std::string key, std::string value);
};

} // namespace radiopropa

#endif /* CRPROPA_MODULE_H */
