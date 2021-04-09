#ifndef CRPROPA_ABSTRACT_OUTPUT_H
#define CRPROPA_ABSTRACT_OUTPUT_H

#include "radiopropa/Module.h"
#include "radiopropa/Variant.h"

#include <bitset>
#include <vector>
#include <string>

namespace radiopropa {

/**
 @class Output
 @brief Configurable output base class.
 */
class Output: public Module {
protected:
	double lengthScale, frequencyScale;
	std::bitset<64> fields;

	struct Property
	{
		std::string name;
		std::string comment;
		Variant defaultValue;
	};
	std::vector<Property> properties;

	bool oneDimensional;
	mutable size_t count;

	void modify();

public:
	enum OutputColumn {
		TrajectoryLengthColumn,
		PropagationTimeColumn,
		ColumnDensityColumn,
		CurrentIdColumn,
		CurrentAmplitudeColumn,
		CurrentFrequencyColumn,
		CurrentPositionColumn,
		CurrentDirectionColumn,
		SourceIdColumn,
		SourceFrequencyColumn,
		SourceAmplitudeColumn,
		SourcePositionColumn,
		SourceDirectionColumn,
		CreatedIdColumn,
		CreatedFrequencyColumn,
		CreatedPositionColumn,
		CreatedDirectionColumn,
		SerialNumberColumn,
		WeightColumn
	};
	enum OutputType {
		Trajectory1D,
		Trajectory3D,
		Event1D,
 		Event3D,
		Everything
	};

	std::string OutputTypeName(OutputType outputtype);
	const std::string outputName;
	
	Output();
	Output(OutputType outputtype);

	void setFrequencyScale(double scale);
	void setLengthScale(double scale);

	void setOutputType(OutputType outputtype);
	void set(OutputColumn field, bool value);
	/// Add a property to output. Default value is required to assign a type in
	/// the output
	void enableProperty(const std::string &property, const Variant& defaultValue, const std::string &comment = "");
	void enable(OutputColumn field);
	void disable(OutputColumn field);
	void enableAll();
	void disableAll();
	void set1D(bool value);
	size_t size() const;

	void process(Candidate *) const;
};

} // namespace radiopropa

#endif // CRPROPA_ABSTRACT_OUTPUT_H
