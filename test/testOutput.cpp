/** Unit tests for Output modules of RadioPropa
    Output
    TextOutput
    ParticleCollector
 */

#include "RadioPropa.h"

#include <string>
#include "gtest/gtest.h"
#include <iostream>

// compare two arrays (intead of using Google Mock)
// https://stackoverflow.com/a/10062016/6819103
template<typename T, size_t size>
::testing::AssertionResult ArraysMatch(const T (&expected)[size], 
		const T (&actual)[size]){
	for (size_t i(0); i < size; ++i){
		if (expected[i] != actual[i]){
			return ::testing::AssertionFailure() << "array[" << i
			<< "] (" << actual[i] << ") != expected[" << i
			<< "] (" << expected[i] << ")";
		}
	}

	return ::testing::AssertionSuccess();
}

namespace radiopropa {

//-- Output

TEST(Output, size) {
	Candidate c;
	Output output;
	for (int it=0; it<5; ++it, output.process(&c));

	EXPECT_EQ(output.size(), 5);
}

//-- TextOutput


TEST(TextOutput, printProperty) {
	Candidate c;
	TextOutput output(Output::Event1D);
	output.disableAll();
	output.enableProperty("foo", 2.0, "Bar");

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	// name in first line of header
	EXPECT_EQ(captured.substr(0, captured.find("\n")),
	          "#\tfoo");
}

TEST(TextOutput, printHeader_Version) {
	Candidate c;
	TextOutput output(Output::Event1D);

	::testing::internal::CaptureStdout();
	output.process(&c);
	std::string captured = testing::internal::GetCapturedStdout();

	// length of the prefix is 19 chars
	size_t version_pos = captured.find("# RadioPropa version: ") + 22;

	EXPECT_EQ(
		captured.substr(
			version_pos,
		      	captured.find("\n", version_pos) - version_pos
			),
	         g_GIT_DESC);
}

//-- ParticleCollector

TEST(ParticleCollector, size) {
	ref_ptr<Candidate> c = new Candidate();
	ParticleCollector output;

	for (int it=0; it<5; ++it, output.process(c));

	EXPECT_EQ(output.size(), 5);
}

TEST(ParticleCollector, fetchItem) {
	ref_ptr<Candidate> c = new Candidate(1, 1*EeV);
	ParticleCollector output;

	output.process(c);

	EXPECT_EQ(output[0], c);
}

TEST(ParticleCollector, reprocess) {
	ref_ptr<Candidate> c = new Candidate(1, 1*EeV);
	ParticleCollector collector;
	ParticleCollector output;

	collector.process(c);
	collector.reprocess(&output);

	EXPECT_EQ(output[0], c);
}

TEST(ParticleCollector, dumpload) {
	ref_ptr<Candidate> c = new Candidate(1, 1.234*EeV);
	c->current.setPosition(Vector3d(1,2,3));
	c->current.setDirection(Vector3d(-1,-1,-1));
	c->setTrajectoryLength(1*Mpc);
	c->current.setAmplitude(Vector3d(2,0,0));

	ParticleCollector input;
	ParticleCollector output;

	for(int i=0; i<=10; ++i){
		input.process(c);
	}

	// Well, it would be nicer if we don't need to receate any file
	input.dump("ParticleCollector_DumpTest.txt");
	output.load("ParticleCollector_DumpTest.txt");

	EXPECT_EQ(input.size(), output.size());
	EXPECT_EQ(output[0]->current.getFrequency(), c->current.getFrequency());
	EXPECT_EQ(output[1]->getTrajectoryLength(), c->getTrajectoryLength());
	EXPECT_EQ(output[3]->current.getAmplitude().x, c->current.getAmplitude().x);
	EXPECT_EQ(output[3]->current.getAmplitude().y, c->current.getAmplitude().y);
	EXPECT_EQ(output[3]->current.getAmplitude().z, c->current.getAmplitude().z);
}

// Just test if the trajectory is on a line for rectilinear propagation
TEST(ParticleCollector, getTrajectory) {
	int pos_x[10];
	int pos_x_expected[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};


	ParticleState p;
	p.setPosition(Vector3d(10, 0, 0));
	p.setDirection(Vector3d(-1, 0, 0));
	ref_ptr<Candidate> c = new Candidate(p);

	ref_ptr<ParticleCollector> output = new ParticleCollector();
	ref_ptr<ParticleCollector> trajectory = new ParticleCollector();
	trajectory->setClone(true);

	ref_ptr<ModuleList> sim = new ModuleList();
	sim->add(new SimplePropagation(1, 1));

	ref_ptr<Observer> obs = new Observer();
        obs->add(new ObserverPoint());
	obs->onDetection(output);
	sim->add(obs);

	sim->run(c);

	output->getTrajectory(sim, 0, trajectory);

	Vector3d pos; int i;

	for (ParticleCollector::iterator itr = trajectory->begin(); itr != trajectory->end(); ++itr){
		pos = (*(itr->get())).current.getPosition();
		pos_x[i] = pos.getX();
		++i;
	}

	EXPECT_TRUE(ArraysMatch(pos_x_expected, pos_x));
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace radiopropa
