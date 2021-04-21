#include "radiopropa/Candidate.h"
#include "radiopropa/Geometry.h"
#include "radiopropa/Units.h"
#include "radiopropa/module/Discontinuity.h"

#include "gtest/gtest.h"

namespace radiopropa {

TEST(testDiscontinuity, limitNextStep) {

	Discontinuity discontinuity(new Plane(Vector3d(0, 0, 0), Vector3d(0, 0, 1)), 1., 1.5);

	ParticleState p;
	p.setPosition(Vector3d(0, 0, -1));
	p.setDirection(Vector3d(0, 0, 1));

	Candidate c;
	c.current = p;
	c.previous = p;
	c.setNextStep(10);

	discontinuity.process(&c);

	EXPECT_EQ(1, c.getNextStep());
}

TEST(testDiscontinuity, testRefractionPerpendicularZ) {
	double n1 = 1.;
	double n2 = 1.5;

	Discontinuity discontinuity(new Plane(Vector3d(0, 0, 0), Vector3d(0, 0, 1)), n1, n2);

	ParticleState pp, cp;
	Candidate c;
	pp.setPosition(Vector3d(0, 0, -1E-8));
	cp.setPosition(Vector3d(0, 0, 1E-8));

	// perpendicular case, no refraction
	pp.setDirection(Vector3d(0, 0, 1));
	cp.setDirection(Vector3d(0, 0, 1));

	c.current = cp;
	c.previous = pp;
	discontinuity.process(&c);

	EXPECT_EQ(c.current.getDirection(), Vector3d(0,0,-1));
	EXPECT_EQ(c.secondaries.size(), 1);
	EXPECT_EQ((c.secondaries[0])->current.getDirection(), Vector3d(0,0,1));
}

TEST(testDiscontinuity, testRefractionPerpendicularX) {
	double n1 = 1.;
	double n2 = 1.5;

	Discontinuity discontinuity(new Plane(Vector3d(0, 0, 0), Vector3d(1, 0, 0)), n1, n2);

	ParticleState pp, cp;
	Candidate c;
	pp.setPosition(Vector3d(-1E-8, 0, 0));
	cp.setPosition(Vector3d(1E-8, 0, 0));

	// perpendicular case, no refraction
	pp.setDirection(Vector3d(1, 0, 0));
	cp.setDirection(Vector3d(1, 0, 0));

	c.current = cp;
	c.previous = pp;
	discontinuity.process(&c);

	EXPECT_EQ(c.current.getDirection(), Vector3d(-1,0,0));
	EXPECT_EQ(c.secondaries.size(), 1);
	EXPECT_EQ((c.secondaries[0])->current.getDirection(), Vector3d(1,0,0));
}

TEST(testDiscontinuity, testRefractionAngleZ) {
	double n1 = 1.;
	double n2 = 1.5;

	Discontinuity discontinuity(new Plane(Vector3d(0, 0, 0), Vector3d(0, 0, 1)), n1, n2);

	ParticleState pp, cp;
	Candidate c;
	pp.setPosition(Vector3d(0, 0, -1E-8));
	cp.setPosition(Vector3d(0, 0, 1E-8));

	Vector3d v;
	v.setRThetaPhi(1., 10 * deg, 0);
	// perpendicular case, no refraction
	pp.setDirection(v);
	cp.setDirection(v);

	c.current = cp;
	c.previous = pp;
	discontinuity.process(&c);

	EXPECT_NEAR(c.current.getDirection().getAngleTo(Vector3d(0,0,1)), 170*deg, 1E-8);
	EXPECT_EQ(c.secondaries.size(), 1);
	EXPECT_NEAR((c.secondaries[0])->current.getDirection().getTheta(), asin(sin(10. * deg) / 1.5), 1e-8);
}




int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace radiopropa
