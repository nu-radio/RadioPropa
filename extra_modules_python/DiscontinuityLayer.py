import radiopropa
import numpy as np
#from ObserverPlane import ObserverPlane 


class DiscontinuityLayer(radiopropa.Module):
    """
    A layer a radiopropa surface.
    Part of a candidate crossing the layer is reflected. 
    Future propagation in scalar field will take care of directivity change automatically???
    """
    def __init__(self, surface, transmisionCoefficient):
        radiopropa.Module.__init__(self)
        self.__surface = surface
        self.__transmisionCoefficient = transmisionCoefficient

    def process(self, candidate):

        currentDistance = self.__surface.distance(candidate.current.getPosition())
        previousDistance = self.__surface.distance(candidate.previous.getPosition())

        if np.sign(currentDistance) == np.sign(previousDistance):
            candidate.limitNextStep(abs(currentDistance))
        else:
            E = candidate.current.getAmplitude()

            # The secondary propagates further, while the candidate is
            # reflected: legacy from CRPropa interface as secondaries have same
            # direction as parents.
            candidate.addSecondary(0, self.__transmisionCoefficient * E, 1)
            candidate.current.setAmplitude((1-self.__transmisionCoefficient) * E)

            normal = self.__surface.normal(candidate.current.getPosition())
            v = candidate.current.getDirection()
            cos_theta = v.dot(normal)
            u = normal * (cos_theta)
            new_direction = v - u*2 #new direction due to reflection of surface
            candidate.current.setDirection(new_direction)
            
            if cos_theta < 0.: candidate.appendReflectionAngle(np.pi - np.arccos(cos_theta))
            else: candidate.appendReflectionAngle(np.arccos(cos_theta))

            # update position slightly to move on correct side of plane
            x = candidate.current.getPosition()
            candidate.current.setPosition(x + new_direction * candidate.getCurrentStep())
            # update position (this is a hack to avoid double scatter)
            # candidate.previous.setPosition(candidate.current.getPosition())

class TransmissiveLayer(radiopropa.Module):
    """
    A layer a radiopropa surface.
    Part of a candidate crossing the layer is reflected. 
    Future propagation in scalar field will take care of directivity change automatically???
    """
    def __init__(self, surface, transmisionCoefficient):
        radiopropa.Module.__init__(self)
        self.__surface = surface
        self.__transmisionCoefficient = transmisionCoefficient

    def process(self, candidate):

        currentDistance = self.__surface.distance(candidate.current.getPosition())
        previousDistance = self.__surface.distance(candidate.previous.getPosition())

        if np.sign(currentDistance) == np.sign(previousDistance):
            candidate.limitNextStep(abs(currentDistance))
        else:
            E = candidate.current.getAmplitude()
            candidate.current.setAmplitude(self.__transmisionCoefficient * E)

class ReflectiveLayer(radiopropa.Module):
    """
    A layer a radiopropa surface.
    Part of a candidate crossing the layer is reflected. 
    Future propagation in scalar field will take care of directivity change automatically???
    """
    def __init__(self, surface, reflectionCoefficient):
        radiopropa.Module.__init__(self)
        self.__surface = surface
        self.__reflectionCoefficient = reflectionCoefficient

    def process(self, candidate):

        currentDistance = self.__surface.distance(candidate.current.getPosition())
        previousDistance = self.__surface.distance(candidate.previous.getPosition())

        if np.sign(currentDistance) == np.sign(previousDistance):
            candidate.limitNextStep(abs(currentDistance))
        else:
            E = candidate.current.getAmplitude()

            candidate.current.setAmplitude(self.__reflectionCoefficient * E)

            normal = self.__surface.normal(candidate.current.getPosition())
            v = candidate.current.getDirection()
            u = normal * (v.dot(normal))
            new_direction = v - u*2 #new direction due to reflection of surface
            candidate.current.setDirection(new_direction)

            # update position slightly to move on correct side of plane
            x = candidate.current.getPosition()
            candidate.current.setPosition(x + new_direction * candidate.getCurrentStep())
            # update position (this is a hack to avoid double scatter)
            # candidate.previous.setPosition(candidate.current.getPosition())