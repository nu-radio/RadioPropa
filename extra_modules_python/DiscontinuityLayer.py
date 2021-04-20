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
        self.__times_reflectedoff ={}

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
            candidate.addSecondary(0,E*self.__transmisionCoefficient, 1)
            candidate.current.setAmplitude(E*np.sqrt(1-self.__transmisionCoefficient**2))

            normal = self.__surface.normal(candidate.current.getPosition())
            v = candidate.current.getDirection()
            cos_theta = v.dot(normal)
            u = normal * (cos_theta)
            new_direction = v - u*2 #new direction due to reflection of surface
            candidate.current.setDirection(new_direction)
            
            #if cos_theta < 0.: candidate.appendReflectionAngle(np.pi - np.arccos(cos_theta))
            #else: candidate.appendReflectionAngle(np.arccos(cos_theta))

            # update position slightly to move on correct side of plane
            x = candidate.current.getPosition()
            candidate.current.setPosition(x + new_direction * candidate.getCurrentStep())
            # update position (this is a hack to avoid double scatter)
            # candidate.previous.setPosition(candidate.current.getPosition())
            if candidate not in self.__times_reflectedoff.keys(): 
                self.__times_reflectedoff[candidate] = 1
            else:
                self.__times_reflectedoff[candidate] += 1

    def get_times_reflectedoff(self, candidate):
        return self.__times_reflectedoff[candidate]

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
            candidate.current.setAmplitude(E*self.__transmisionCoefficient)

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
        self.__times_reflectedoff ={}

    def process(self, candidate):

        currentDistance = self.__surface.distance(candidate.current.getPosition())
        previousDistance = self.__surface.distance(candidate.previous.getPosition())

        if np.sign(currentDistance) == np.sign(previousDistance):
            candidate.limitNextStep(abs(currentDistance))
        else:
            E = candidate.current.getAmplitude()

            candidate.current.setAmplitude(E*self.__reflectionCoefficient)

            normal = self.__surface.normal(candidate.current.getPosition())
            v = candidate.current.getDirection()
            cos_theta = v.dot(normal)
            u = normal * (cos_theta)
            new_direction = v - u*2 #new direction due to reflection of surface
            candidate.current.setDirection(new_direction)

            #if cos_theta < 0.: candidate.appendReflectionAngle(np.pi - np.arccos(cos_theta))
            #else: candidate.appendReflectionAngle(np.arccos(cos_theta))

            # update position slightly to move on correct side of plane
            x = candidate.current.getPosition()
            candidate.current.setPosition(x + new_direction * candidate.getCurrentStep())
            # update position (this is a hack to avoid double scatter)
            # candidate.previous.setPosition(candidate.current.getPosition())
            if candidate not in self.__times_reflectedoff.keys(): 
                self.__times_reflectedoff[candidate] = 1
            else:
                self.__times_reflectedoff[candidate] += 1

    def get_times_reflectedoff(self, candidate):
        return self.__times_reflectedoff[candidate]