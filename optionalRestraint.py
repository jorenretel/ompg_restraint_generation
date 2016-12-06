'''This module contains classes representing restraints.
   These are not the classes from which the restraint
   objects stem that are saved in the CCPN project. These
   classes are just used to be able to be able to easily
   filter restraints and contributions to restraints
   before generating the more permanent object stored
   in the CCPN data model.

'''

#pylint: disable=invalid-name


from ccpnmr.analysis.core.ConstraintBasic import getNoeDistance, getIntensityDistanceTable, getDistancesFromIntensity
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims


class OptionalRestraint(object):
    """All information for a restraint belonging to
       one peak. Depending on the applied criteria,
       might not end up in the final restraint set
       saved in the CCPN project.

    """
    def __init__(self, peak, peak_assignment_options, restraint_options,
                 labelling=None, intensityValue=None):
        super(OptionalRestraint, self).__init__()
        self.peak = peak
        self.intensityValue = intensityValue
        self.targetValue = None
        self.upperLimit = None
        self.lowerLimit = None
        self.error = None
        self.contributions = []

        for assignment, restraint_resonances, label in zip(peak_assignment_options,
                                                           restraint_options,
                                                           labelling):

            self.add_contribution(assignment, restraint_resonances, label)




    def calculate_distance(self, intensityScale=None,
                           distanceFunction=None, params=None,
                           intensityType='volume'):
        '''Uses distanceFunction to calculate the distance bases on
           the intensity of the peak.

        '''

        if not distanceFunction:
            if params:
                distanceFunction = lambda val:getNoeDistance(val, params)
            else:
                spectrum = self.peak.peakList.dataSource
                noeDistClasses = getIntensityDistanceTable(spectrum)
                distanceFunction = lambda val:getDistancesFromIntensity(noeDistClasses, val)

        intensity = self.peak.findFirstPeakIntensity(intensityType=intensityType)
        intensityValue = abs(intensity.value)

        target, lower, upper = distanceFunction(intensityValue/float(intensityScale))

        self.intensityValue = intensityValue
        self.lowerLimit = lower
        self.upperLimit = upper
        self.targetValue = target
        self.error = abs(lower - upper)

    def add_contribution(self, assignment, restraint_resonances, labelling):
        '''Add one contribution to the restraint.'''

        contribution = OptionalRestraintContribution(assignment,
                                                     restraint_resonances,
                                                     labelling)
        self.contributions.append(contribution)

    def apply_round_tolerance(self, tolerances):
        '''Evaluates dimensions corresponding to a through-bond
           transfer together and uses the deviation between the
           expected peak position (based on the assigned chemical
           shifts) and the actual peak position to make the
           tolerance 'round' instead of square. All possible
           contributions to the restraints that deviate from
           the expected peak position within twice the deviation
           of the least deviating contribution or deviate less
           than half of the maximum deviation still within the
           tolerance square (2**0.5/2) are kept. All other
           options are removed.

        '''

        if not self.contributions:
            return

        tolerances = dict([(dim, maxtol) for dim, mintol, maxtol, multi in tolerances])
        shiftList = self.peak.peakList.dataSource.experiment.shiftList
        peakDims = self.peak.sortedPeakDims()
        spectrum = self.peak.peakList.dataSource

        for dataDims in getOnebondDataDims(spectrum):
            distances = []
            for contrib in self.contributions:
                sum_of_squares = 0
                for dataDim in dataDims:
                    dim = dataDim.dim
                    peakDim = peakDims[dim-1]
                    resonance = contrib.assignment[dim-1]
                    tolerance = tolerances[dataDim]
                    v1 = resonance.findFirstShift(parentList=shiftList).value
                    v2 = peakDim.value
                    delta_squared = ((v1-v2)/tolerance)**2
                    sum_of_squares += delta_squared
                distances.append(sum_of_squares**0.5)

            closest = min(distances)
            for distance, contrib in zip(distances, self.contributions):
                if distance > 2*closest and distance > 2**0.5/2:
                    self.contributions.remove(contrib)



class OptionalRestraintContribution(object):
    """A Contribution to an optional restraint."""
    def __init__(self, assignment, restraint_resonances, labelling):
        super(OptionalRestraintContribution, self).__init__()
        self.assignment = assignment
        self.restraint_resonances = restraint_resonances
        self.labelling = labelling
        self.symmetry = 0.0
