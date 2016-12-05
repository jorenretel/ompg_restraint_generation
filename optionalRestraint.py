
from ccpnmr.analysis.core.ConstraintBasic import getNoeDistance, getIntensityDistanceTable, getDistancesFromIntensity
from ccpnmr.analysis.core.ExperimentBasic import getOnebondDataDims


class OptionalRestraint(object):
    """docstring for OptionalRestraint"""
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

        contribution = OptionalRestraintContribution(assignment,
                                                     restraint_resonances,
                                                     labelling)
        self.contributions.append(contribution)

    def apply_round_tolerance(self, tolerances):

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
    """docstring for OptionalRestraintContribution"""
    def __init__(self, assignment, restraint_resonances, labelling):
        super(OptionalRestraintContribution, self).__init__()
        self.assignment = assignment
        self.restraint_resonances = restraint_resonances
        self.labelling = labelling
        self.symmetry = 0.0
