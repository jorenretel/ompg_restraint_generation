
from ccpnmr.analysis.core.ConstraintBasic import getNoeDistance, getIntensityDistanceTable, getDistancesFromIntensity


class OptionalRestraint(object):
	"""docstring for OptionalRestraint"""
	def __init__(self, peak, peak_assignment_options, restraint_options,
		         labelling=None, distance=None, intensity=None):
		super(OptionalRestraint, self).__init__()
		self.peak = peak
		self.peak_assignment_options = peak_assignment_options
		self.restraint_options = restraint_options
		self.labelling = labelling
		self.distance = distance
		self.intensity = intensity
		self.targetValue = None
		self.upperLimit = None
		self.lowerLimit = None
		self.error = None


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



