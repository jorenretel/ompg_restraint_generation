
import restraintGeneration
reload(restraintGeneration)
from restraintGeneration import make_optional_restraint_set, record_constraint_list, ambiguity_info
import symmetryFilter
reload(symmetryFilter)
from symmetryFilter import symmetry_filter
import gitStamp
reload(gitStamp)
from gitStamp import get_git_hash


def create_ompg_restraints(argServer):
    '''Creates the distance restraints used for the structure
    determination of Outer Membrane Protein G (OmpG).

     '''

    # I stamp the created restraint sets with
    # the current git hash, so I know which version
    # of this script and parameters were used.
    stamp = 'git: ' + get_git_hash()

    project = argServer.getProject()
    nmrProject = project.findFirstNmrProject()
    #molSystem = project.findFirstMolSystem(name='MS1')
    #chain = molSystem.findFirstChain()

    # All parameters for shift-matching hNhhNH
    expName = '250_HNN_30_10_sp035'
    peakListSerial = 2
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)
    dataDim1, dataDim2, dataDim3 = spectrum.sortedDataDims()
    tolerances = [(dataDim1, 0.07, 0.07, 1.0), (dataDim2, 0.4, 0.4, 1.0), (dataDim3, 0.4, 0.4, 1.0)]
    chemShiftRanges = [(dataDim1, '1H', 0.0,13.0), (dataDim2, '15N', 100.0, 140.0),(dataDim3, '15N', 100.0, 140.0)]
    distanceFunction = DistanceFunctionProtonDetected([(2500000, 3.1, 1.0, 3.5), (0.0, 4.4, 1.0, 5.5)])
    optionalRestraintsNNH = make_optional_restraint_set(peakList, tolerances,
    	                                                chemShiftRanges,
    	                                                distanceFunction=distanceFunction,
    	                                                intensityType='height')
    toleranceString = str(tolerances[0][1]) + '_' + str(tolerances[1][1]) + '_' + str(tolerances[2][1])
    NNH_name = 'NNH_peaklist{}_shiftmatch_{}_symmetry_filetered'.format(peakListSerial, toleranceString)

    # All parameters for shift-matching hNHH
    expName = '248_07HHN'
    peakListSerial = 2
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)
    dataDim1, dataDim2, dataDim3 = spectrum.sortedDataDims()
    tolerances = [(dataDim1, 0.07, 0.07, 1.0), (dataDim2, 0.1, 0.1, 1.0), (dataDim3, 0.4, 0.4, 1.0)]
    chemShiftRanges = [(dataDim1, '1H', 0.0,13.0), (dataDim2, '1H', 0.0, 13.0),(dataDim3, '15N', 100.0, 140.0)]
    distanceFunction = DistanceFunctionProtonDetected([(1500000, 3.1, 1.0, 3.5), (0.0, 4.4, 1.0, 5.5)])
    optionalRestraintsHHN = make_optional_restraint_set(peakList, tolerances,
    	                                                chemShiftRanges,
    	                                                distanceFunction=distanceFunction,
    	                                                intensityType='height')
    toleranceString = str(tolerances[0][1]) + '_' + str(tolerances[1][1]) + '_' + str(tolerances[2][1])
    HHN_name = 'HHN_peaklist{}_shiftmatch_{}_symmetry_filetered'.format(peakListSerial, toleranceString)

    # Use redundancy from four peaks in hNhhNH and hNHH to select most
    # likely restraint items:
    print 'NNH before:'
    ambiguity_info(optionalRestraintsNNH)
    print 'HHNH before:'
    ambiguity_info(optionalRestraintsHHN)

    symmetry_filter([optionalRestraintsHHN, optionalRestraintsNNH], cutoff_fraction=1.0)

    print 'NNH after:'
    ambiguity_info(optionalRestraintsNNH)
    print 'NNH after:'
    ambiguity_info(optionalRestraintsHHN)

    nmrConstraintStore = project.newNmrConstraintStore(nmrProject=nmrProject)

    # Make the optional proton-proton restraints into real ones and save them:
    HHN = record_constraint_list(optionalRestraintsHHN, nmrConstraintStore)
    NNH = record_constraint_list(optionalRestraintsNNH, nmrConstraintStore)
    HHN.name = HHN_name
    NNH.name = NNH_name
    NNH.details = stamp
    HHN.details = stamp

    # Long mixing time experiments:
    shift_match_2D_CC(nmrProject, '029_PDSD400_13C', 4, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '063_PDSD400_2C', 5, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '045_PDSD400_13TEMPQANDSG', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '113_PDSD400_2TEMPQANDSG', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '110_PDSD400_2SLYGWAFVL', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '160_PDSD500_GAFY', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4, details=stamp)

    # Medium mixing time spectra:
    shift_match_2D_CC(nmrProject, '076_PDSD200_2C', 2, nmrConstraintStore, (3.0, 1.5, 5.5), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '004_PDSD200_13C', 2, nmrConstraintStore, (3.0, 1.5, 5.5), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '112_PDSD150_2TEMPQANDSG', 2, nmrConstraintStore, (3.0, 1.5, 5.5), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '046_PDSD150_13TEMPQANDSG', 2, nmrConstraintStore, (3.0, 1.5, 5.5), 0.4, details=stamp)
    shift_match_2D_CC(nmrProject, '106_PDSD150_2SLYGWAFVL', 2, nmrConstraintStore, (3.0, 1.5, 5.5), 0.4, details=stamp)



def shift_match_2D_CC(nmrProject, expName, peakListSerial,
                      nmrConstraintStore, bucket, tolerance, details=None):
    '''Just a function
    '''

    # One length distance bucket
    dist_function_13c = lambda x: bucket
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)
    dataDim1, dataDim2 = spectrum.sortedDataDims()
    tol = tolerance
    tolerances = [(dataDim1, tol, tol, 1.0), (dataDim2, tol, tol, 1.0)]
    chemShiftRanges = [(dataDim1, '13C', 0.0, 200.0), (dataDim2, '13C', 0.0, 200.0)]
    constraints = make_optional_restraint_set(peakList, tolerances,
                                              chemShiftRanges,
                                              labelling=True,
                                              distanceFunction=dist_function_13c)

    print expName
    ambiguity_info(constraints)
    constraints =record_constraint_list(constraints, nmrConstraintStore)
    toleranceString = str(tolerances[0][1]) + '_' + str(tolerances[1][1])
    constraints.name = '{}_{}_shiftmatch_{}'.format(expName, peakListSerial, toleranceString)
    constraints.details = details

    return constraints


class DistanceFunctionProtonDetected(object):

    def __init__(self, distance_bins):

        distance_bins.sort(reverse=True)
        self.intensity_distance_bins = distance_bins

    def __call__(self, intensity):

        for threshold, dist, minDist, maxDist in self.intensity_distance_bins:
            if intensity >= threshold:
                return (dist, minDist, maxDist)