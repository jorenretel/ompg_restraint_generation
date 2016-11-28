
from collections import defaultdict
from itertools import combinations, product
from ccpnmr.analysis.core.MoleculeBasic import areResonancesBound
from ccpnmr.analysis.core.AssignmentBasic import findMatchingPeakDimShifts
from ccpnmr.analysis.core.ConstraintBasic import getPeakDimTolerance, isResidueInRange
from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance
from labeling_simple import getExperimentResonanceSetFractions




def test(argServer):
    project = argServer.getProject()
    nmrProject = project.findFirstNmrProject()
    nmrConstraintStore = project.newNmrConstraintStore(nmrProject=nmrProject)

    #HHN
    expName = '248_07HHN'
    peakListSerial = 2
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)

    dataDim1, dataDim2, dataDim3 = spectrum.sortedDataDims()
    tolerances = [(dataDim1, 0.07, 0.07, 1.0), (dataDim2, 0.1, 0.1, 1.0), (dataDim3, 0.4, 0.4, 1.0)]
    chemShiftRanges = [(dataDim1, '1H', 0.0,13.0), (dataDim2, '1H', 0.0, 13.0),(dataDim3, '15N', 100.0, 140.0)]

    for peak in peakList.peaks:
        print peak
        if peak_is_diagonal(peak, tolerances):
            continue
        if peak_has_poor_merit(peak, 1.0):
            continue
        if peak_is_out_of_chemical_shift_ranges(peak, chemShiftRanges):
            continue
        print 'assignment status ', peak_is_fully_assigned(peak)

        print find_peak_assignment_options(peak, tolerances, chemShiftRanges)


def find_peak_assignment_options(peak, tolerances, chemShiftRanges,
                                 aliasing=True, onlyAssignedResonances=True,
                                 onlyDimensionalAssignmentWhenPresent=False,
                                 labelling=True, minLabelFraction=0.1):

    experiment = peak.peakList.dataSource.experiment
    if labelling is True:
        labelling = experiment

    tolerances = dict([(dim,(mintol, maxtol, multi)) for dim, mintol, maxtol, multi in tolerances])

    chemShiftRangeDict = defaultdict(list)
    for dim, iso, minshift, maxshift in chemShiftRanges:
        chemShiftRangeDict[dim].append([minshift, maxshift])
    chemShiftRanges = chemShiftRangeDict

    dimensional_options = []
    for peakDim in peak.sortedPeakDims():
        tolerance = tolerances[peakDim.dataDim]
        chemShiftRange = chemShiftRanges[peakDim.dataDim]
        resonances = find_dimensional_assignment_options(peakDim, tolerance,
                                                         chemShiftRange,
                                                         aliasing= aliasing,
                                                         onlyAssignedResonances=onlyAssignedResonances,
                                                         onlyDimensionalAssignmentWhenPresent=onlyDimensionalAssignmentWhenPresent)
        #if not resonances:
        #    return []
        dimensional_options.append(resonances)

    peak_assignment_options = []
    for resonances in product(*dimensional_options):
        # 1. experiment
        if not resonances_fit_experiment(resonances, experiment):
            continue
        # 2. labelling
        if labelling:
            colabelling = getExperimentResonanceSetFractions(labelling, resonances)
            if colabelling < minLabelFraction:
                continue
        peak_assignment_options.append(resonances)

    return peak_assignment_options


def find_dimensional_assignment_options(peakDim, tolerance, chemShiftRange,
                                        aliasing=True, onlyAssignedResonances=True,
                                        onlyDimensionalAssignmentWhenPresent=False):

    tolerance = getPeakDimTolerance(peakDim, *tolerance)
    if onlyDimensionalAssignmentWhenPresent:
        resonances = resonances_assigned_to_peak_dimension(peakDim, onlyAssignedResonnces=onlyAssignedResonances)
        if resonances:
            return resonances

    resonances = shift_match_resonances(peakDim, tolerance, chemShiftRange,
                                        aliasing=aliasing,
                                        onlyAssignedResonances=onlyAssignedResonances)
    return resonances


def resonances_assigned_to_peak_dimension(peakDim, onlyAssignedResonances=True):

    resonance_assignments = []
    for peakDimContrib in peakDim.peakDimContribs:
        resonance = peakDimContrib.resonance
        resonanceSet = resonance.resonanceSet
        if resonanceSet or not onlyAssignedShifts:
            resonance_assignments.append(resonance)
    return resonance_assignments


def shift_match_resonances(peakDim, tolerance, chemShiftRange,
                           aliasing, onlyAssignedResonances=True):

    shifts = findMatchingPeakDimShifts(peakDim, chemShiftRange,
                                       tolerance=tolerance,
                                       aliasing=aliasing,
                                       findAssigned=onlyAssignedResonances)
    resonances = [shift.resonance for shift in shifts if shift.resonance]
    return resonances


def peak_is_diagonal(peak, tolerances):
    '''Only recognizes peaks that are completely diagonal
       in the sense that all dimensions measuring the same
       isotope have the same chemical shift. There are
       spectra imaginable producing peaks that could
       be encoding a 'transfer to self' without all
       dimensions of the same isotope having the same
       shift.

    '''

    dimensions_by_isotope = group_peak_dimensions_by_isotope(peak)
    tolerance_dict = {}
    for dim, mintol, maxtol, multiplyer in tolerances:
        tolerance_dict[dim] = maxtol

    for dimensions in dimensions_by_isotope.values():
        for dim1, dim2 in combinations(dimensions, 2):
            delta = abs(dim1.value - dim2.value)
            tol = max(tolerance_dict[dim1.dataDim], tolerance_dict[dim2.dataDim])
            if delta > tol:
                return False
    return True


def group_peak_dimensions_by_isotope(peak):

    dimensions_by_isotope = defaultdict(list)
    for dim in peak.sortedPeakDims():
        isotope = dim.dataDim.expDim.findFirstExpDimRef().isotopeCodes[0]
        dimensions_by_isotope[isotope].append(dim)
    return dimensions_by_isotope


def peak_has_poor_merit(peak, minMerit):
    return peak.figOfMerit < minMerit


def peak_is_out_of_chemical_shift_ranges(peak, chemShiftRanges):

    range_dict = {}
    for dim, isotope, minshift, maxshift in chemShiftRanges:
        range_dict[dim] = (minshift, maxshift)
    for dim in peak.sortedPeakDims():
        minshift, maxshift = range_dict[dim.dataDim]
        if not minshift < dim.value < maxshift:
            True
    return False


def peak_is_fully_assigned(peak):

    for dim in peak.sortedPeakDims():
        if not dim.peakDimContribs:
            return False
    return True


def resonances_fit_experiment(resonances, experiment):
    '''Check whether a combination of resonances fit the
       experiment type. Resonances are expected to be given
       ordered by peakDim.dim numbers (which directly corresponds
       to dataDim.dim).

    '''

    dataDims = experiment.findFirstDataSource().sortedDataDims()
    transferDict = {}
    for expTransfer in experiment.expTransfers:
        key = tuple([expDimRef.expDim for expDimRef in expTransfer.expDimRefs])
        transferDict[key] = expTransfer

    for dataDim1, dataDim2 in combinations(dataDims, 2):
        resonance1 = resonances[dataDim1.dim - 1]
        resonance2 = resonances[dataDim2.dim - 1]
        expTransfer = transferDict.get((dataDim1.expDim, dataDim2.expDim)) or transferDict.get((dataDim2.expDim, dataDim1.expDim))
        if not expTransfer:
            continue
        if not transfer_is_possible([resonance1, resonance2], expTransfer):
            return False
    return True


def transfer_is_possible(resonances, expTransfer):
    '''Check whether a transfer is possible between two
       resonances. Not very complete, should be made
       more specific.

    '''

    transferType = expTransfer.transferType
    if transferType == 'onebond' and not areResonancesBound(*resonances):
        return False
    #elif transferType == 'Jcoupling':
    return True
