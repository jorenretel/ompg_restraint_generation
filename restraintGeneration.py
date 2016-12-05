
#pylint: disable=line-too-long, invalid-name

from collections import defaultdict
from itertools import combinations, product
from ccpnmr.analysis.core.MoleculeBasic import areResonancesBound, getBoundAtoms
from ccpnmr.analysis.core.AssignmentBasic import findMatchingPeakDimShifts, getBoundResonances, assignAtomsToRes
from ccpnmr.analysis.core.ConstraintBasic import getPeakDimTolerance, getMeanPeakIntensity, getFixedResonance #isResidueInRange
from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpaceDataDims, getIndirectDataDims, getIndirectThroughSpaceIsotopes
from labeling_simple import getExperimentResonanceSetFractions
import optionalRestraint
reload(optionalRestraint)
from optionalRestraint import OptionalRestraint


def record_constraint_list(optional_restraints, constraintSet):
    '''Take a list of OptionalRestraints and make
       more permanent restraints based on them that
       get stored in a constraintSet in CCPN.
       Done like this because the CCPN constraints
       do not have information about the direction
       of magnetization transfer.

    '''

    peak = optional_restraints[0].peak
    peakList = peak.peakList
    spectrum = peakList.dataSource
    experiment = spectrum.experiment

    distConstraintList = constraintSet.newDistanceConstraintList()
    distConstraintList.addExperimentSerial(experiment.serial)
    newConstraint = distConstraintList.newDistanceConstraint

    for optional_restraint in optional_restraints:

        constraint = newConstraint(weight=1.0,
                                   origData=optional_restraint.intensityValue,
                                   targetValue=optional_restraint.targetValue,
                                   upperLimit=optional_restraint.upperLimit,
                                   lowerLimit=optional_restraint.lowerLimit,
                                   error=optional_restraint.error)


        peak = optional_restraint.peak
        peakList = peak.peakList
        spectrum = peakList.dataSource
        experiment = spectrum.experiment

        constraint.newConstraintPeakContrib(experimentSerial=experiment.serial,
                                            dataSourceSerial=spectrum.serial,
                                            peakListSerial=peakList.serial,
                                            peakSerial=peak.serial)

        for contrib in optional_restraint.contributions:
            res0, res1 = contrib.restraint_resonances
            fixedResonance0 = getFixedResonance(constraintSet, res0)
            fixedResonance1 = getFixedResonance(constraintSet, res1)
            constraint.newDistanceConstraintItem(resonances=[fixedResonance0, fixedResonance1])

    return distConstraintList



def make_optional_restraint_set(peakList, tolerances, chemShiftRanges,
                                aliasing=True, onlyAssignedResonances=True,
                                onlyDimensionalAssignmentWhenPresent=False,
                                labelling=None, minLabelFraction=0.1,
                                structure=None, maxDist=None, scale=False,
                                intensityType='volume', ignoreDiagonals=True,
                                ignoreAssignedPeaks=True, round_tolerance=False,
                                distanceFunction=None, params=None, minMerit=0.0):

    optional_restraints = []
    for peak in list(peakList.peaks):

        if ignoreDiagonals and peak_is_diagonal(peak, tolerances):
            continue
        if peak_has_poor_merit(peak, minMerit):
            continue
        if peak_is_out_of_chemical_shift_ranges(peak, chemShiftRanges):
            continue
        if ignoreAssignedPeaks and peak_is_fully_assigned(peak):
            continue

        optional_restraint = create_optional_restraint(peak, tolerances, chemShiftRanges,
                                                       aliasing=aliasing, onlyAssignedResonances=onlyAssignedResonances,
                                                       onlyDimensionalAssignmentWhenPresent=onlyDimensionalAssignmentWhenPresent,
                                                       labelling=labelling, minLabelFraction=minLabelFraction,
                                                       structure=structure, maxDist=maxDist, round_tolerance=round_tolerance)
        if not optional_restraint:
            continue
        optional_restraints.append(optional_restraint)

    intensityScale = 1.0
    if scale:
        intensityScale = getMeanPeakIntensity([restraint.peak for restraint in optional_restraints], intensityType=intensityType)

    for optional_restraint in optional_restraints:
        optional_restraint.calculate_distance(intensityScale=intensityScale,
                                              distanceFunction=distanceFunction,
                                              params=params,
                                              intensityType=intensityType)

    return optional_restraints





def create_optional_restraint(peak, tolerances, chemShiftRanges,
                              aliasing=True, onlyAssignedResonances=True,
                              onlyDimensionalAssignmentWhenPresent=False,
                              labelling=None, minLabelFraction=0.1,
                              structure=None, maxDist=20.0,
                              round_tolerance=False):


    assignments, options, label = find_peak_assignment_options(peak, tolerances, chemShiftRanges,
                                                               aliasing=aliasing, onlyAssignedResonances=onlyAssignedResonances,
                                                               onlyDimensionalAssignmentWhenPresent=onlyDimensionalAssignmentWhenPresent,
                                                               labelling=labelling, minLabelFraction=minLabelFraction,
                                                               structure=structure, maxDist=maxDist)


    if not assignments or not options:
        return None

    new_restraint = OptionalRestraint(peak=peak,
                                      peak_assignment_options=assignments,
                                      restraint_options=options,
                                      labelling=label)

    if round_tolerance:
        new_restraint.apply_round_tolerance(tolerances)

    return new_restraint


def find_peak_assignment_options(peak, tolerances, chemShiftRanges,
                                 aliasing=True, onlyAssignedResonances=True,
                                 onlyDimensionalAssignmentWhenPresent=False,
                                 labelling=None, minLabelFraction=0.1,
                                 structure=None, maxDist=20.0):

    experiment = peak.peakList.dataSource.experiment
    if labelling is True:
        labelling = experiment

    dimensional_options = find_peak_dimensional_assignment_options(peak,
                                                                   tolerances,
                                                                   chemShiftRanges,
                                                                   aliasing,
                                                                   onlyAssignedResonances,
                                                                   onlyDimensionalAssignmentWhenPresent)

    peak_assignment_options = []
    restraint_options = []
    labelling_fractions = []
    for resonances in product(*dimensional_options):
        # 1. experiment
        if not resonances_fit_experiment(resonances, experiment):
            continue

        through_space_pairs = through_space_resonances_combinations(resonances, peak)
        for through_space_pair in through_space_pairs:

            all_resonances = set(resonances + through_space_pair)

            # 2. labelling
            colabelling = 1.0
            if labelling:
                colabelling = getExperimentResonanceSetFractions(labelling, all_resonances)
                if colabelling < minLabelFraction:
                    continue
            # 3. Distance in structure
            if structure:
                if not within_distance_on_structure(structure, through_space_pair, maxDist):
                    continue

            restraint_options.append(through_space_pair)
            peak_assignment_options.append(resonances)
            labelling_fractions.append(colabelling)

    return peak_assignment_options, restraint_options, labelling_fractions


def find_peak_dimensional_assignment_options(peak, tolerances, chemShiftRanges,
                                             aliasing=True,
                                             onlyAssignedResonances=True,
                                             onlyDimensionalAssignmentWhenPresent=False):

    tolerances = dict([(dim, (mintol, maxtol, multi)) for dim, mintol, maxtol, multi in tolerances])

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
                                                         aliasing=aliasing,
                                                         onlyAssignedResonances=onlyAssignedResonances,
                                                         onlyDimensionalAssignmentWhenPresent=onlyDimensionalAssignmentWhenPresent)

        dimensional_options.append(resonances)

    return dimensional_options


def find_dimensional_assignment_options(peakDim, tolerance, chemShiftRange,
                                        aliasing=True, onlyAssignedResonances=True,
                                        onlyDimensionalAssignmentWhenPresent=False):
    '''Find all resonances that can be assigned to a peak dimension.
       If onlyAssignedResonances is True, only resonances are considered
       that have an assignment to an atomset.
       When onlyDimensionalAssignmentWhenPresent is True, no shiftmatching
       is performed if there are already assignments to the peak dimension
       present.

    '''

    tolerance = getPeakDimTolerance(peakDim, *tolerance)
    if onlyDimensionalAssignmentWhenPresent:
        resonances = resonances_assigned_to_peak_dimension(peakDim, onlyAssignedResonances=onlyAssignedResonances)
        if resonances:
            return resonances

    resonances = shift_match_resonances(peakDim, tolerance, chemShiftRange,
                                        aliasing=aliasing,
                                        onlyAssignedResonances=onlyAssignedResonances)
    return resonances


def resonances_assigned_to_peak_dimension(peakDim, onlyAssignedResonances=True):
    '''Get the resonances assigned to a peak dimension.
       If onlyAssignedResonances is True, only resonances
       are considered that have an assignment to an atomset.

    '''

    resonance_assignments = []
    for peakDimContrib in peakDim.peakDimContribs:
        resonance = peakDimContrib.resonance
        resonanceSet = resonance.resonanceSet
        if (resonanceSet and resonanceSet.atomSets) or not onlyAssignedResonances:
            resonance_assignments.append(resonance)
    return resonance_assignments


def shift_match_resonances(peakDim, tolerance, chemShiftRange,
                           aliasing, onlyAssignedResonances=True):
    '''Shiftmatch a peak dimension. If onlyAssignedResonances
       is True, only resonances are considered that have an
       assignment to an atomset.
    '''

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
    '''Returns a dict {isotopeCode:[peakDims]}.'''

    dimensions_by_isotope = defaultdict(list)
    for dim in peak.sortedPeakDims():
        isotope = dim.dataDim.expDim.findFirstExpDimRef().isotopeCodes[0]
        dimensions_by_isotope[isotope].append(dim)
    return dimensions_by_isotope


def peak_has_poor_merit(peak, minMerit):
    return peak.figOfMerit < minMerit


def peak_is_out_of_chemical_shift_ranges(peak, chemShiftRanges):
    '''Whether a peak is within the chemical shift ranges.'''

    range_dict = {}
    for dim, isotope, minshift, maxshift in chemShiftRanges:
        range_dict[dim] = (minshift, maxshift)
    for dim in peak.sortedPeakDims():
        minshift, maxshift = range_dict[dim.dataDim]
        if not minshift < dim.value < maxshift:
            return True
    return False


def peak_is_fully_assigned(peak):
    '''True when there are no unassigned
       peak dimensions.

    '''

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


def within_distance_on_structure(structure, resonances, maxDist):
    '''Check whether the distance between two resonances
       is within the maximum allowed distance.

    '''
    atomSets1 = resonances[0].resonanceSet.atomSets
    atomSets2 = resonances[1].resonanceSet.atomSets
    distance = getAtomSetsDistance(atomSets1, atomSets2,
                                   structure, method='noe')
    if distance > maxDist:
        return False
    else:
        return True


def through_space_resonances_combinations(resonances, peak):
    '''For a peak and a set of resonances that could be
       assigned to the sorted peak dimensions, return
       which pairs of resonances could be actual resonances
       between which magnetization is transferred in the
       through-space transfer (and therefor are the resonances
       restrained). When an indirect transfer is involved,
       these can be resonances not included in the any
       of the dimensional assignments.

    '''

    spectrum = peak.peakList.dataSource
    experiment = spectrum.experiment

    distDataDims = sorted(getThroughSpaceDataDims(spectrum), key=lambda dataDim: dataDim.dim)
    indirectDataDims = [set(dims) for dims in getIndirectDataDims(spectrum)]

    if len(distDataDims) != 2:
        print 'Need 2 through-space data connection, this experiment contains {}.'.format(len(distDataDims))
        return

    if not set(distDataDims) in indirectDataDims:
        return [tuple([resonances[distDataDim.dim-1] for distDataDim in distDataDims])]

    isotopesDict = getIndirectThroughSpaceIsotopes(experiment)

    through_space_resonances = []

    for dataDim in distDataDims:
        expDimRef = dataDim.expDim.sortedExpDimRefs()[0]
        resonance = resonances[dataDim.dim-1]
        indirectIsotope = isotopesDict[expDimRef][1]

        if not indirectIsotope:
            through_space_resonances.append([resonance])
            continue
        through_space_resonances.append(get_bound_resonances_of_isotope(resonance, indirectIsotope))

    return list(product(*through_space_resonances))




def get_bound_resonances_of_isotope(resonance, isotope):
    '''For a resonance, get all resonances of an
       isotope type bound to this resonance.

    '''

    isotopeCode = '{}{}'.format(isotope.massNumber, isotope.chemElement.symbol)

    # Use getBoundResonance to get from e.g. Cga to Hga* and not Hgb*
    resonancesA = set(x for x in getBoundResonances(resonance, recalculate=True)
                      if x.isotopeCode == isotopeCode
                      and x.resonanceSet)

    # get covalently bound atomSts
    atoms = set()
    for atomSet in resonance.resonanceSet.atomSets:
        atoms.update(getBoundAtoms(atomSet.findFirstAtom()))

    atomSets = set(a.atomSet for a in atoms if a.atomSet and \
                   a.chemAtom.chemElement is isotope.chemElement)


    if resonancesA:
        # remove covalently impossible resonances
        resonanceSets = set(y for x in atomSets for y in x.resonanceSets)
        resonancesA = set(x for x in resonancesA
                          if x.resonanceSet in resonanceSets)

    if not resonancesA:
        nmrProject = resonance.parent
        # make new resonances for covanlently bound atoms
        for atomSet in atomSets:
            resonanceB = nmrProject.newResonance(isotopeCode=isotopeCode)
            assignAtomsToRes([atomSet,], resonanceB)
            resonancesA.add(resonanceB)

    return resonancesA


def ambiguity_info(optional_restraints):
    '''prints a list of length 10, representing
       the amount of restraints with ambiguity
       1-9 and >=10.

    '''

    ambiguities = [0]*10

    for restraint in optional_restraints:
        ambiguity = len(restraint.contributions)
        if ambiguity > 10:
            ambiguity = 10
        ambiguities[ambiguity-1] += 1

    #for i, count in enumerate(ambiguities):
    #    print 'ambiguity: {} {}'.format(i+1, count)
    print ambiguities










