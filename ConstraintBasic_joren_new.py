
from ccpnmr.analysis.core.ConstraintBasic import *
from filter_restraints_ompg import select_most_symmetric_contributions
import labeling_simple
reload(labeling_simple)
from labeling_simple import getExperimentResonanceSetFractions


def main_test(argServer):
    project = argServer.getProject()
    nmrProject = project.findFirstNmrProject()
    nmrConstraintStore = project.newNmrConstraintStore(nmrProject=nmrProject)

    #HHN
    expName = '248_07HHN'
    peakListSerial = 2
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)
    #distanceFunction = DistanceFunctionProtonDetected(peakList, [(0.25, 3.1, 1.0, 3.5), (0.0, 4.4, 1.0, 5.5)])
    distanceFunction = DistanceFunctionProtonDetected(peakList, [(1500000, 3.1, 1.0, 3.5), (0.0, 4.4, 1.0, 5.5)])

    dataDim1, dataDim2, dataDim3 = spectrum.sortedDataDims()
    tolerances = [(dataDim1, 0.07, 0.07, 1.0), (dataDim2, 0.1, 0.1, 1.0), (dataDim3, 0.4, 0.4, 1.0)]
    chemShiftRanges = [(dataDim1, '1H', 0.0,13.0), (dataDim2, '1H', 0.0, 13.0),(dataDim3, '15N', 100.0, 140.0)]
    HHN = makeAmbigDistConstraints_modified(peakList, tolerances,
                                            chemShiftRanges,
                                            constraintSet=nmrConstraintStore,
                                            testOnly=False,
                                            distanceFunction=distanceFunction,
                                            round_tolerances=False)
    toleranceString = str(tolerances[0][1]) + '_' + str(tolerances[1][1]) + '_' + str(tolerances[2][1])
    HHN.name = 'HHN_peaklist{}_shiftmatch_{}'.format(peakListSerial, toleranceString)
    HHN.details = 'symmetry_filtered'

    #NNH
    expName = '250_HNN_30_10_sp035'
    peakListSerial = 2
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)
    #distanceFunction = DistanceFunctionProtonDetected(peakList, [(0.25, 3.1, 1.0, 3.5), (0.0, 4.4, 1.0, 5.5)])
    distanceFunction = DistanceFunctionProtonDetected(peakList, [(2500000, 3.1, 1.0, 3.5), (0.0, 4.4, 1.0, 5.5)])
    dataDim1, dataDim2, dataDim3 = spectrum.sortedDataDims()
    tolerances = [(dataDim1, 0.07, 0.07, 1.0), (dataDim2, 0.4, 0.4, 1.0), (dataDim3, 0.4, 0.4, 1.0)]
    chemShiftRanges = [(dataDim1, '1H', 0.0,13.0), (dataDim2, '15N', 100.0, 140.0),(dataDim3, '15N', 100.0, 140.0)]
    NNH = makeAmbigDistConstraints_modified(peakList, tolerances,
                                            chemShiftRanges,
                                            constraintSet=nmrConstraintStore,
                                            testOnly=False,
                                            distanceFunction=distanceFunction,
                                            round_tolerances=False)
    toleranceString = str(tolerances[0][1]) + '_' + str(tolerances[1][1]) + '_' + str(tolerances[2][1])
    NNH.name = 'NNH_peaklist{}_shiftmatch_{}'.format(peakListSerial, toleranceString)
    NNH.details = 'symmetry_filtered'


    simple_shift_match_pdsd(nmrProject, '029_PDSD400_13C', 4, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4)
    simple_shift_match_pdsd(nmrProject, '063_PDSD400_2C', 5, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4)
    simple_shift_match_pdsd(nmrProject, '045_PDSD400_13TEMPQANDSG', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4)
    simple_shift_match_pdsd(nmrProject, '113_PDSD400_2TEMPQANDSG', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4)
    simple_shift_match_pdsd(nmrProject, '110_PDSD400_2SLYGWAFVL', 2, nmrConstraintStore, (4.0, 1.5, 7.0), 0.4)


    # keep those contributions that appear most
    # symmetrically in spectra.
    #select_most_symmetric_contributions([HHN, NNH], 1.0)


def simple_shift_match_pdsd(nmrProject, expName, peakListSerial,
                            nmrConstraintStore, bucket, tolerance):

    # One length distance bucket
    dist_function_13c = lambda x: bucket
    experiment = nmrProject.findFirstExperiment(name=expName)
    spectrum = experiment.findFirstDataSource()
    peakList = spectrum.findFirstPeakList(serial=peakListSerial)
    dataDim1, dataDim2 = spectrum.sortedDataDims()
    tol = tolerance
    tolerances = [(dataDim1, tol, tol, 1.0), (dataDim2, tol, tol, 1.0)]
    chemShiftRanges = [(dataDim1, '13C', 0.0, 200.0), (dataDim2, '13C', 0.0, 200.0)]
    constraints = makeAmbigDistConstraints_modified(peakList, tolerances,
                                                    chemShiftRanges,
                                                    constraintSet=nmrConstraintStore,
                                                    labelling=True,
                                                    distanceFunction=dist_function_13c,
                                                    testOnly=False,
                                                    round_tolerances=False,
                                                    minLabelFraction=0.25)
    toleranceString = str(tolerances[0][1]) + '_' + str(tolerances[1][1])
    constraints.name = '{}_{}_shiftmatch_{}'.format(expName, peakListSerial, toleranceString)

    return constraints


class DistanceFunctionProtonDetected(object):

    def __init__(self, peakList, distance_bins, intensityType='height'):

        self.peakList = peakList
        self.intensityType = intensityType
        peak_intensities = [abs(peak.findFirstPeakIntensity(intensityType=intensityType).value) for peak in peakList.peaks]
        max_intensity = max(peak_intensities)
        min_intensity = min(peak_intensities)

        distance_bins.sort(reverse=True)
        self.intensity_distance_bins = []

        for fraction, dist, minDist, maxDist in distance_bins:

            if fraction > 1.0:
                # If the fraction is bigger than 1, it is not a fraction
                # but a specified cutoff intensity.
                threshold_intensity = fraction
            else:
                threshold_intensity = fraction * (max_intensity - min_intensity) + min_intensity

            self.intensity_distance_bins.append((threshold_intensity, dist, minDist, maxDist))

        print 'threshhold'
        print peakList
        print self.intensity_distance_bins

    def __call__(self, peak):

        intensity = peak.findFirstPeakIntensity(intensityType=self.intensityType).value

        for threshold, dist, minDist, maxDist in self.intensity_distance_bins:
            if intensity >= threshold:
                return (dist, minDist, maxDist)



def makeAmbigDistConstraints_modified(peakList, tolerances, chemShiftRanges, constraintSet=None,
                             testOnly=False, labelling=None, minLabelFraction=0.1,
                             distanceFunction=None, residueRanges=None, minMerit=0.0, progressBar=None,
                             intensityType='volume', ignoreDiagonals=True, doAliasing=True,
                             structure=None, maxDist=None,
                             scale=None, params=None, peakCategories=None,
                             round_tolerances=False):
  """
  Makes a constraint list with constraints by matching known shifts
  in given range and within specified tolerances to a NOESY peak
  list. Optional labelling scheme/mixture and labelling threshold to filter
  according to a given set of residue isotopomers.
  Constraints will be put in a new NMR constraint store object
  if none is specified. Peaks are catergorised into various lists.
  A structure and max distance can be used to filter contributions.
  The scale option is the value by which peak intensities are scaled
  for NOE table/function lookup. Params relate to the generic NOE
  distance function if neither these nor a distance function is
  specified a lookuptable is used

  .. describe:: Input

  .PeakList (NOESY), List of (Nmr.DataDim, Float, Float, Float) (shift tolerances),
  List of (Nmr.DataDim, String (isotope code), Float, Float) (chem shift ranges)
  Nmr.NmrConstraintStore, Boolean (test only),
  ChemCompLabel.LabelingScheme or True (automatic from experiment MolLabel), Float,
  Function (to get distances from NOEs),  Nmr.PeakIntensity.intensityType,
  List of (List of Nmr.DataDims, Nmr.Chain, Integer, Integer) (residue ranges),
  Float (min Peak.figOfMerit), ProgressBar (Analysis popup),
  Boolean, Boolean, MolStructure.StructureEnsemble, Float
  Float, List of Floats (distance function parameters)
  Dict (for Category Name:Lists of Nmr.Peaks)

  .. describe:: Output

  NmrConstraint.DistanceConstraintList
  """

  from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance

  if peakCategories is None:
    peakCategories = {}

  assignedPeaks = peakCategories['Assigned'] = []
  diagonalPeaks = peakCategories['Diagonal'] = []
  unmatchedPeaks = peakCategories['Unmatchable'] = []
  poorMeritPeaks = peakCategories['Poor Merit'] = []
  outOfRangePeaks = peakCategories['Out of range'] = []
  distalPeaks = peakCategories['Too Distal'] = []

  #peaks = peakList.peaks
  spectrum = peakList.dataSource
  experiment = spectrum.experiment
  nmrProject = experiment.nmrProject
  distDataDims = getThroughSpaceDataDims(spectrum)
  distIndices  = [dd.dim-1 for dd in distDataDims]

  if len(distDataDims) != 2:
    return

  distDim1, distDim2 = distDataDims

  if labelling is True:
    labelling = experiment

  if not residueRanges:
    residueRanges = None

  bondedDims = {}
  for dataDim1, dataDim2 in getOnebondDataDims(spectrum):
    bondedDims[dataDim1] = dataDim2
    bondedDims[dataDim2] = dataDim1

  if testOnly:
    distConstraintList = None
  else:
    if not constraintSet:
      constraintSet = makeNmrConstraintStore(experiment.topObject)
    distConstraintList = constraintSet.newDistanceConstraintList()
    distConstraintList.addExperimentSerial(experiment.serial)
    newConstraint = distConstraintList.newDistanceConstraint

  tolDict = {}
  for (dataDim,minT,maxT,multi) in tolerances:
    tolDict[dataDim] = (minT,maxT,multi)

  chemShiftRangesDict = {}
  for (dataDim, iso, minShift, maxShift) in chemShiftRanges:
    if chemShiftRangesDict.get(dataDim) is None:
      chemShiftRangesDict[dataDim] = []

    chemShiftRangesDict[dataDim].append([minShift, maxShift])
    if tolDict.get(dataDim) is None:
      msg = 'No tolerance set for dataDim %s of dataSource %s' % (dataDim,spectrum)
      raise Exception(msg)

  # go through peaks
  # if not assigned in all the Hydrogen dims or H + bonded dim

  # Check for indirect transfers
  indirectDims = {}
  for dataDims in getIndirectDataDims(spectrum):
    if set(dataDims) == set(distDataDims):
      isotopesDict = getIndirectThroughSpaceIsotopes(experiment)

      for dataDim in dataDims:
        expDimRef = dataDim.expDim.sortedExpDimRefs()[0]
        indirectDims[dataDim.dim] = isotopesDict[expDimRef]

  workingPeaks = []
  for peak in peakList.sortedPeaks():
    # filter out diagonals
    if ignoreDiagonals:
      peakDims = peak.sortedPeakDims()
      peakDim1 = peakDims[distIndices[0]]
      peakDim2 = peakDims[distIndices[1]]
      ppm1 = peakDim1.realValue
      ppm2 = peakDim2.realValue

      delta = abs(ppm1-ppm2)

      if (delta <= tolDict[distDim1][0] ) or (delta <= tolDict[distDim2][0]):
        dataDimA = bondedDims.get(distDim1)
        dataDimB = bondedDims.get(distDim2)

        if dataDimA and dataDimB :
          peakDimA = peak.findFirstPeakDim(dataDim=dataDimA)
          peakDimB = peak.findFirstPeakDim(dataDim=dataDimB)
          ppmA = pnt2ppm(peakDimA.position,peakDimA.dataDimRef)
          ppmB = pnt2ppm(peakDimB.position,peakDimB.dataDimRef)

          delta2 = abs(ppmA-ppmB)
          if (delta2 <= tolDict[dataDimA][0] ) or (delta2 <= tolDict[dataDimB][0]):
            diagonalPeaks.append(peak)
            continue

        else:
          diagonalPeaks.append(peak)
          continue

    if peak.figOfMerit < minMerit:
      poorMeritPeaks.append(peak)
      continue

    workingPeaks.append(peak)

  mean = getMeanPeakIntensity(workingPeaks, intensityType=intensityType)
  if scale: # neither Zero nor None
    mean = scale

  if not mean:
    msg  = 'Cannot make restraints: peak %s is zero on average.' % intensityType
    msg += ' Maybe intensities are missing or the peak list is empty'
    showWarning('Failure', msg)
    return


  if not distanceFunction:
    if params:
      distanceFunction = lambda val:getNoeDistance(val, params)
    else:
      noeDistClasses = getIntensityDistanceTable(spectrum)
      distanceFunction = lambda val:getDistancesFromIntensity(noeDistClasses,val)

  for peak in workingPeaks:
    #print 'newpeak', peak

    if progressBar:
      progressBar.increment()

    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if not intensity:
      continue
    intensityValue = abs(intensity.value)

    outOfRange = 0

    unassignedPeakDims = []
    peakDims = peak.sortedPeakDims()

    outOfShiftRange = False
    for peakDim in peakDims:
      inRange = isShiftInRange(peakDim.realValue,chemShiftRangesDict[peakDim.dataDim])
      if chemShiftRanges and not inRange:
        outOfShiftRange = True
        break


    if outOfShiftRange:
      unmatchedPeaks.append(peak)
      continue

    #n = 0
    for i in distIndices:
      peakDim = peakDims[i]
      if not peakDim.peakDimContribs:
        unassignedPeakDims.append( peakDim )


    # filter out assigned peaks
    if not unassignedPeakDims:
      assignedPeaks.append(peak)
      continue

    peakResonances = []
    for i in distIndices:
      resonances = []
      peakDim = peakDims[i]
      dataDim = peakDim.dataDim

      if peakDim in unassignedPeakDims:
        #isotope    = dataDim.expDim.findFirstExpDimRef().isotopeCodes[0]
        #peakDimPos = peakDim.position + (peakDim.numAliasing*dataDim.numPointsOrig)
        (minT,maxT,multi) = tolDict[dataDim]
        tolerance  = getPeakDimTolerance(peakDim,minT,maxT,multi)

        bondedDim = bondedDims.get(peakDim.dataDimRef.dataDim)

        if bondedDim:
          # check that both bonded dim possibilities are within tolerances

          shifts = findMatchingPeakDimShifts(peakDim,
                                             chemShiftRangesDict[dataDim],
                                             tolerance=tolerance,
                                             aliasing=doAliasing,
                                             findAssigned=True)

          if shifts:
            for peakDim2 in peakDims:
              if peakDim2.dataDimRef.dataDim is bondedDim:

                shifts2 = []
                if peakDim2.peakDimContribs:
                  for contrib in peakDim2.peakDimContribs:
                    shift = contrib.resonance.findFirstShift(parentList=experiment.shiftList)
                    if shift:
                      shifts2.append(shift)

                else:
                  dataDim2    = peakDim2.dataDim
                  (minT,maxT,multi) = tolDict[dataDim2]
                  tolerance2  = getPeakDimTolerance(peakDim2,minT,maxT,multi)

                  shifts2 = findMatchingPeakDimShifts(peakDim2,
                                                      chemShiftRangesDict[dataDim2],
                                                      tolerance=tolerance2,
                                                      aliasing=doAliasing,
                                                      findAssigned=True)

                for shift in shifts:
                  resonance = shift.resonance
                  shiftDeviation = abs(shift.value - peakDim.value)

                  for shift2 in shifts2:
                    resonance2 = shift2.resonance
                    shiftDeviation2 = abs(shift2.value - peakDim2.value)
                    eucledianDist = ((shiftDeviation/tolerance)**2 + (shiftDeviation2/tolerance2)**2)**0.5

                    if areResonancesBound(resonance, resonance2):
                      if residueRanges:
                        residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
                        if isResidueInRange(residue, residueRanges, dataDim):
                          resonances.append((resonance, resonance2, [], eucledianDist))
                        else:
                          outOfRange += 1
                      else:
                        #print eucledianDist
                        resonances.append((resonance, resonance2, [], eucledianDist))
                      break

        else:

          shifts = findMatchingPeakDimShifts(peakDim,
                                             chemShiftRangesDict[dataDim],
                                             tolerance=tolerance,
                                             aliasing=doAliasing,
                                             findAssigned=True)

          for shift in shifts:
            resonance = shift.resonance

            if residueRanges:
              residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
              if isResidueInRange(residue, residueRanges, dataDim):
                resonances.append((resonance, None, [], 0.0))
              else:
                outOfRange += 1
            else:
              resonances.append((resonance, None, [], 0.0))

      else:
        # this dim is assigned
        for contrib in peakDim.peakDimContribs:
          resonance = contrib.resonance
          resonanceSet = resonance.resonanceSet

          if resonanceSet:
            if residueRanges:
              residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
              if isResidueInRange(residue, residueRanges, dataDim):
                resonances.append((resonance, None, [], 0.0))
              else:
                outOfRange += 1
            else:
              resonances.append((resonance, None, [], 0.0))

      # Deal with indirect transfers

      if peakDim.dim in indirectDims and indirectDims[peakDim.dim][1]:
        isotopeA, isotopeB = indirectDims[peakDim.dim]
        chemElement = isotopeB.chemElement

        for resonance, bound, indirect, eucledianDist in resonances:

          isotopeCode = '%d%s' % (isotopeB.massNumber, chemElement.symbol)

          # Use getBoundResonance to get from e.g. Cga to Hga* and not Hgb*
          resonancesA = set(x for x in getBoundResonances(resonance, recalculate=True)
                            if x.isotopeCode == isotopeCode
                            and x.resonanceSet)

          # get covalently bound atomSts
          atoms = set()
          for atomSet in resonance.resonanceSet.atomSets:
            atoms.update(getBoundAtoms(atomSet.findFirstAtom()))

          atomSets = set(a.atomSet for a in atoms if a.atomSet and \
                         a.chemAtom.chemElement is chemElement)


          if resonancesA:
            # remove covalently impossible resonances
            resonanceSets = set(y for x in atomSets for y in x.resonanceSets)
            resonancesA = set(x for x in resonancesA
                              if x.resonanceSet in resonanceSets)

          if not resonancesA:
            # make new resonances for covanlently bound atoms
            for atomSet in atomSets:
              resonanceB = nmrProject.newResonance(isotopeCode=isotopeCode)
              assignAtomsToRes([atomSet,], resonanceB)
              resonancesA.add(resonanceB)

          indirect.extend(resonancesA)

      # Store resonances for this dim
      if round_tolerances:
        resonances = select_closest_resonance_pair(resonances)
      else:
        resonances = [res[:3] for res in resonances]
        #for resonance, bound, indirect, eucledianDist in resonances:


      peakResonances.append(resonances)

    if peakResonances[0] and peakResonances[1]:
      distal = False

      resonancePairs = set()
      for resonance0, bound0, indirect0 in peakResonances[0]:
        resonanceSet0 = resonance0.resonanceSet

        for resonance1, bound1, indirect1 in peakResonances[1]:
          #print resonance0, resonance1, bound0, bound1, indirect0, indirect1
          if resonance1 is resonance0:
            continue

          if labelling:

            all_resonances = [resonance0, resonance1, bound0, bound1] + indirect0 + indirect1
            all_resonances = [res for res in all_resonances if res]

            colabelling = getExperimentResonanceSetFractions(labelling, all_resonances)

            if colabelling < minLabelFraction:
              continue


          if structure and resonanceSet0 and (maxDist is not None):
            resonanceSet1 = resonance1.resonanceSet

            if resonanceSet1:
              atomSets0 = list(resonanceSet0.atomSets)
              atomSets1 = list(resonanceSet1.atomSets)
              dist = getAtomSetsDistance(atomSets0, atomSets1, structure, method='noe')

              if dist > maxDist:
                distal = True
                continue

          resonances0 = indirect0 or [resonance0,]
          resonances1 = indirect1 or [resonance1,]

          for resonanceA in resonances0:
            for resonanceB in resonances1:
              if resonanceA is not resonanceB:
                resonancePairs.add(frozenset([resonanceA, resonanceB]))

      if not resonancePairs:
        unmatchedPeaks.append(peak)

        if distal:
          distalPeaks.append(peak)

      elif not testOnly:
        resonances0 = [(resonance,indirect) for resonance, bound, indirect in peakResonances[0]]
        resonances1 = [(resonance,indirect) for resonance, bound, indirect in peakResonances[1]]
        #(dist,minDistL,maxDistL) = getDistMinMax(intensityValue, mean, resonances0, resonances1, distanceFunction, labelling=labelling)
        dist,minDistL,maxDistL = distanceFunction(peak)
        #dist = 3.2
        #minDistL = 1.0
        #maxDistL = 5.0
        error = abs(maxDistL - minDistL)
        constraint  = newConstraint(weight=1.0, origData=intensityValue,
                                    targetValue=dist, upperLimit=maxDistL,
                                    lowerLimit=minDistL, error=error)

        constraint.newConstraintPeakContrib(experimentSerial=experiment.serial,
                                            dataSourceSerial=spectrum.serial,
                                            peakListSerial=peakList.serial,
                                            peakSerial=peak.serial)

        for resonance0, resonance1 in resonancePairs:
          fixedResonance0 = getFixedResonance(constraintSet,resonance0)
          fixedResonance1 = getFixedResonance(constraintSet,resonance1)
          constraint.newDistanceConstraintItem(resonances=[fixedResonance0,fixedResonance1])

    elif outOfRange > 1:
      outOfRangePeaks.append(peak)
    else:
      unmatchedPeaks.append(peak)

  return distConstraintList



def select_closest_resonance_pair(resonances):
    '''Select only those resonances in a round restraint bucket where the
       distance from the peak to the expected peak position is either
       smaller than:
           2 * the closest assignment option
       or
           half the radius of the (normalized) bucket (sqrt(2)/2)
       the second rule is to prevent throwing away close options just
       because the very closest option is just really close.

    '''
    print 'round restraints'
    closest_resonances = []
    if not resonances:
        return resonances
    resonances = sorted(resonances,key=lambda x: x[3])
    closest = max(resonances[0][3], 2**0.5/4)

    for resonance, bound, indirect, eucledianDist in resonances:

        #print eucledianDist

        if eucledianDist <= 2*closest:
            closest_resonances.append((resonance, bound, indirect))

    #print len(resonances)
    #print len(closest_resonances)
    return closest_resonances
