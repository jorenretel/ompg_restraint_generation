

from itertools import product


def group_by_residue(atoms):
    '''Groups and iterable of atoms by residue.
       Returns a dict where the keys are molResidues
       and the values are the sets of atoms.

    '''

    residueDict = {}
    for atom in atoms:
        molResidue = atom.residue.molResidue
        if molResidue in residueDict:
            residueDict[molResidue].add(atom)
        else:
            residueDict[molResidue] = set([atom])
    return residueDict


def group_by_molecule(atoms):
    '''Groups and iterable of atoms by molecule.
       Returns a dict where the keys are molecules
       and the values are the sets of atoms.

    '''

    moleculeDict = {}
    for atom in atoms:
        molecule = atom.residue.molResidue.molecule
        if molecule in moleculeDict:
            moleculeDict[molecule].add(atom)
        else:
            moleculeDict[molecule] = set([atom])
    return moleculeDict


def check_atoms_isotopes(atoms, isotopeCodes):
    '''An atom can never be of two isotopes at
       the same time.
    '''

    atomDict = {}
    for atom, isotopeCode in zip(atoms, isotopeCodes):
        if atom in atomDict and atomDict[atom] != isotopeCode:
            return False
    return True


def getAtomsForResonance(resonance):

    atoms = set()
    resonanceSet = resonance.resonanceSet
    atomSets = resonanceSet.atomSets
    for atomSet in atomSets:
        atoms |= set(atomSet.atoms)
    return list(atoms)


def getExperimentResonanceSetFractions(experiment, resonances):
    '''Get the colabelling (simulatanuous labelling) of a set of resonances.
       inputs:  experiment (Nmr.Experiment)
                resonances (iterable of Resonances)
       returns: fraction of colabelling (Float between 0 and 1)

       DISCLAIMER: works only for assigned resonances at the moment.

    '''

    atomGroups = [getAtomsForResonance(resonance) for resonance in resonances]
    isotopes = [resonance.isotopeCode for resonance in resonances]
    fractions = []
    for atoms in product(*atomGroups):
        fractions.append(getExperimentAtomSetFractions(experiment, atoms, isotopes))
    fraction = max(fractions)
    return fraction


def getExperimentAtomSetFractions(experiment, atoms, isotopes):
    '''Get simulatanuous labelling of a set of atoms to specific isotopes.
       inputs:  experiment (Nmr.Experiment)
                atoms (list of Atoms)
                isotopes (list of str isotopeCodes corresponding to
                list of atoms)
       returns: fraction of colabelling (Float between 0 and 1)

    '''

    if not check_atoms_isotopes(atoms, isotopes):
        return 0.0

    atom_isotope_dict = dict((atom, isotope) for atom, isotope in zip(atoms, isotopes))
    atoms_by_molecule = group_by_molecule(atoms)
    labelledMixtures = experiment.labeledMixtures
    colabelling = 1.0

    for molecule, atoms in atoms_by_molecule.items():
        atoms = tuple(atoms)
        isotopes = [atom_isotope_dict[atom] for atom in atoms]
        for mixture in labelledMixtures:
            if mixture.labeledMolecule.molecule is molecule:
                fraction = getMixtureAtomSetFractions(mixture, atoms, isotopes)
                colabelling *= fraction
                break

    return colabelling


def getMixtureAtomSetFractions(labeledMixture, atoms, isotopes):

    molLabelFractions = labeledMixture.molLabelFractions
    molWeightSum = sum([x.weight for x in molLabelFractions])
    colabelling = 0.0

    # Now I want to loop over all different labelling patterns that are
    # present in the sample (denoted with capital letter in the analysis
    # GUI)
    for mlf in molLabelFractions:

        fraction = getMolLabelFractionAtomSetFractions(mlf, atoms, isotopes)
        colabelling += fraction * mlf.weight / float(molWeightSum)

    return colabelling


def getMolLabelFractionAtomSetFractions(molLabelFraction, atoms, isotopes):

    atom_isotope_dict = dict((atom, isotope) for atom, isotope in zip(atoms, isotopes))
    atoms_by_residue = group_by_residue(atoms)
    molLabel = molLabelFraction.molLabel
    colabelling = 1.0

    # Loop over involved residues
    for molResidue, atoms in atoms_by_residue.items():
        atoms = tuple(atoms)
        isotopes = [atom_isotope_dict[atom] for atom in atoms]
        resLabel = molLabel.findFirstResLabel(resId=molResidue.serial)
        fraction = getResLabelAtomSetFractions(resLabel, atoms, isotopes)
        colabelling *= fraction
    return colabelling


def getResLabelAtomSetFractions(resLabel, atoms, isotopes):

    resLabelFractions = resLabel.resLabelFractions
    rlfWeightSum = sum([x.weight for x in resLabelFractions])
    colabelling = 0.0

    for rlf in resLabelFractions:

        fraction = getResLabelFractionAtomSetFractions(rlf, atoms, isotopes)
        colabelling += fraction * rlf.weight / float(rlfWeightSum)

    return colabelling


def getResLabelFractionAtomSetFractions(resLabelFraction, atoms, isotopes):

    isotopomers = resLabelFraction.isotopomers
    isoWeightSum = sum([x.weight for x in isotopomers])

    colabelling = 0.0

    for isotopomer in isotopomers:
        isotopomerLabeling = getIsotopomerAtomSetFractions(isotopomer,
                                                           atoms,
                                                           isotopes)
        colabelling += isotopomerLabeling * isotopomer.weight / float(isoWeightSum)

    return colabelling


def getIsotopomerAtomSetFractions(isotopomer, atoms, isotopes):

    colabelling = 1.0

    for atom, isotope in zip(atoms, isotopes):

        name = atom.name
        subType = atom.chemAtom.subType
        atLabels = isotopomer.findAllAtomLabels(name=name,
                                                subType=subType)
        sumWeight = sum([x.weight for x in atLabels])

        atLabel = isotopomer.findFirstAtomLabel(name=name,
                                                subType=subType,
                                                isotopeCode=isotope)
        atomLabelling = atLabel.weight / float(sumWeight)
        colabelling *= atomLabelling

    return colabelling

