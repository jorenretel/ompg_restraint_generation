
from ccpnmr.analysis.core import ConstraintBasic
from ccpnmr.analysis.core import AssignmentBasic
#import interaction_matrix
#reload(interaction_matrix)
#from interaction_matrix import create_interaction_matrix


def main(argServer):
    project = argServer.getProject()
    #maker = RestraintMaker(project)
    #maker.generate_constraints(['HhNH_406050'])

    restraintFilter = RestraintFilter()
    #restraintFilter.filter(project, 25, 1)
    #restraintFilter.filter(project, 25, 2)
    HHN = getConstraintList(project, 11, 3)
    NNH = getConstraintList(project, 11, 4)
    #filter_constraint_lists_against_each_other([HHN, NNH])
    #make_interaction_matrix([HHN, NNH])
    make_interaction_matrix([HHN, NNH])



class RestraintMaker(object):
    """Making restraints for OmpG"""
    def __init__(self, project):
        super(RestraintMaker, self).__init__()
        self.project = project
        self.nmrProject = project.findFirstNmrProject()
        self.nmrConstraintStore = self.project.newNmrConstraintStore(nmrProject=self.nmrProject)

    def generate_constraints(self, experiment_names):

        for exp_name in experiment_names:
            experiment = self.nmrProject.findFirstExperiment(name=exp_name)
            print experiment


class RestraintFilter(object):

    def __init__(self):

        self.topology = OmpgTopology()

    def filter(self, project, store_serial=None, list_serial=None):

        nmrProject = project.findFirstNmrProject()
        if store_serial:
            nmrConstraintStores = [nmrProject.findFirstNmrConstraintStore(serial=store_serial)]
        else:
            nmrConstraintStores = nmrProject.sortedNmrConstraintStores()

        for store in nmrConstraintStores:
            if store_serial and list_serial:
                constraintLists = [store.findFirstConstraintList(serial=list_serial)]
            else:
                constraintLists = store.sortedConstraintLists()

            for constraintList in constraintLists:
                self.filter_list(constraintList)
                #self.filter_list_symmetry(constraintList)

    def filter_list(self, constraintList):

        print 'Filtering list by topology:'
        print constraintList

        impossible_items = set()
        for constraint in constraintList.sortedConstraints():

            for item in constraint.items:
                residues = get_residue_pair_for_item(item)
                if not self.constraint_is_possible(residues):
                    impossible_items.add(item)
        for item in impossible_items:
            deleteItem(item)

    def filter_list_symmetry(self, constraintList):

        constraint_matrix = make_constraint_matrix(constraintList)
        delete_items = []

        for residues, items in constraint_matrix.items():

            resA, resB = residues

            if not (resB, resA) in constraint_matrix:
                delete_items.extend(items)

        for item in delete_items:
            deleteItem(item)

    def constraint_is_possible(self, residues):

        resA, resB = residues
        return self.topology.residues_in_vicinity(resA.seqCode, resB.seqCode)


def filter_constraint_lists_against_each_other(constraintLists,
                                               full_symmetry=True):

    print 'Filtering constraint lists against each other:'
    print constraintLists

    if not constraintLists:
        return

    delete_items = []
    matrices = []
    for constraintList in constraintLists:
        matrices.append(make_constraint_matrix(constraintList))

    for matrix in matrices:
        for residues, items in matrix.items():
            #if not present_in_all_matrices(residues,
            #                               matrices,
            #                               full_symmetry=full_symmetry):
            #    delete_items.extend(items)
            if not fraction_present_in_matrices(residues, matrices) >= 0.5:
                delete_items.extend(items)

    for item in delete_items:
        deleteItem(item)


def present_in_all_matrices(residue_pair, matrices, full_symmetry=True):

    resA, resB = residue_pair
    reversed_pair = (resB, resA)

    if full_symmetry:
        for matrix in matrices:
            if residue_pair not in matrix or reversed_pair not in matrix:
                return False

    else:
        for matrix in matrices:
            if residue_pair not in matrix and reversed_pair not in matrix:
                return False

    return True

#def make_interaction_matrix(constraintLists):
def select_most_symmetric_contributions(constraintLists, cutoff_fraction):

    print 'selecting most symmetric contributions'
    matrices = []
    residue_combinations = set()
    interaction_matrix = {}
    for constraintList in constraintLists:
        matrix = make_constraint_matrix(constraintList)
        matrices.append(matrix)
        for residues in matrix.keys():
            resA, resB = residues
            if resB.seqCode > resA.seqCode:
                residue_combinations.add(residues)
            else:
                residue_combinations.add((resB, resA))

    for constraintList in constraintLists:
        print constraintList
        for constraint in constraintList.sortedConstraints():
            print '---------'
            print constraint
            print 'constraint: ', constraint.serial
            item_symmetries = []
            print constraint
            for item in constraint.items:
                residues = get_residue_pair_for_item(item)
                a, b = residues
                seqCodes = (a.seqCode, b.seqCode)
                fraction = fraction_present_in_matrices(residues, matrices)
                item_symmetries.append((fraction, item))
                #if fraction > 0.7:
                #    print seqCodes, fraction
            item_symmetries.sort(reverse=True)
            #print item_symmetries
            max_fraction = item_symmetries[0][0]
            if max_fraction >= cutoff_fraction:
                for fraction, item in  item_symmetries:
                    if fraction > max_fraction / 2.0:
                        print 'retaining: ', get_residue_pair_for_item(item)
                        print fraction
                    else:
                        print 'deleting: ', get_residue_pair_for_item(item)
                        print fraction
                        deleteItem(item)
                        #item.weight = 0.0


    #for residues in residue_combinations:
    #    a, b = residues
    #    seqCodes = (a.seqCode, b.seqCode)
    #    fraction = fraction_present_in_matrices(residues, matrices)

    #    interaction_matrix[seqCodes] = fraction

    #for key, fraction in interaction_matrix.items():
    #    if fraction == 1 :
    #        print 4, key

    #    elif fraction >= 0.7 :
    #        print 3, key






def fraction_present_in_matrices(residue_pair, matrices):

    resA, resB = residue_pair
    reversed_pair = (resB, resA)
    total_possible = float(2*len(matrices))    # A sees B and B sees A
    number = 0

    for matrix in matrices:
        if residue_pair in matrix:
            number += 1
        if reversed_pair in matrix:
            number += 1
    #if number > 2:
    #    print residue_pair, number
    return number / total_possible


def make_constraint_matrix(constraintList):

    constraint_matrix = {}
    for constraint in constraintList.sortedConstraints():
        for item in constraint.items:
            residues = get_residue_pair_for_item(item)
            if residues in constraint_matrix:
                constraint_matrix[residues].append(item)
            else:
                constraint_matrix[residues] = [item]
    return constraint_matrix


def getConstraintList(project, store_serial, list_serial):

    nmrProject = project.findFirstNmrProject()
    print nmrProject.nmrConstraintStores
    nmrConstraintStore = nmrProject.findFirstNmrConstraintStore(serial=store_serial)
    constraintList = nmrConstraintStore.findFirstConstraintList(serial=list_serial)
    return constraintList


class StructureElement(object):
    """docstring for StructureElement"""
    def __init__(self, topology_type, start, stop):
        super(StructureElement, self).__init__()
        self.start = start
        self.stop = stop
        self.topology_type = topology_type
        self.neighboring = set()

    @property
    def residue_numbers(self):
        return range(self.start, self.stop+1)

    def __str__(self):
        return 'structural element {}, residues {} to {}'.format(self.topology_type,
                                                                 self.start,
                                                                 self.stop)


class OmpgTopology(object):

    def __init__(self):
        se = StructureElement
        self.domains = [se('in', 1, 6),
                        se('tm', 7, 15),
                        se('out', 16, 28),
                        se('tm', 29, 37),
                        se('in', 38, 43),
                        se('tm', 44, 50),
                        se('out', 51, 65),
                        se('tm', 66, 78),
                        se('in', 79, 83),
                        se('tm', 84, 92),
                        se('out', 93, 112),
                        se('tm', 113, 123),
                        se('in', 124, 126),
                        se('tm', 127, 135),
                        se('out', 136, 151),
                        se('tm', 152, 162),
                        se('in', 163, 165),
                        se('tm', 166, 178),
                        se('out', 179, 191),
                        se('tm', 192, 202),
                        se('in', 203, 204),
                        se('tm', 205, 215),
                        se('out', 216, 234),
                        se('tm', 235, 245),
                        se('in', 246, 248),
                        se('tm', 249, 259),
                        se('out', 260, 270),
                        se('tm', 271, 281)]

        for i, domain in enumerate(self.domains):

            Ndomains = len(self.domains)

            if domain.topology_type == 'in':
                neighbors = [-4, -3, -1, 1, 3, 4]
                domain.neighboring = set([self.domains[(i+n)%Ndomains] for n in neighbors])
            elif domain.topology_type == 'tm':
                neighbors = [-3, -2, -1, 1, 2, 3]
                domain.neighboring = set([self.domains[(i+n)%Ndomains] for n in neighbors])
            elif domain.topology_type == 'out':
                neighbors = [-4, -3, -1, 1, 3, 4]
                domain.neighboring = set([self.domains[(i+n)%Ndomains] for n in neighbors])

    def domain_for_residue_number(self, residue_number):

        for domain in self.domains:
            if domain.start <= residue_number <= domain.stop:
                return domain

    def residues_in_vicinity(self, residue_number1, residue_number2):

        delta = abs(residue_number1 - residue_number2)

        if delta == 0 or delta == 1:
            return True

        domainA = self.domain_for_residue_number(residue_number1)
        domainB = self.domain_for_residue_number(residue_number2)

        if domainA in domainB.neighboring:
            return True
        else:
            return False


def deleteItem(item):

    constraintSet = item.constraint.parent.nmrConstraintStore
    project = constraintSet.parent
    project.__dict__['override'] = True
    constraint = item.constraint
    constraint_got_deleted = False
    item.delete()
    if not constraint.items:
        constraint.delete()
        constraint_got_deleted = True

    constraintSet.checkAllValid()
    constraintSet.__dict__['isModified'] = True

    for func in item._notifies.get('delete', []):
        func(item)

    if constraint_got_deleted:
        for func in constraint._notifies.get('delete', []):
            func(constraint)


def get_residue_pair_for_item(item):
    '''This function for now assumes
       that one of the dimensions of the experiment is the one
       involved in through space transfer, which is not always true.

    '''

    #print item

    constraint = item.constraint
    peak = constraint.findFirstPeak()
    acqDim = get_aquired_dim(peak)
    #peaksDims = peak.sortedPeakDims()
    #shiftList = peak.peakList.dataSource.experiment.shiftList
    #residues = []

    resonanceA, resonanceB = [resonance.resonance for resonance in item.resonances]

    #print AssignmentBasic.makeResonanceGuiName(resonanceA)
    #print AssignmentBasic.makeResonanceGuiName(resonanceB)

    acquired_isotope = get_isotope(acqDim)
    acquired_shift = acqDim.value

    if not resonanceA.findFirstShift():
        ordered_resonances = [resonanceA, resonanceB]
    elif not resonanceB.findFirstShift():
        ordered_resonances = [resonanceB, resonanceA]

    elif resonanceA.isotopeCode == acquired_isotope and resonanceB.isotopeCode == acquired_isotope:

        deltaA = abs(resonanceA.findFirstShift().value - acquired_shift)
        deltaB = abs(resonanceB.findFirstShift().value - acquired_shift)
        if deltaA <= deltaB:
            ordered_resonances = [resonanceB, resonanceA]
        else:
            ordered_resonances = [resonanceA, resonanceB]

    elif resonanceA.isotopeCode == acquired_isotope:
        ordered_resonances = [resonanceB, resonanceA]

    elif resonanceA.isotopeCode == acquired_isotope:
        ordered_resonances = [resonanceA, resonanceB]
    else:
        raise NotImplementedError



    #dimensionsA = []
    #dimensionsB = []
    #for peakDim in peaksDims:
    #    dim_isotope = get_isotope(peakDim)
    #    if dim_isotope == resonanceA.isotopeCode:
    #        delta = abs(peakDim.value - resonanceA.findFirstShift().value)
    #        dimensionsA.append((delta, peakDim.dim))
    #    if dim_isotope == resonanceB.isotopeCode:
    #        delta = abs(peakDim.value - resonanceB.findFirstShift().value)
    #        dimensionsB.append((delta, peakDim.dim))
    #dimensionsA.sort()
    #dimensionsB.sort()
    #favA = dimensionsA[0][1]
    #favB = dimensionsB[0][1]
    #if favA == favB:
    #    if dimensionsA[0][0] <= dimensionsB[0][0]:
    #        favB = dimensionsB[1][1]     # second choice
    #    else:
    #        favA = dimensionsA[1][1]
    #ordered_resonances = [tup[1] for tup in sorted([(favA, resonanceA), (favB, resonanceB)])]

    residues = [resonance.getResonanceGroup().residue for resonance in ordered_resonances]

    return tuple(residues)


def get_isotope(peakDim):
    isotope = peakDim.dataDimRef.expDimRef.findFirstIsotope()
    code = '{}{}'.format(isotope.massNumber, isotope.chemElement.symbol)
    return code

def get_aquired_dim(peak):
    for peakDim in peak.sortedPeakDims():
        if peakDim.dataDim.expDim.isAcquisition:
            return peakDim


if __name__ == '__main__':
    main(None)
