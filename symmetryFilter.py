#pylint: disable=invalid-name

def symmetry_filter(optionalRestraintLists, max_symmetry=None,
                    cutoff_fraction=1.0):
    '''Use multiple restraint lists based on different experiments
       and search for contributions that are supported by multiple
       peaks.
       optionalRestraintLists: list of lists of
                     optionalRestraint.OptionalRestraint
       max_symmetry: integer describing the maximum amount of
                     peaks/restraints describe the same interaction
                     between two resonances. This depends a bit on
                     the experiment. In a hNHH on a deuterated and
                     backexchanged protein for instance this
                     would be 2 (two peaks for every interaction
                     between two amide groups). In the combination
                     of an hNhhNH and hNHH, this would be 4.
                     If set to None, the amount of restraint sets
                     (experiments) times 2 is used.
       cutoff_fraction: minimum fraction of the symmetry for which
                     a contribution of a restraint causes less
                     likely contributions to be removed. 1.0
                     means that only contributions with
                     max_symmetry cause this to happen (i.e.
                     every single expected peak for the
                     interaction is there).

    '''

    if max_symmetry is None:
        max_symmetry = float(len(optionalRestraintLists) * 2)

    interaction_sets = create_interaction_sets(optionalRestraintLists)

    for optionalRestraintList in optionalRestraintLists:
        for optionalRestraint in optionalRestraintList:
            remove_less_symmetric_options(optionalRestraint, interaction_sets,
                                          max_symmetry, cutoff_fraction)



def create_interaction_sets(optionalRestraintLists):
    '''Make a set of tuples (resonance, resonance) for each
       optionalRestraintList, containing all interactions between
       resonances in the restraint set.

    '''

    interaction_sets = []
    for optionalRestraintList in optionalRestraintLists:
        interaction_set = set()
        for optionalRestraint in optionalRestraintList:
            for contrib in optionalRestraint. contributions:
                interaction_set.add(contrib.restraint_resonances)
        interaction_sets.append(interaction_set)
    return interaction_sets


def remove_less_symmetric_options(optionalRestraint, interaction_sets,
                                  max_symmetry, cutoff_fraction):
    '''Removing restaint items that have lower symmetry.'''


    contribs = optionalRestraint.contributions
    for contrib in contribs:
        option = contrib.restraint_resonances
        for interaction_set in interaction_sets:
            if option in interaction_set:
                contrib.symmetry += 1
            if (option[1], option[0]) in interaction_set:
                contrib.symmetry += 1

    max_symmetry_fraction = max([c.symmetry for c in contribs])/max_symmetry

    if max_symmetry_fraction < cutoff_fraction:
        return

    for contrib in contribs[:]:

        symmetry_fraction = contrib.symmetry / max_symmetry
        if symmetry_fraction <= max_symmetry_fraction/ 2.0:
            contribs.remove(contrib)
