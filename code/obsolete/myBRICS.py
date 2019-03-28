#
# rdkit version: 2018.09.1dev1 / 2018.09.1pre.132


#  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Greg Landrum, Nov 2008
""" Implementation of the BRICS algorithm from Degen et al. ChemMedChem *3* 1503-7 (2008)

"""
from __future__ import print_function
import sys, re, random
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.six import iteritems, iterkeys, next
# from rdkit.six.moves import range

# These are the definitions that will be applied to fragment molecules:
environs = {
    'L1': '[C;D3]([#0,#6,#7,#8])(=O)',
    #
    # After some discussion, the L2 definitions ("N.pl3" in the original
    # paper) have been removed and incorporated into a (almost) general
    # purpose amine definition in L5 ("N.sp3" in the paper).
    #
    # The problem is one of consistency.
    #    Based on the original definitions you should get the following
    #    fragmentations:
    #      C1CCCCC1NC(=O)C -> C1CCCCC1N[2*].[1*]C(=O)C
    #      c1ccccc1NC(=O)C -> c1ccccc1[16*].[2*]N[2*].[1*]C(=O)C
    #    This difference just didn't make sense to us. By switching to
    #    the unified definition we end up with:
    #      C1CCCCC1NC(=O)C -> C1CCCCC1[15*].[5*]N[5*].[1*]C(=O)C
    #      c1ccccc1NC(=O)C -> c1ccccc1[16*].[5*]N[5*].[1*]C(=O)C
    #
    # 'L2':'[N;!R;!D1;!$(N=*)]-;!@[#0,#6]',
    # this one turned out to be too tricky to define above, so we set it off
    # in its own definition:
    # 'L2a':'[N;D3;R;$(N(@[C;!$(C=*)])@[C;!$(C=*)])]',
    'L3': '[O;D2]-;!@[#0,#6,#1]',
    'L4': '[C;!D1;!$(C=*)]-;!@[#6]',
    # 'L5':'[N;!D1;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]',
    'L5': '[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]',
    'L6': '[C;D3;!R](=O)-;!@[#0,#6,#7,#8]',
    'L7a': '[C;D2,D3]-[#6]',
    'L7b': '[C;D2,D3]-[#6]',
    '#L8': '[C;!R;!D1]-;!@[#6]',
    'L8': '[C;!R;!D1;!$(C!-*)]',
    'L9': '[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]',
    'L10': '[N;R;$(N(@C(=O))@[C,N,O,S])]',
    'L11': '[S;D2](-;!@[#0,#6])',
    'L12': '[S;D4]([#6,#0])(=O)(=O)',
    'L13': '[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]',
    'L14': '[c;$(c(:[c,n,o,s]):[n,o,s])]',
    'L14b': '[c;$(c(:[c,n,o,s]):[n,o,s])]',
    'L15': '[C;$(C(-;@C)-;@C)]',
    'L16': '[c;$(c(:c):c)]',
    'L16b': '[c;$(c(:c):c)]',
}
reactionDefs = (
    # L1
    [
        ('1', '3', '-'),
        ('1', '5', '-'),
        ('1', '10', '-'),
    ],

    # L3
    [
        ('3', '4', '-'),
        ('3', '13', '-'),
        ('3', '14', '-'),
        ('3', '15', '-'),
        ('3', '16', '-'),
    ],

    # L4
    [
        ('4', '5', '-'),
        ('4', '11', '-'),
    ],

    # L5
    [
        ('5', '12', '-'),
        ('5', '14', '-'),
        ('5', '16', '-'),
        ('5', '13', '-'),
        ('5', '15', '-'),
    ],

    # L6
    [
        ('6', '13', '-'),
        ('6', '14', '-'),
        ('6', '15', '-'),
        ('6', '16', '-'),
    ],

    # L7
    [
        ('7a', '7b', '='),
    ],

    # L8
    [
        ('8', '9', '-'),
        ('8', '10', '-'),
        ('8', '13', '-'),
        ('8', '14', '-'),
        ('8', '15', '-'),
        ('8', '16', '-'),
    ],

    # L9
    [
        ('9', '13', '-'),  # not in original paper
        ('9', '14', '-'),  # not in original paper
        ('9', '15', '-'),
        ('9', '16', '-'),
    ],

    # L10
    [
        ('10', '13', '-'),
        ('10', '14', '-'),
        ('10', '15', '-'),
        ('10', '16', '-'),
    ],

    # L11
    [
        ('11', '13', '-'),
        ('11', '14', '-'),
        ('11', '15', '-'),
        ('11', '16', '-'),
    ],

    # L12
    # none left

    # L13
    [
        ('13', '14', '-'),
        ('13', '15', '-'),
        ('13', '16', '-'),
    ],

    # L14
    [
        ('14', '14', '-'),  # not in original paper
        ('14', '15', '-'),
        ('14', '16', '-'),
    ],

    # L15
    [
        ('15', '16', '-'),
    ],

    # L16
    [
        ('16', '16', '-'),  # not in original paper
    ], )
import copy
smartsGps = copy.deepcopy(reactionDefs)
for gp in smartsGps:
    for j, defn in enumerate(gp):
        # reaction definitions
        g1, g2, bnd = defn
        # get environments in the reaction
        r1 = environs['L' + g1]
        r2 = environs['L' + g2]
        # remove letters?
        g1 = re.sub('[a-z,A-Z]', '', g1)
        g2 = re.sub('[a-z,A-Z]', '', g2)
        # convert reaction into SMARTS string
        sma = '[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]' % (r1, bnd, r2, g1, g2)
        gp[j] = sma

# check reactions
for gp in smartsGps:
    for defn in gp:
        try:
            t = Reactions.ReactionFromSmarts(defn)
            t.Initialize()
        except Exception:
            print(defn)
            raise

# convert environments to molecules
environMatchers = {}
for env, sma in iteritems(environs):
    environMatchers[env] = Chem.MolFromSmarts(sma)

# add pattern to reaction definitions
bondMatchers = []
for i, compats in enumerate(reactionDefs):
    tmp = []
    for i1, i2, bType in compats:
        e1 = environs['L%s' % i1]
        e2 = environs['L%s' % i2]
        patt = '[$(%s)]%s;!@[$(%s)]' % (e1, bType, e2)
        patt = Chem.MolFromSmarts(patt)
        tmp.append((i1, i2, bType, patt))
    bondMatchers.append(tmp)

# convert reactions to rdkit reactions and add reverse reactions
reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])
reverseReactions = []
for i, rxnSet in enumerate(smartsGps):
    for j, sma in enumerate(rxnSet):
        rs, ps = sma.split('>>')
        sma = '%s>>%s' % (ps, rs)
        rxn = Reactions.ReactionFromSmarts(sma)
        labels = re.findall(r'\[([0-9]+?)\*\]', ps)
        rxn._matchers = [Chem.MolFromSmiles('[%s*]' % x) for x in labels]
        reverseReactions.append(rxn)


# Paulas custom BRICS function
# should return fragments as molecules and bonds as atom numbers
# minFragmentSize = minimum number of heavy atoms in a fragment
def myBRICS(mol, minFragmentSize=1, silent=True, returnMols=False):

    global reactions
    # what is the 1 ?
    mSmi = Chem.MolToSmiles(mol, 1)

    allNodes = set()
    allBonds = set()

    activePool = {mSmi: mol}
    allNodes.add(mSmi)
    foundMols = {mSmi: mol}
    # iterate over environments (several reactions per environment)
    for gpIdx, reactionGp in enumerate(reactions):
        newPool = {}
        # iterate over all fragments ? (to break them further apart)
        while activePool:
            matched = False
            nSmi = next(iterkeys(activePool))
            mol = activePool.pop(nSmi)
            # iterate over reactions for this environment
            for rxnIdx, reaction in enumerate(reactionGp):
                if not silent:
                    print('--------')
                    print(smartsGps[gpIdx][rxnIdx])
                # return reaction products (break molecule)
                ps = reaction.RunReactants((mol, ))
                # if a BRICS bond was present and the molecule was broken
                if ps:
                    if not silent:
                        print(nSmi, '->', len(ps), 'products')
                    # iterates over product sequences
                    for prodSeq in ps:
                        print([Chem.MolToSmiles(frag) for frag in prodSeq])
                        seqOk = True
                        # we want to disqualify small fragments, so sort the product sequence by size
                        tSeq = [(prod.GetNumHeavyAtoms(), idx) for idx, prod in enumerate(prodSeq)]
                        tSeq.sort()
                        # iterates over fragments from this reaction
                        for nAts, idx in tSeq:
                            prod = prodSeq[idx]
                            try:
                                # sanitize fragment
                                Chem.SanitizeMol(prod)
                            except Exception:
                                continue
                            # SMILES of fragment
                            pSmi = Chem.MolToSmiles(prod, 1)
                            # check length of fragment
                            if nAts < minFragmentSize:
                                seqOk = False
                                break
                            prod.pSmi = pSmi

                        # [(nAtoms, fragment)] for each fragment from this reaction
                        ts = [(x, prodSeq[y]) for x, y in tSeq]
                        prodSeq = ts
                        if seqOk:
                            matched = True
                            for nAts, prod in prodSeq:
                                pSmi = prod.pSmi
                                # print('\t',nats,pSmi)
                                if pSmi not in allNodes:
                                    activePool[pSmi] = prod
                                    allNodes.add(pSmi)
                                    foundMols[pSmi] = prod
            if not matched:
                newPool[nSmi] = mol
        activePool = newPool
    if not returnMols:
        res = set(activePool.keys())
    else:
        res = activePool.values()
    return res


m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
myBRICS(m, minFragmentSize=2)


def FindBRICSBonds(mol, randomizeOrder=False, silent=True):
    """ returns the bonds in a molecule that BRICS would cleave

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(FindBRICSBonds(m))
    >>> res
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4'))]

    a more complicated case:
    >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
    >>> res = list(FindBRICSBonds(m))
    >>> res
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

    we can also randomize the order of the results:
    >>> random.seed(23)
    >>> res = list(FindBRICSBonds(m,randomizeOrder=True))
    >>> sorted(res)
    [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

    Note that this is a generator function :
    >>> res = FindBRICSBonds(m)
    >>> res
    <generator object ...>
    >>> next(res)
    ((3, 2), ('3', '4'))

    >>> m = Chem.MolFromSmiles('CC=CC')
    >>> res = list(FindBRICSBonds(m))
    >>> sorted(res)
    [((1, 2), ('7', '7'))]

    make sure we don't match ring bonds:
    >>> m = Chem.MolFromSmiles('O=C1NCCC1')
    >>> list(FindBRICSBonds(m))
    []

    another nice one, make sure environment 8 doesn't match something connected
    to a ring atom:
    >>> m = Chem.MolFromSmiles('CC1(C)CCCCC1')
    >>> list(FindBRICSBonds(m))
    []

    """
    letter = re.compile('[a-z,A-Z]')
    indices = list(range(len(bondMatchers)))
    bondsDone = set()
    if randomizeOrder:
        random.shuffle(indices, random=random.random)

    envMatches = {}
    # iterate over BRICS environments
    for env, patt in iteritems(environMatchers):
        # check if environment is present in given molecule
        envMatches[env] = mol.HasSubstructMatch(patt)
    for gpIdx in indices:
        if randomizeOrder:
            compats = bondMatchers[gpIdx][:]
            random.shuffle(compats, random=random.random)
        else:
            compats = bondMatchers[gpIdx]
        # iterate over BRICS bonds
        for i1, i2, bType, patt in compats:
            # check if bond is present in given molecule
            if not envMatches['L' + i1] or not envMatches['L' + i2]:
                continue
            # find index of bond if present
            matches = mol.GetSubstructMatches(patt)
            i1 = letter.sub('', i1)
            i2 = letter.sub('', i2)
            for match in matches:
                if match not in bondsDone and (match[1], match[0]) not in bondsDone:
                    bondsDone.add(match)
                    # get atom indices and environment types of broken bond
                    yield (((match[0], match[1]), (i1, i2)))


def BreakBRICSBonds(mol, bonds=None, sanitize=True, silent=True):
    """ breaks the BRICS bonds in a molecule and returns the results

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2=BreakBRICSBonds(m)
    >>> Chem.MolToSmiles(m2,True)
    '[3*]O[3*].[4*]CC.[4*]CCC'

    a more complicated case:
    >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
    >>> m2=BreakBRICSBonds(m)
    >>> Chem.MolToSmiles(m2,True)
    '[16*]c1ccccc1.[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O'


    can also specify a limited set of bonds to work with:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2 = BreakBRICSBonds(m,[((3, 2), ('3', '4'))])
    >>> Chem.MolToSmiles(m2,True)
    '[3*]OCC.[4*]CCC'

    this can be used as an alternate approach for doing a BRICS decomposition by
    following BreakBRICSBonds with a call to Chem.GetMolFrags:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> m2=BreakBRICSBonds(m)
    >>> frags = Chem.GetMolFrags(m2,asMols=True)
    >>> [Chem.MolToSmiles(x,True) for x in frags]
    ['[4*]CCC', '[3*]O[3*]', '[4*]CC']

    """
    if not bonds:
        # bonds = FindBRICSBonds(mol)
        res = Chem.FragmentOnBRICSBonds(mol)
        if sanitize:
            Chem.SanitizeMol(res)
        return res
    eMol = Chem.EditableMol(mol)
    nAts = mol.GetNumAtoms()

    dummyPositions = []
    for indices, dummyTypes in bonds:
        ia, ib = indices
        obond = mol.GetBondBetweenAtoms(ia, ib)
        bondType = obond.GetBondType()
        eMol.RemoveBond(ia, ib)

        da, db = dummyTypes
        atoma = Chem.Atom(0)
        atoma.SetIsotope(int(da))
        atoma.SetNoImplicit(True)
        idxa = nAts
        nAts += 1
        eMol.AddAtom(atoma)
        eMol.AddBond(ia, idxa, bondType)

        atomb = Chem.Atom(0)
        atomb.SetIsotope(int(db))
        atomb.SetNoImplicit(True)
        idxb = nAts
        nAts += 1
        eMol.AddAtom(atomb)
        eMol.AddBond(ib, idxb, bondType)
        if mol.GetNumConformers():
            dummyPositions.append((idxa, ib))
            dummyPositions.append((idxb, ia))
    res = eMol.GetMol()
    if sanitize:
        Chem.SanitizeMol(res)
    if mol.GetNumConformers():
        for conf in mol.GetConformers():
            resConf = res.GetConformer(conf.GetId())
            for ia, pa in dummyPositions:
                resConf.SetAtomPosition(ia, conf.GetAtomPosition(pa))
    return res


def BRICSDecompose(mol, allNodes=None, minFragmentSize=1, onlyUseReactions=None, silent=True,
                   keepNonLeafNodes=False, singlePass=False, returnMols=False):
    """ returns the BRICS decomposition for a molecule

    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(BRICSDecompose(m))
    >>> sorted(res)
    ['[14*]c1ccccn1', '[16*]c1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]']

    >>> res = list(BRICSDecompose(m,returnMols=True))
    >>> res[0]
    <rdkit.Chem.rdchem.Mol object ...>
    >>> smis = [Chem.MolToSmiles(x,True) for x in res]
    >>> sorted(smis)
    ['[14*]c1ccccn1', '[16*]c1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]']

    nexavar, an example from the paper (corrected):
    >>> m = Chem.MolFromSmiles('CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1')
    >>> res = list(BRICSDecompose(m))
    >>> sorted(res)
    ['[1*]C([1*])=O', '[1*]C([6*])=O', '[14*]c1cc([16*])ccn1', '[16*]c1ccc(Cl)c([16*])c1', '[16*]c1ccc([16*])cc1', '[3*]O[3*]', '[5*]NC', '[5*]N[5*]', '[8*]C(F)(F)F']

    it's also possible to keep pieces that haven't been fully decomposed:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[3*]OCCC', '[3*]O[3*]', '[4*]CC', '[4*]CCC']

    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True))
    >>> sorted(res)
    ['CCCOCc1cccc(-c2ccccn2)c1', '[14*]c1ccccn1', '[16*]c1cccc(-c2ccccn2)c1', '[16*]c1cccc(COCCC)c1', '[16*]c1cccc([16*])c1', '[3*]OCCC', '[3*]OC[8*]', '[3*]OCc1cccc(-c2ccccn2)c1', '[3*]OCc1cccc([16*])c1', '[3*]O[3*]', '[4*]CCC', '[4*]C[8*]', '[4*]Cc1cccc(-c2ccccn2)c1', '[4*]Cc1cccc([16*])c1', '[8*]COCCC']

    or to only do a single pass of decomposition:
    >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
    >>> res = list(BRICSDecompose(m,singlePass=True))
    >>> sorted(res)
    ['CCCOCc1cccc(-c2ccccn2)c1', '[14*]c1ccccn1', '[16*]c1cccc(-c2ccccn2)c1', '[16*]c1cccc(COCCC)c1', '[3*]OCCC', '[3*]OCc1cccc(-c2ccccn2)c1', '[4*]CCC', '[4*]Cc1cccc(-c2ccccn2)c1', '[8*]COCCC']

    setting a minimum size for the fragments:
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True,minFragmentSize=2))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[3*]OCCC', '[4*]CC', '[4*]CCC']
    >>> m = Chem.MolFromSmiles('CCCOCC')
    >>> res = list(BRICSDecompose(m,keepNonLeafNodes=True,minFragmentSize=3))
    >>> sorted(res)
    ['CCCOCC', '[3*]OCC', '[4*]CCC']
    >>> res = list(BRICSDecompose(m,minFragmentSize=2))
    >>> sorted(res)
    ['[3*]OCC', '[3*]OCCC', '[4*]CC', '[4*]CCC']


    """
    global reactions
    mSmi = Chem.MolToSmiles(mol, 1)

    if allNodes is None:
        allNodes = set()

    if mSmi in allNodes:
        return set()

    activePool = {mSmi: mol}
    allNodes.add(mSmi)
    foundMols = {mSmi: mol}
    for gpIdx, reactionGp in enumerate(reactions):
        newPool = {}
        while activePool:
            matched = False
            nSmi = next(iterkeys(activePool))
            mol = activePool.pop(nSmi)
            for rxnIdx, reaction in enumerate(reactionGp):
                if onlyUseReactions and (gpIdx, rxnIdx) not in onlyUseReactions:
                    continue
                if not silent:
                    print('--------')
                    print(smartsGps[gpIdx][rxnIdx])
                # return reaction products (break molecule)
                ps = reaction.RunReactants((mol, ))
                # if a BRICS bond was present and the molecule was broken
                if ps:
                    if not silent:
                        print(nSmi, '->', len(ps), 'products')
                    # iterates over product sequences
                    for prodSeq in ps:
                        seqOk = True
                        # we want to disqualify small fragments, so sort the product sequence by size
                        tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx) for idx, prod in enumerate(prodSeq)]
                        tSeq.sort()
                        # iterates over fragments
                        for nats, idx in tSeq:
                            prod = prodSeq[idx]
                            try:
                                # sanitize fragment
                                Chem.SanitizeMol(prod)
                            except Exception:
                                continue
                            # SMILES of fragment
                            pSmi = Chem.MolToSmiles(prod, 1)
                            # check length of fragment
                            if minFragmentSize > 0:
                                nDummies = pSmi.count('*')
                                if nats - nDummies < minFragmentSize:
                                    seqOk = False
                                    break
                            prod.pSmi = pSmi

                        # [nAtoms, fragment]
                        ts = [(x, prodSeq[y]) for x, y in tSeq]
                        prodSeq = ts
                        if seqOk:
                            matched = True
                            for nats, prod in prodSeq:
                                pSmi = prod.pSmi
                                # print('\t',nats,pSmi)
                                if pSmi not in allNodes:
                                    if not singlePass:
                                        activePool[pSmi] = prod
                                    allNodes.add(pSmi)
                                    foundMols[pSmi] = prod
            if singlePass or keepNonLeafNodes or not matched:
                newPool[nSmi] = mol
        activePool = newPool
    if not (singlePass or keepNonLeafNodes):
        if not returnMols:
            res = set(activePool.keys())
        else:
            res = activePool.values()
    else:
        if not returnMols:
            res = allNodes
        else:
            res = foundMols.values()
    return res
