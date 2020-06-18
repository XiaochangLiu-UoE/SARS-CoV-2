from oddt import fingerprints, interactions, toolkit
import numpy as np
import csv, os


def extractMOL(sdffile, *num):
    if not os.path.isfile(sdffile):
        raise IOError('No such file: "%s"' % sdffile)
    num = num if num else -1
    count = 0
    mols = {}
    firstline = True
    with open(sdffile, 'r', encoding='utf-8') as o:
        mol = []
        for line in o:
            mol.append(line)
            if firstline:
                name = line
                firstline = False
            if 'M  CHG' in line:
                mol.append('M  END\n')
            if line == '$$$$\n':
                firstline = True
                mols[name] = mol
                mol = []
                count += 1
                if count == num:
                    break
                else:
                    continue
    with open('extracted.sdf', 'w+', encoding='utf-8') as w:
        for molecule in mols.values():
            for line in molecule:
                w.write(line)


def detectInteractions(protein, ligand):
    # interactions to be detected
    interactions_dict = {'hbonds', 'halogenbonds',
                         'pi_stacking', 'salt_bridges',
                         'hydrophobic_contacts', 'pi_cation', 'pi_metal'}
    for contact in interactions_dict:
        exec('%s = set()' % contact)
    # prepare protein and ligand objects using oddt
    suffix = [protein.split('.')[1], ligand.split('.')[1]]
    protein = next(toolkit.readfile(suffix[0], protein))
    protein.protein = True
    ligand = list(toolkit.readfile(suffix[1], ligand))
    # to store ligands' interactions with protein
    profiles = []
    # start to detect, ligand by ligand
    for lig in ligand:
        # to store each ligand's interactions with protein
        profile = {}
        """
        start to detect interactions
        for each ligand, detect interaction at atom level, each contacted residue's atom will be recorded
        for the binding site, the contact will stay at residue level
        """
        for contact in interactions_dict:
            # this function will return contacting atoms
            # and a boolean array indicating whether the interaction is strict
            exec('contact_atoms = interactions.%s(protein, lig)' % contact)
            if contact in ['honds', 'halogenbonds']:
                profile[contact] = locals()['contact_atoms'][0]['resnum'][locals()['contact_atoms'][2]]
            elif contact == 'pi_stacking':
                # face to face; edge to face
                profile[contact] = np.concatenate((locals()['contact_atoms'][0]['resnum'][locals()['contact_atoms'][2]],
                                                   locals()['contact_atoms'][0]['resnum'][locals()['contact_atoms'][3]])
                                                  , axis=None)
            elif contact in ['salt_bridges', 'hydrophobic_contacts']:
                profile[contact] = locals()['contact_atoms'][0]['resnum']
            else:
                # ring; cation or metal
                profile[contact] = np.concatenate((locals()['contact_atoms'][0]['resnum'][locals()['contact_atoms'][2]],
                                                   locals()['contact_atoms'][1]['resnum'][locals()['contact_atoms'][2]])
                                                  , axis=None)
            locals()[contact].update(profile[contact]) if profile[contact].size != 0 else exec('pass')
        profiles.append(profile)
    bindingsite = {}
    for contact in interactions_dict:
        bindingsite[contact] = locals()[contact]
    return bindingsite, profiles


def makeFP(bindingsite, ligand):
    FP = np.array([])
    for contact in bindingsite:
        contact_array = np.array(list(bindingsite[contact]))
        contact_FP = np.zeros(contact_array.size)
        for resi in ligand[contact]:
            for index, residue in enumerate(contact_array):
                if resi == residue:
                    contact_FP[index] += 1
        FP = np.concatenate((FP, contact_FP), axis=None)
    return FP


def calSimMatrix(bindingsite, ligand, similarity='tanimoto', heading=None):

    def calEUDist(a, b):
        import math
        dist = 0
        for i in range(len(a)):
            dist += (a[i] - b[i]) ** 2
        return math.sqrt(dist)

    SimiMatrix = []
    for i in range(len(ligand)):
        if heading:
            each_result = [heading[i]]
        else:
            each_result = [i + 1]
        ref = makeFP(bindingsite=bindingsite, ligand=ligand[i])
        for j in range(len(ligand)):
            exec('score = 0')
            query = makeFP(bindingsite, ligand[j])
            if similarity == 'euclidean distance':
                locals()['score'] = calEUDist(a=ref, b=query)
            else:
                #score = eval('fingerprints.%s(a=ref, b=query)' % similarity)
                exec('locals()["score"] = fingerprints.%s(a=ref, b=query)' % similarity)
            each_result.append(locals()['score'])
        SimiMatrix.append(each_result)
    return SimiMatrix


def makeHeatmap(matrix, filename, heading=None):
    x_axis = [' ']
    if heading:
        x_axis.extend(heading)
    else:
        x_axis.extend([i for i in range(1, len(matrix) + 1)])
    with open(filename, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(x_axis)
        for row in matrix:
            writer.writerow(row)


def clusterMOL(bindingsite, ligand, sdffile, similarity='tanimoto'):
    output = {}
    mols = []
    with open(sdffile, 'r', encoding='utf-8') as o:
        temp = []
        for line in o:
            temp.append(line)
            if line == '$$$$\n':
                mols.append(temp)
                temp = []
    for i in range(len(ligand)):
        output[i] = []
        ref = makeFP(bindingsite=bindingsite, ligand=ligand[i])
        for j in range(len(ligand)):
            query = makeFP(bindingsite=bindingsite, ligand=ligand[j])
            exec('score = fingerprints.%s(a=ref, b=query)' % similarity)
            if locals()['score'] > 0.5:
                output[i].append(j)
    for i in output:
        with open('Cluster/'+str(i+1)+'.sdf', 'w+', encoding='utf-8') as w:
            for j in output[i]:
                for line in mols[j]:
                    w.write(line)
        print(str(i+1) + ' done')



if __name__ == '__main__':
    #extractMOL('myresults.sdf')
    bindingsite, ligands = detectInteractions(protein='protein_fixed.pdb', ligand='extracted.sdf')
    simimatrix = calSimMatrix(bindingsite=bindingsite, ligand=ligands)
    makeHeatmap(matrix=simimatrix, filename='heatmap_improved_2.csv')
    #clusterMOL(bindingsite=bindingsite, ligand=ligands, sdffile='extracted.sdf')
    # use fingerpints built in oddt to conduct the same experiment
    """protein = next(toolkit.readfile('pdb', 'protein_fixed.pdb'))
    protein.protein = True
    ligands = list(toolkit.readfile('sdf', 'extracted.sdf'))
    rows = []
    for i in range(len(ligands)):
        row = [i+1]
        ref = fingerprints.SimpleInteractionFingerprint(ligand=ligands[i], protein=protein)
        for j in range(len(ligands)):
            query = fingerprints.SimpleInteractionFingerprint(ligand=ligands[j], protein=protein)
            score = fingerprints.tanimoto(a=ref, b=query)
            #score = fingerprints.similarity_SPLIF(reference=ref, query=query)
            row.append(score)
        rows.append(row)
    axis = [' ']
    axis.extend([i for i in range(1, 91)])
    with open('SIF_heatmap.csv', 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(axis)
        for row in rows:
            writer.writerow(row)"""
