import os
from openbabel import pybel
from InteractionFP.classifyMols import detectInteractions, calSimMatrix, makeHeatmap


def countSpace(string):
    space = []
    nspace = 0
    for i, s in enumerate(string):
        if s.isspace():
            nspace += 1
        else:
            if nspace == 0:
                continue
            space.append(nspace)
            nspace = 0
    return space


def extractLigand(pdbfile):
    protein = []
    ligand = []
    with open(pdbfile, 'r', encoding='utf-8') as o:
        for line in o:
            words = line.split()
            if 'LIG' not in words:
                protein.append(line)
            elif 'HETATM' in words and 'LIG' in words:
                ligand.append(line)
        ligand = list(sorted(ligand, key=lambda x: x.split()[2][1:]))
    return protein, ligand


def makeFile(lines, filename):
    with open(filename, 'w+', encoding='utf-8') as w:
        for line in lines:
            w.write(line)


if __name__ == '__main__':
    path = 'C:\\*\\*\\COVID-19\\Mpro_All_PDBs - ver 2020-03-24\\'
    pdbfiles = os.listdir(path)
    heading = []
    proteins = []
    ligands = []
    for pdbfile in pdbfiles:
        heading.append(pdbfile.strip('.pdb'))
        protein, ligand = extractLigand(pdbfile=path+pdbfile)
        if pdbfiles.index(pdbfile) == -1:
            ligand.append('ENDMDL\n')
            ligand.append('END\n')
        else:
            ligand.append('ENDMDL\n')
        ligands.extend(ligand)
        proteins.append(protein)
    makeFile(proteins[0], 'protein.pdb')
    makeFile(ligands, 'temp.pdb')
    ligand_sdf = pybel.Outputfile('sdf', 'ligands.sdf', overwrite=True)
    for mol in pybel.readfile('pdb', 'temp.pdb'):
        ligand_sdf.write(mol)
    ligand_sdf.close()
    bindingsite, ligands = detectInteractions(protein='protein.pdb', ligand='ligands.sdf')
    simimatrix = calSimMatrix(bindingsite=bindingsite, ligand=ligands, similarity='dice')
    makeHeatmap(matrix=simimatrix, filename='sar-cov-2_heatmap_dice.csv')
