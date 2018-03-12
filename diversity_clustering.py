'''
Select a diverse subset of compounds from a library
'''
# Import
import sys, os, argparse
import pandas as pd
import numpy as np
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import PandasTools
from rdkit.SimDivFilters import rdSimDivPickers

# compute fingerprints
def getFingerprints(ftype, radius, fps):
    if ftype == 'MACCS':
        return [Chem.rdMolDescriptors.GetMACCSKeysFingerprint(x) for x in fps]
    elif ftype == 'circular':
        return [Chem.rdMolDescriptors.GetMorganFingerprint(x, radius) for x in fps]

def pickSubset(algorithm, fps, size, link=None, seed=None):
    if algorithm == 'MaxMin':
        picker = rdSimDivPickers.MaxMinPicker()
        return picker.LazyPick(distij, len(fps), size, seed=seed)
    elif algorithm == 'hierarchical':
        # compute distance matrix (lower triangle) as 1D array
        dists = []
        for i in range(1,len(fps)):
            sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
            dists.extend([1-x for x in sims])
        dists=np.array(dists)
        if link == 'single':
            linkage = rdSimDivPickers.SLINK
        elif link == 'complete':
            linkage = rdSimDivPickers.CLINK
        elif link == 'centroid':
            linkage = rdSimDivPickers.CENTROID
        elif link == 'average':
            linkage = rdSimDivPickers.UPGMA
        elif link == 'ward':
            linkage = rdSimDivPickers.WARD
        picker = rdSimDivPickers.HierarchicalClusterPicker(linkage)
        return picker.Pick(dists, len(fps), size)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    description='Select a diverse subset of compounds from a library',
    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", metavar="filename", type=str, required=True,
    help='Input file containing SMILES')
    parser.add_argument("-p", "--position", metavar="number", type=int, required=False,
    help='Position of the SMILES column (from 1 to ncolumns). Default: 1', default=1)
    parser.add_argument("--no_header", action='store_true', required=False,
    help='Absence of a header in the input file', default=False)
    parser.add_argument("-o", "--output", metavar="filename", type=str, required=True,
    help='Output file containing the diverse subset')
    parser.add_argument("-s", "--size", metavar="number", type=int, required=True,
    help='Size of the diverse subset (number of molecules)')
    parser.add_argument("-a", "--algorithm", metavar="{MaxMin, hierachical}", type=str, required=False,
    help='Algorithm used for clustering. Default: MaxMin',
    choices=['MaxMin', 'hierarchical'], default='MaxMin')
    parser.add_argument("-l", "--link", metavar="{single, complete, average, centroid, ward}", type=str, required=False,
    help='Linkage method (only used for the hierarchical clustering). Default: complete',
    choices=['single', 'complete', 'average', 'centroid', 'ward'], default='complete')
    parser.add_argument("-f", "--fingerprint", metavar="{MACCS, circular}", type=str, required=False,
    help='Type of fingerprints used for clustering. Default: MACCS',
    choices=['MACCS', 'circular'], default='MACCS')
    parser.add_argument("-r", "--radius", metavar="number", type=int, required=False,
    help='Radius for Morgan circular fingerprints', default=None)
    parser.add_argument("--seed", metavar='number', type=int, required=False,
    help='Random seed (only used for the MaxMin algorithm)', default=-1)
    # Parse arguments from command line
    args = parser.parse_args()
    extension = os.path.splitext(args.input)[1]

    # read input file
    if extension == '.xlsx':
        db = pd.read_excel(args.input, header=None if args.no_header else 0)
    else:
        db = pd.read_csv(args.input, header=None if args.no_header else 0)

    # create molecule by reading SMILES
    smilesCol = db.columns[args.position-1]
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    PandasTools.AddMoleculeColumnToFrame(db, smilesCol=smilesCol, molCol='RDmol')
    lg.setLevel(RDLogger.ERROR)

    # compute fingerprints
    fps = getFingerprints(args.fingerprint, args.radius, db['RDmol'])

    # define distance
    def distij(i,j, fps=fps):
        return 1-DataStructs.TanimotoSimilarity(fps[i],fps[j])

    # pick N molecules
    pickIndices = pickSubset(args.algorithm, fps, size=args.size, link=args.link, seed=args.seed)

    if extension == '.xlsx':
        db.iloc[pickIndices].to_excel(args.output, engine='openpyxl',
        columns=db.columns.drop('RDmol'), index=False)
    else:
        db.iloc[pickIndices].to_csv(args.output,
        columns=db.columns.drop('RDmol'), index=False)
