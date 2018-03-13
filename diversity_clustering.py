'''
Select a diverse subset of compounds from a library
'''
# Import
import sys, os, argparse
import pandas as pd
import numpy as np
from concurrent import futures
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import PandasTools
from rdkit.SimDivFilters import rdSimDivPickers

def getMACCS_fps(mol):
    return Chem.rdMolDescriptors.GetMACCSKeysFingerprint(mol)

def getCircular_fps(mol, radius):
    return Chem.rdMolDescriptors.GetMorganFingerprint(mol, radius)

# Compute fingerprints
def getFingerprints(mols, args):
    sys.stdout.write('Computing {} fingerprints\n'.format(args.fingerprint))
    # uses a pool of processes to execute calls asynchronously
    with futures.ProcessPoolExecutor(max_workers=args.cpu) as executor:
        if args.fingerprint == 'MACCS':
            fps = executor.map(getMACCS_fps, mols)
        elif args.fingerprint == 'circular':
            fps = executor.map(getCircular_fps, mols, args.radius)
    fps = [fp for fp in fps]
    return fps

# Compute distance matrix as 1D array
def computeDistanceMatrix(fps, args):
    sys.stdout.write('Computing distance matrix\n')
    distmatrix = []
    # uses a pool of processes to execute calls asynchronously
    with futures.ProcessPoolExecutor(max_workers=args.cpu) as executor:
        # get results in order of submission using map
        for distances in executor.map(getDistances, range(1,len(fps))):
            distmatrix.extend(distances)
    return np.array(distmatrix)

# Pick a subset of molecules from a distance matrix
def pickSubset(distmatrix, nfps, args):
    sys.stdout.write('Picking a subset of molecules using {} clustering algorithm\n'.format(args.algorithm))
    if args.algorithm == 'MaxMin':
        picker = rdSimDivPickers.MaxMinPicker()
        return picker.Pick(distmatrix, nfps, args.size, seed=args.seed)
    elif args.algorithm == 'hierarchical':
        if args.link == 'single':
            linkage = rdSimDivPickers.SLINK
        elif args.link == 'complete':
            linkage = rdSimDivPickers.CLINK
        elif args.link == 'centroid':
            linkage = rdSimDivPickers.CENTROID
        elif args.link == 'average':
            linkage = rdSimDivPickers.UPGMA
        elif args.link == 'ward':
            linkage = rdSimDivPickers.WARD
        picker = rdSimDivPickers.HierarchicalClusterPicker(linkage)
        return picker.Pick(distmatrix, nfps, args.size)


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
    parser.add_argument("--cpu", metavar='number', type=int, required=False,
    help='Number of CPU cores to use', default=None)
    # Parse arguments from command line
    args = parser.parse_args()
    extension_input = os.path.splitext(args.input)[1]
    extension_output = os.path.splitext(args.output)[1]

    # read input file
    sys.stdout.write('Reading input file {}\n'.format(args.input))
    if extension_input in ['.xlsx', '.xls']:
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
    fps = getFingerprints(db['RDmol'], args)
    nfps = len(fps)
    db.drop(columns=['RDmol'], inplace=True)

    # Compute Soergel distances for a given compound
    def getDistances(i):
        return [1-x for x in DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])]

    # compute distance matrix
    distmatrix = computeDistanceMatrix(fps, args)

    # pick N molecules
    pickIndices = pickSubset(distmatrix, nfps, args)

    if extension_output in ['.xlsx', '.xls']:
        db.iloc[pickIndices].to_excel(args.output, engine='openpyxl', index=False)
    else:
        db.iloc[pickIndices].to_csv(args.output, index=False)
