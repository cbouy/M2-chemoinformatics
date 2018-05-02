import re, argparse

def get_proteins(*args, residues=None):
    """Create proteins from alignment input file"""
    proteins = []

    for input_file in args:
        # read alignment file
        with open(input_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            # parse file as expected by modeller
            header = re.search(r'^>P1;(.+)$', line)
            info   = re.search(r'^(.*:.*:.*:.*:.*:.*:.*:.*:.*:.*)$', line)
            end    = re.search(r'^.*\*$', line)

            if header:
                protein = {'sequence':''}
                protein['code'] = header.group(1)
            elif info:
                data = info.group(1).split(':')
                protein['method']       = data[0]
                protein['PDB']          = data[1]
                protein['begin']        = data[2]
                protein['begin_chain']  = data[3]
                protein['end']          = data[4]
                protein['end_chain']    = data[5]
                protein['name']         = data[6]
                protein['source']       = data[7]
                protein['resolution']   = data[8]
                protein['r-factor']     = data[9]
            else:
                protein['sequence'] += line.replace('\n','')
                if end:
                    protein['sequence'] = protein['sequence'][:-1]
                    proteins.append(protein)

    # if list of residues was passed
    if residues:
        for protein in proteins:
            # restrict the sequence to only contain these residues
            seq = ""
            for r in residues: # r = 1-36 for residue 1 to 36
                begin, end = [int(x) for x in r.split('-')]
                cut = protein['sequence'][begin-1:end]
                seq += cut
            protein['sequence'] = seq

    return proteins


def compare_sequences(target, structure, ignore_gaps=True):
    """Compare a template sequence and a target"""
    assert len(target['sequence']) == len(structure['sequence'])
    match    = 0
    mismatch = 0

    for i in range(len(target['sequence'])):
        if (target['sequence'][i] == '-') and (structure['sequence'][i] == '-'):
            continue
        if (target['sequence'][i] == '-') or (structure['sequence'][i] == '-'):
            if ignore_gaps:
                continue
            else:
                mismatch += 1
                continue
        if target['sequence'][i] == structure['sequence'][i]:
            match += 1
        else:
            mismatch += 1

    return match, mismatch

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute sequence identity between target proteins and templates from a PIR alignment file")
    parser.add_argument('-i', '--input', nargs='+', help='PIR alignment files for modeller', required=True)
    parser.add_argument('-r', '--resid', metavar='begin-end', nargs='*', help='Residues id used to compute the sequence identity')
    parser.add_argument('--ignoregaps', action="store_true", help="Don't count gaps as mismatch", default=False)
    parser.add_argument('-o', '--output', help='Output file')
    args = parser.parse_args()

    proteins = get_proteins(*args.input, residues=args.resid)
    structures = []
    sequences  = []
    results    = []

    for protein in proteins:
        if 'sequence' in protein['method']:
            sequences.append(protein)
        else:
            structures.append(protein)

    print(' '*9,''.join(['{:^9}'.format(structure['code']) for structure in structures]))
    print(' '*9,''.join(['{:―^9}'.format('') for structure in structures]))
    for target in sequences:
        temp = []
        print('{:^9}│'.format(target['code']), end='')
        for structure in structures:
            match, mismatch = compare_sequences(target, structure, ignore_gaps=args.ignoregaps)
            pid = 100*match/(match + mismatch)
            temp.append(pid)
            print('{:^ 9.2f}'.format(pid), end='')
        results.append(temp)
        print()

    if args.output:
        with open(args.output, 'w') as f:
            header = 'Sequence,' + ','.join([structure['code'] for structure in structures]) + '\n'
            f.write(header)
            for target, pids in zip(sequences, results):
                line = target['code'] + ',' + ','.join(['{:.3f}'.format(pid) for pid in pids]) + '\n'
                f.write(line)
