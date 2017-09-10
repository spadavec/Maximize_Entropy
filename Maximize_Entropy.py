import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import csv
import numpy
from rdkit.Chem import DataStructs
from random import randint
from operator import itemgetter
import math
import csv
import argparse



def get_index_totals(fps):
    """ Counts the number of '1'(e.g. 'on') bits at every FP index of a given mol """

    indexes = [x for x in range(1024)]
    idx_counts = {key: 0 for key in indexes}
    ct = 0
    for fp in fps:
        for i,x in enumerate(fp):
            if x is '1':
                ct += 1
                idx_counts[i] += 1

    return idx_counts


def check_addition(train_fps, test_fps):
     """ Checks the change in bit diversity of a given fp (test_fps) against
        a set of known fps (train_fps) """
     sums = 0

     # Get the sum vector [1,0,12,23,....0,0,2]
     train_sums = sum_fps(train_fps)

     for i,x in enumerate(train_sums):
         if x == 0 and test_fps[i] == 1:
             sums += 1

     return sums


def sum_fps(train_fps):
    """ Generates a "master" fps which sums all train_fps at each index """
    master_fps = numpy.zeros(1024)

    for i,x in enumerate(master_fps):
        index_sum = 0
        for fp in train_fps:
            if fp[i] == 1:
                index_sum += 1

        master_fps[i] = index_sum

    return master_fps




def read_csv(filename):
    """ Presumes input csv file is in format of ids/smiles/ic50 """

    smiles = []
    ic50 = []
    ids = []

    with open(filename) as ifile:
        data = csv.reader(ifile, delimiter=',')

        for line in data:
            ids.append(line[0])
            smiles.append(line[1])
            ic50.append(float(line[2]))

    return ids, smiles, ic50


def generate_fingerprints(smiles):
    """ Generates ECFP4 circular fingerprints in approrpate format """
    temp = []
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    bit_info = {}

    fp = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024, bitInfo=bit_info) for x in mols]

    for arr in fp:
            temp.append([int(x) for x in arr])
    fps = numpy.array(temp, dtype=float)

    return fps


def write_csv(train,test, output_name):


    with open('{}-train.csv'.format(output_name), 'w') as ofile:
        writer = csv.writer(ofile)
        writer.writerows(train)
    with open('{}-test.csv'.format(output_name), 'w') as ofile2:
        writer=csv.writer(ofile2)
        writer.writerows(test)



def end_condition(train, test, target_diversity, output_name):
    """ Checks if end condition is met--in this case it waits until a threshold bit diversity is reached. If met, writes output files """

    temp = 0
    train_smiles = [x[1] for x in train]
    train_fps = generate_fingerprints(train_smiles)
    sums = sum_fps(train_fps)

    for i,x in enumerate(sums):
        if x != 0:
            temp += 1

    fraction = float(temp)/1024

    if fraction < target_diversity:
        return False
    else:
        write_csv(train,test, output_name)
        return True

def main():
    test = []
    train = []
    iteration_num = 1

    # Argument parser
    parser = argparse.ArgumentParser(description='argument parser for shannon entropy search')
    parser.add_argument('--input', type=str, dest="input_file")
    parser.add_argument('--output', type=str, dest="output_file")
    parser.add_argument('--cutoff', type=float, dest="diversity_cutoff")

    args = parser.parse_args()
    input_csv_name = args.input_file
    diversity_cutoff = args.diversity_cutoff
    output_file = args.output_file

    print ""
    print "Grabbing data..."

    # Grab data
    ids, smiles, ic50 = read_csv(input_csv_name)

    # Create training set
    # For now, grabs just the last compound in the list
    train.append([ids[-1],smiles[-1],ic50[-1]])

    # Create the test set
    num = len(smiles)
    for i in range(num-1):
        test.append([ids[i],smiles[i], ic50[i]])

    print "Begginning search..."

    while not end_condition(train, test, diversity_cutoff, output_file):
        se_idxs = []
        se = []
        print ""
        print "-"*25
        print "Iteration #{}".format(iteration_num)
        print "-"*25

        # Generate all necessary training data
        train_smiles = [x[1] for x  in train]

        print train_smiles
        train_fps = generate_fingerprints(train_smiles)

        test_smiles = [x[1] for x in test]

        # Generate intial actives
        test_fps = generate_fingerprints(test_smiles)

        additional_bits = []

        # Get all of the entropies of the test molecules
        for x in test_fps:
            additional_bits.append(check_addition(train_fps,x))

        max_index = additional_bits.index(max(additional_bits))

        train.append(test[max_index])
        del test[max_index]

        iteration_num += 1



if __name__ == "__main__":
    main()
