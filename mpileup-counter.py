# Based on https://github.com/Niknafs/NGSTools/blob/master/baseParser.py.

import os
import sys
import argparse
import time
from multiprocessing import Pool
import tqdm


# Args
parser = argparse.ArgumentParser(description='Quantify alt alleles.')
parser.add_argument('--mpileup', type=str, help='Mpileup File Input')
parser.add_argument('--output', type=str, help='Counts File Output')
args = parser.parse_args()


# Parser Class
class parseString(object):
    
    def __init__(self, ref, string):
        self.ref = ref.upper()
        self.string = string.upper()
        self.types = {'A':0,'G':0,'C':0,'T':0,'-':[],'*':0,'+':[],'X':[]}
        self.process()
        
    def process(self):
        # remove end of read character
        self.string = self.string.replace('$','')
        while self.string != '':
            if self.string[0] == '^':
                # skip two characters when encountering '^' as it indicates
                # a read start mark and the read mapping quality
                self.string = self.string[2:]
            elif self.string[0] == '*':
                self.types['*'] += 1
                # skip to next character
                self.string = self.string[1:]
            
            elif self.string[0] in ['.',',']:
                if (len(self.string)== 1) or (self.string[1] not in ['+','-']):
                    # a reference base
                    self.types[self.ref] += 1
                    self.string = self.string[1:]
                elif self.string[1] == '+': 
                    insertionLength = int(self.string[2])
                    insertionSeq = self.string[3:3+ insertionLength]
                    self.types['+'].append(insertionSeq)
                    self.string = self.string[3+insertionLength:]
                elif self.string[1] == '-':
                    deletionLength = int(self.string[2])
                    deletionSeq = self.string[3:3+deletionLength]
                    self.types['-'].append(deletionSeq)
                    self.string = self.string[3+deletionLength:]
                    
            elif self.string[0] in self.types and\
                 ((len(self.string)==1) or (self.string[1] not in ['-','+'])):
                # one of the four bases
                self.types[self.string[0]] += 1
                self.string = self.string[1:]
            else:
                # unrecognized character
                # or a read that reports a substitition followed by an insertion/deletion
                self.types['X'].append(self.string[0])
                self.string = self.string[1:]
        return
    def __repr__(self):
        types = self.types
        return '\t'.join(
            list(map(str,[types['A'], types['C'], types['G'],types['T'],types['*']])) + \
            list(map(','.join, [types['-'],types['+'],types['X']]))
            )
        
# Single function for multiprocessing
def process_line(line):
    toks = line.strip('\n').split('\t')
    ref = toks[2].upper()
    cov = toks[3]
    out_str = '\t'.join([toks[0], toks[1],ref, cov]) + '\t' + \
        parseString(ref, toks[4]).__repr__()

    return out_str

def main():
    with open(args.mpileup, 'r') as f:
        with Pool() as pool:
            out_rows = pool.map(
                process_line, 
                f.readlines()
                )

    with open(args.output, 'w+') as f:
        f.write('\t'.join([
            "chrom", "pos", "ref", "cov", "A", 
            "C", "G", "T", "*", "-", "+", "X"
            ]))

        f.write('\n'.join(out_rows))
        f.write('\n')

if __name__ == '__main__':
    start = time.time()
    main()
    time_taken = time.time() - start

    print(f"Mpileup Counter Done! Time: {time_taken:.2f}s")
