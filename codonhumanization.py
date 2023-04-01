
import os
import random
import numpy as np
import re
import random

from Bio.SeqUtils import CodonUsage as CU
from Bio.SeqUtils import GC
from Bio.Seq import Seq
import matplotlib.pyplot as plt


pt = os.path.abspath(__file__).replace('codonhumanization.py', '')
class CodonHumanizer():
    def __init__(self):
        model_path = pt + './data/xgboost.model'
        dna2vec_map = pt + './data/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'

        # dna2vec
        file_path = dna2vec_map
        f = open()
        lines = f.readlines()
        f.close()
        self.w2v = dict()
        for line in lines[1:]:
            s = line.split()
            key = s[0]
            value = s[1:]
            value = np.array(list(map(float, value)))
            self.w2v[key] = value
        
        # model
        self.model.load_model(model_path)

        # codon table
        self.codon_table = {
        'F':['TTT', 'TTC'],
        'L':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'I':['ATT', 'ATC', 'ATA'],
        'M':['ATG'],
        'V':['GTT', 'GTC', 'GTA', 'GTG'],
        'S':['TCT', 'TCC', 'TCA', 'TCG'],
        'P':['CCT', 'CCC', 'CCA', 'CCG'],
        'T':['ACT', 'ACC', 'ACA', 'ACG'],
        'A':['GCT', 'GCC', 'GCA', 'GCG'],
        'Y':['TAT', 'TAC'],
        '*':['TAA', 'TAG', 'TGA'],
        'H':['CAT', 'CAC'],
        'Q':['CAA', 'CAG'],
        'N':['AAT', 'AAC'],
        'K':['AAA', 'AAG'],
        'D':['GAT', 'GAC'],
        'E':['GAA', 'GAG'],
        'C':['TGT', 'TGC'],
        'W':['TGG'],
        'R':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S':['AGT', 'AGC'],
        'G':['GGT', 'GGC', 'GGA', 'GGG']
        }

    def make_random_seq(self, length=900):
        new_seq = []
        base = ['A', 'T', 'G', 'C']
        for _ in range(length):
            new_seq.append(random.choice(base))
        return ''.join(new_seq)
    
    def silent_mutation(self, seq):
        codons = re.findall('\w\w\w', seq)
        c = random.randint(0, len(codons)-1) # make point mutation on n_th codon
        codon_list = self.codon_table[str(Seq(codons[c]).translate())]
        codons[c] = random.choice(codon_list) # point mutation
        return ''.join(codons)

    def make_kmer(self, seq):
        n = len(seq)
        bags = []
        for i in range(len(seq)):
            if i+3>n: break
            bags.append(seq[i:i+3])
        return bags

    def vectorizer(self, seq):
        seq = seq.strip()
        kmers = self.make_kmer(seq)
        v = np.zeros(100)
        for kmer in kmers:
            v += self.w2v[kmer]
        return v

    def evolution(self, sequence, n=100, k=100, a=10, early_stopper=0.9999):
        seq_init = sequence.strip()
        if len(seq_init)%3 != 0:
            print("[Warning] Input sequence length must be divisible by 3. Some sequences are removed.")
            print(f'Original sequence length: {len(seq_init)}')
            print(f'Trimmed sequence length: {len(seq_init)-len(seq_init)%3}')
            seq_init = seq_init[:len(seq_init)-len(seq_init)%3]
        
        print(f"[Parameters] n={n}, k={k}, a={a}, early_stopper={early_stopper}")
        probs = []
        pool = []
        pool.append(list())
        for _ in range(k):
            pool[-1].append(seq_init)
        probs.append(self.model.predict_proba([self.vectorizer(seq_init)])[0][1])
        print("[Generation 0] ", probs[0])

        for r in range(n):
            new_pool = list()
            new_pool_probs = []
            print(f'[Generation {r+1}]')
            for seq in pool[-1]:
                original = self.model.predict_proba([self.vectorizer(seq)])[0][1]
                for _ in range(a):
                    new_seq = self.silent_mutation(seq)
                    improved = self.model.predict_proba([self.vectorizer(new_seq)])[0][1]
                    if original < improved:
                        new_pool.append(new_seq)
                        new_pool_probs.append(improved)
            result_list = [i for _,i in sorted(zip(new_pool_probs,new_pool), reverse=True)]
            probs.append(sum(new_pool_probs)/len(new_pool_probs))
            print(f'result: {probs[-1]}')
            pool.append(new_pool[:k])
            if probs[-1] > early_stopper:
                print("Early Stopping..")
                break
        
        # Calculation average CAI for each generation
        codon_index = []
        for p in pool:
            indexs = []
            for seq in p:
                myIndex = CU.CodonAdaptationIndex()
                v = myIndex.cai_for_gene(seq)
                indexs.append(v)
            codon_index.append(sum(indexs)/len(indexs))
        
        plt.plot(range(len(probs)), probs, label='Humanity')
        plt.plot(range(len(codon_index)), codon_index, label='CAI')
        plt.xlabel('Generation')
        plt.legend()
        plt.show()

        return new_pool