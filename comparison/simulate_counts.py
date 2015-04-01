import numpy as np
import numpy.random as nr
from rpy2.robjects import *


class Counts:

    def __init__(self, counts):
        """set first paralog counts"""
        self.tot1 = np.array(counts[0], dtype = 'f')
        self.unambig1 = np.array([max(x, 0.1) for x in counts[1]], dtype ='f')
        self.ambig1 = np.array(counts[2], dtype = 'f')

    def add_paralog(self, counts):
        """set second paralog counts, computes distributions"""
        self.tot2 = np.array(counts[0], dtype = 'f')
        self.unambig2 = np.array(counts[1], dtype = 'f')
        self.ambig2 = np.array(counts[2], dtype = 'f')
        self.compute_distr()
        self.compute_distr_reps()

    def compute_distr(self):
        """assignes ambiguous counts to paralogs, based on unambiguous ratio"""
        tot_ambig = self.ambig1 + self.ambig2
        tot_unambig = self.unambig1 + self.unambig2
        unambig1_ratio = self.unambig1 / tot_unambig
        ambig_to1 = unambig1_ratio * tot_ambig
        ambig_to2 = tot_ambig - ambig_to1
        self.distr1 = self.unambig1 + ambig_to1
        self.distr2 = self.unambig2 + ambig_to2

    def compute_distr_reps(self):
        """the unambiguous ratio is from reps average, this time"""
        tot_ambig = self.ambig1 + self.ambig2
        tot_unambig = self.unambig1 + self.unambig2
        unambig1_ratio = np.mean(self.unambig1 / tot_unambig)
        ambig_to1 = unambig1_ratio * tot_ambig
        ambig_to2 = tot_ambig - ambig_to1
        self.distr_rep1 = self.unambig1 + ambig_to1
        self.distr_rep2 = self.unambig2 + ambig_to2


class SimParalog:

    def __init__(self):
        """default options: reads are 100bp long, 3 replicates"""
        self.seed = nr.seed()
        self.ambig_range = []
        self.unambig_range = []
        self.len_reads = 100
        self.mmatches = []
        self.len_cDNA = 0
        self.num_reps = 3

    def set_paralogs(self, len_cDNA, mmatches):
        """set length of cDNA and its mismatches."""
        self.len_cDNA = len_cDNA
        self.mmatches = mmatches

    def set_ranges(self):
        """set ambiguous and unambiguous ranges."""
        for mmatch in self.mmatches:
            start = max(1, mmatch - self.len_reads + 1)
            for k in range(start, mmatch + 1):
                self.unambig_range.append(k)
        for k in range(1, self.len_cDNA + 1):
            if k not in self.unambig_range:
                self.ambig_range.append(k)
        
    def sim_smn(self):
        """set SMN1-SMN2 cDNA length and mismatches."""
        self.set_paralogs(1536, [814, 1098])
        self.set_ranges()

    def get_ambig_bases_ratio(self):
        """get fraction of ambiguous bases."""
        n_ambig = float(len(self.ambig_range))
        return n_ambig / (len(self.unambig_range) + n_ambig)

    def get_unambig_reads(self, rand_poss):
        """return number of unambiguous reads, given rand_poss array."""
        num_unambig_reads = 0
        for rand_pos in rand_poss:
            if rand_pos in self.unambig_range:
                num_unambig_reads += 1
            elif rand_pos not in self.ambig_range:
                assert False, 'error! ranges are not correct'
        return num_unambig_reads

    def sim_unambig_counts(self, num_sims):
        """simulate num_sims read positions and return num of unambiguous."""
        rand_poss = nr.randint(1, self.len_cDNA + 1, num_sims)
        return self.get_unambig_reads(rand_poss)

    def sim_nbinom(self, mean):
        """return rand_nbinom, given mean.
        dispersion / mean trend from UPF1_KD_rescue."""
        disp = 1 / (0.008968 + 2.070553 / mean)
        rand_nbinom = r.rnbinom(1, mu = mean, size = disp)[0]
        return rand_nbinom

    def sim_counts(self, mean):
        """compute gene counts, given a mean and setup."""
        tot_counts = [self.sim_nbinom(mean) for k in range(self.num_reps)]
        unambig_counts = [self.sim_unambig_counts(counts) for counts in tot_counts]
        ambig_counts = [tot_counts[k] - unambig_counts[k] for k in range(self.num_reps)]
        return tot_counts, unambig_counts, ambig_counts

    def run_sim(self, mean1, mean2):
        """run a simulation of 2 paralogs with mean1 and mean2."""
        counts = Counts(self.sim_counts(mean1))
        counts.add_paralog(self.sim_counts(mean2))
        return counts

x = SimParalog()
x.sim_smn()
true_fc = 2
unambig_err = 0
distr_err = 0
distr_rep_err = 0
countsA = x.run_sim(1000,2000)
countsB = x.run_sim(2000,2000)
unambig_fc = np.mean(countsB.unambig1 / countsA.unambig1)
distr_fc = np.mean(countsB.distr1 / countsA.distr1)
distr_rep_fc = np.mean(countsB.distr_rep1 / countsA.distr_rep1)
print 'unambig1 ' + ' '.join(str(int(round(x))) for x in countsA.unambig1.tolist() + countsB.unambig1.tolist())
print 'unambig2 ' + ' '.join(str(int(round(x))) for x in countsA.unambig2.tolist() + countsB.unambig2.tolist())
print 'distr1 ' + ' '.join(str(int(round(x))) for x in countsA.distr1.tolist() + countsB.distr1.tolist())
print 'distr2 ' + ' '.join(str(int(round(x))) for x in countsA.distr2.tolist() + countsB.distr2.tolist())
