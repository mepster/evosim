from scipy.stats import binom
import matplotlib.pyplot as plt
import numpy as np
from helpers import *
import math
import random
import multiprocessing
import copy
import os

import warnings
warnings.filterwarnings('ignore')

class Environment:
    def __init__(self, args, idx, sA, sB):
        self.args = args
        self.idx = idx
        self.sA = sA
        self.sB = sB
        self.active_env = [] # pairs of (gen, env)

    def __repr__(self):
        ret = f"env {self.label()}:\tsA={self.sA}\tsB={self.sB}"
        return ret

    def label(self):
        return "A" if self.idx == 0 else "B"


    def swap(self, gen): # swap which allele is favored
        self.idx = (self.idx+1)%2
        tmp = self.sA
        self.sA = self.sB
        self.sB = tmp

        # log data for plotting
        self.active_env.append((gen, "A" if self.idx == 0 else "B"))

class Clade:
    # a clade contains two genotypes: A and B, with the same fixed mutation rate
    def __init__(self, args, m):
        # initialize one clade
        self.args = args
        self.N = args['N']
        aToB = args['aToB']
        numClades = args['numClades']
        minMu = args['minMu']

        self.m = m # m is 0, 1, 2, ... for clades labelled M2, M3, M4, ...
        self.mu = 10.0**(-1.0 * (m+minMu)) # mutation rate
        self.label = f"M{m+minMu}" # clade labels M2, M3, M4, ... for printing

        fracA = 1.0 / numClades * aToB  # initial fraction of A genotype
        fracB = 1.0 / numClades * (1.0 - aToB) # initial fraction of B genotype
        if self.N < 0:  # infinite pop - floating point "frequency" of each genotype
            self.nA = fracA
            self.nB = fracB
        else: # finite pop - integer count of each genotype
            self.nA, self.nB = binom.rvs(n=round(self.N), p=[fracA, fracB])

        # for logging of the counts of the 2 genotypes within the clade over time
        self.countsA = []
        self.countsB = []
        self.counts = []

    def update_gen_data(self):
        # log the counts of the 2 genotypes within the clade
        self.countsA.append(self.nA) # A genotype
        self.countsB.append(self.nB) # B genotype
        self.counts.append(self.nA+self.nB) # total

    def __repr__(self):
        # for printing clade info in pretty format
        ret = f"\t\tclade {self.label}:"
        if self.N < 0:
            ret = ret + f"\t(freqA={self.nA:.2e}\t freqB={self.nB:.2e})\t(mu:{self.mu:4.1e})"
        else:
            ret = ret + f"\t(n={self.nA+self.nB:6.0f}\t nA={self.nA:6.0f}\t nB={self.nB:6.0f})\t(mu:{self.mu:4.1e})"
        return ret

    def select(self, env):
        # apply selection, in frequency space
        # nA and nB in each clade are non-integer after this, and sumN != N
        self.nA = self.nA*(1.0+env.sA)
        self.nB = self.nB*(1.0+env.sB)

    def mutate(self):
        # apply mutation, in frequency space
        # nA and nB in each clade are non-integer after this, and sumN != N
        AToB = self.nA*self.mu
        AtoA = self.nA*(1.0-self.mu)
        BtoA = self.nB*self.mu
        BtoB = self.nB*(1.0-self.mu)
        #print(f"AtoA:{AtoA} AToB:{AToB} BtoB:{BtoB} BtoA:{BtoA}")
        self.nA = AtoA + BtoA
        self.nB = AToB + BtoB

    def normalize(self, sumN, N):
        # convert from frequency space back to integer counts
        fracA = float(self.nA) / sumN # fraction of A
        fracB = float(self.nB) / sumN # fraction of B
        if N < 0: # this means it's an infinite pop
            # infinite pop - don't convert to integer
            self.nA = fracA # total N of this type
            self.nB = fracB # total N of this type
        else:
            # finite pop - convert to integer
            # save computing binom if frac is 0.0 anyway, just set count to 0
            self.nA = 0 if fracA == 0.0 else binom.rvs(n=round(N), p=fracA) # total N of this type
            self.nB = 0 if fracB == 0.0 else binom.rvs(n=round(N), p=fracB) # total N of this type


class Pop:
    # a population contains multiple clades
    def __init__(self, args, idx):
        # initialize one population
        self.args = args
        self.N = args['N']
        self.env = args['env']
        numClades = args['numClades']
        self.label = f"{idx}"

        # the list of clades in this population
        self.clades = []
        for m in range(0, numClades):
            clade = Clade(args, m=m)
            self.clades.append(clade)
        # args["numClades"] = 1
        # clade = Clade(numClades-3, args)
        # self.clades.append(clade)

    def update_gen_data(self):
        # for each clade, log its data for this generation
        for clade in self.clades:
            clade.update_gen_data()
        # update the Pop's own data here, if any

    def __repr__(self):
        # for printing population info in pretty format
        ret = f"\tpop {self.label}\tsumN: {self.sumN():.2f}\t({self.env})\n"
        # this includes printing each clade
        for idx, clade in enumerate(self.clades):
            ret = ret + f"{clade}"
            
            if idx < len(self.clades)-1:
                ret = ret + "\n"
        return ret

    def sumN(self): # counts individuals in this population
        sumN = 0.0
        for clade in self.clades:
            sumN = sumN + (clade.nA + clade.nB)
        assert sumN > 0.0, f"uh oh, sumN for a population should be > 0. Probably N {self.N} is too small."
        return sumN

    def select(self):
        # nA and nB in each clade are frequencies after this, and sum != N
        for clade in self.clades:
            clade.select(self.env)

    def mutate(self):
        # nA and nB in each clade are frequencies after this, and sum != N
        for clade in self.clades:
            clade.mutate()

    def normalize(self):
        # convert from frequency space back to integer counts
        sumN = self.sumN()
        for clade in self.clades:
            clade.normalize(sumN, self.N)

    def one_gen_mutation_mode(self):
        self.update_gen_data()
        self.select() # converts from int to freq
        if self.args.includeMutation:
            self.mutate() # works on freq
        self.normalize() # converts from freq to int


class Grid():
    # A grid contains multiple populations
    # The point is that you can have migration between populations
    def __init__(self, args):
        # initialize the grid of populations (there is only one grid per run)
        self.args = args
        self.numPops = args['numPops']
        self.numClades = args['numClades']
        self.mode = args['mode'] # 'mutation' or 'migration
        self.printStats = args['printStats']

        s = args['s']
        aToB = args['aToB']

        self.pops = []
        for idx in range(self.numPops):
            if idx%2 == 0:
                args['env'] = Environment(args, idx=1, sA=0.0, sB=s)  # immediately swaps at gen 0
                args['aToB'] = aToB
            else:
                args['env'] = Environment(args, idx=0, sA=s, sB=0.0) # immediately swaps at gen 0
                #args['aToB'] = aToB
                args['aToB'] = 1.0-aToB # reverse initial freq in odd-numbered pops
            pop = Pop(args, idx)
            self.pops.append(pop)

    def __repr__(self):
        # for printing the grid info in pretty format
        # this includes printing each population
        ret = ""
        for idx, pop in enumerate(self.pops):
            ret = ret + f"{pop}"
            
        if idx < len(self.pops)-1:
            ret = ret + "\n"
        return ret

    def migrate(self):
        # apply migration between populations - in frequency space
        if self.numPops == 1: # there's just one population, skip migration
            return

        # For each pop, accumulate inbound and outbound counts of each type, in frequency space
        # There are numPops * numClades * numTypes distinct counts: For each pop, numClades clades. For each clade, numTypes types (A and B).
        numTypes = 2 # A and B
        current = np.zeros(shape = (self.numPops, self.numClades, numTypes)) # current residents
        mig = np.zeros(shape = (self.numPops, self.numClades, numTypes)) # migration rate
        for j, pop in enumerate(self.pops):
            for i, clade in enumerate(pop.clades):
                current[j][i][0] = clade.nA
                current[j][i][1] = clade.nB
                mig[j][i][0] = clade.mu
                mig[j][i][1] = clade.mu
        #print("current\n", current)
        #print("mig\n", mig)

        # go is the frequencies of how many leave each pop
        go = np.multiply(current, mig)
        #print("go\n", go)

        # final is current frequencies, minus those that go
        final = current - go

        # this fraction go to each other (i!=j) pop
        go = go/float(self.numPops-1.0)
        for i in range(self.numPops):
            for j in range(self.numPops):
                if i!=j:
                    final[i] = final[i] + go[j]

        #print("final\n", final)

        # copy final back into pop
        for j, pop in enumerate(self.pops):
            for i, clade in enumerate(pop.clades):
                clade.nA = final[j][i][0]
                clade.nB = final[j][i][1]

    def one_gen(self, gen, epochGen):
        # one generation of evolution for the whole grid
        # swap environment either:
        #   at gen==0
        #   or, if not args.stochasticEnv, then every epochGen generations *exactly*.
        #   or, if args.stochasticEnv, then every epochGen generations *in expectation*.
        if gen == 0 or \
            (not self.args.stochasticEnv and (gen % epochGen == 0)) or \
            (self.args.stochasticEnv and (random.random() < 1./epochGen)):

            if self.printStats:
                print(f"gen:{gen} epoch:{int(gen / epochGen)}")# envt:{envt}")
                #print(f"swap environment to {'A' if self.pops[0].env.idx == 1 else 'B'}")
            
            for pop in self.pops:
                # swap environment in each pop
                pop.env.swap(gen)

            if self.printStats:
                print(self) # start of every epoch

        # the mode is fixed for each run - either 'mutation' or 'migration'
        if self.mode == 'mutation':
            for pop in self.pops:
                pop.one_gen_mutation_mode()
        elif self.mode == 'migration':
            # migration mode is tricker because migration happens between populations
            for pop in self.pops:
                pop.update_gen_data() #log data
                pop.select() # converts from int to freq
            if self.args.includeMutation:
                self.migrate() # works on freq
            for pop in self.pops:
                pop.normalize() # converts from freq to int
        else:
            raise Exception(f"invalid mode: {self.mode}")


default_args = DotAccessibleDict({'numPops': 1, 'N': 1e5, 'minMu':2, 'numClades': 1, 'aToB': 0.5, 'T':1e3, 's':0.1, \
                'numEpochs':1, 'numRuns':1, 'includeMutation': True, 'stochasticEnv': False, \
                'plotLog': False, 'plotAB': True, 'plotKFit': False, 'hideB': False, 'printStats': False})

def round_to_int(x, name):
    rounded = round(x)
    if rounded != x:
        print(f"note: rounding {name} to {rounded}!")
    return rounded

def one_run(args):
    grid = Grid(args)

    for gen in range(int(args.numEpochs*args.T)):
        grid.one_gen(gen, args.T)
    
    run = DotAccessibleDict()
    run.args = args
    run.grid = grid
    
    return run

def check_args(default_args, override_args):
    #print("default_args:", default_args)
    #print("override_args:", override_args)

    unknown_keys = set(override_args.keys()) - set(default_args.keys())
    if unknown_keys:
        raise KeyError(f"Unknown argument keys: {sorted(unknown_keys)}")

    args = DotAccessibleDict(default_args.copy())
    args.update(override_args)
    #print("args:", args)

    # check constraints on arguments
    assertx(args.N <0 or args.N >=10, "N must be -1 (which means infinite pop) or >=10")
    assertx(args.N <= 1e9, "N must be <= 1e9")
    assertx(args.minMu >=2 and args.minMu <=8, "minMu must be between 2 and 8")
    assertx(args.numClades >=1 and args.minMu+args.numClades <=9, "numClades must be >= 1, and minMu+numClades must be <=9")
    assertx(args.numPops >= 1, "numPops must be 1 for mutation mode, or >=1 for migration mode")
    assertx((args.aToB >= 0 and args.aToB <= 1.0), "aToB must be between 0.0 and 1.0")

    args['mode'] = 'mutation' if args['numPops'] == 1 else 'migration' # might change this later
    args.N = round_to_int(args.N, "N")
    args.numPops = round_to_int(args.numPops, "numPops")
    args.minMu = round_to_int(args.minMu, "minMu")
    args.numClades = round_to_int(args.numClades, "numClades")
    args.s = args.s
    args.T = args.T = round_to_int(args.T, "T")
    args.numEpochs = round_to_int(args.numEpochs, "numEpochs")
    args.numRuns = round_to_int(args.numRuns, "numRuns")
    
    return args


#timeit
def evosim_run(override_args=DotAccessibleDict()):
    global default_args
    # numPops: number of populations
    # N: maximum population size (N <0 means infinite pop)
    # minMu: minimum mutator allele (2 means M2, 3 means M
    # numClades: number of clades (mutator alleles). Mutation rates are 10^-minMu, 10^-(minMu+1), ... 10^-(minMu+numClades-1)
    # aToB: initial fraction of A genotype
    # T: number of generations per epoch
    # numEpochs: number of epochs
    # s: selection coefficient
    # includeMutation: if True, perform mutation step
    # plotLog: if True, use log scale for plots
    # plotAB: if True, plot A and B genotypes separately; if False, plot total (A+B)
    # printStats: if True, print grid stats at start of each epoch
    # numRuns: number of independent simulation runs to perform

    args = check_args(default_args, override_args)
    N, numPops, minMu, numClades, s, T, numEpochs, aToB, numRuns = \
        args.N, args.numPops, args.minMu, args.numClades, args.s, args.T, args.numEpochs, args.aToB, args.numRuns
        
    if N>0:
        Ns = abs(N*s)
        print(f"Ns:{Ns} is low < 10.0 (selection is WEAK relative to drift)" if Ns < 10.0 else f"Ns:{Ns} is high >= 10.0 (selection is STRONG relative to drift)")
    else:
        print(f"Infinite population size N (there is no drift)")
    sT = abs(s*T)
    print(f"sT:{sT} is low < 20.0 (fixation WON'T usually happen in one epoch)" if sT < 20.0 else f"sT:{sT} is high >= 20.0 (fixation CAN happen in one epoch (ignoring mutation))")
    
    if args.mode == 'migration':
        print(f"Running in MIGRATION mode with {numPops} populations")
        
    # numRuns must be 1 with infinite pop, because runs are identical
    if N < 0:
        args.numRuns = 1
        print("Note: setting numRuns to 1 because N < 0 (infinite population)")
    
    ## run numRuns independent simulations
    runs = []
    if args.numRuns == 1:
        # fallback to sequential for a single run
        for r in range(args.numRuns):
            if args.printStats:
                print(f"=== Starting run {r+1} of {args.numRuns} ===")
            run = one_run(copy.deepcopy(args))
            runs.append(run)
    else:
        nproc = min(args.numRuns, os.cpu_count() or 1)
        if True: #args.printStats:
            print(f"Using {nproc} processes for {args.numRuns} runs.")
        # deepcopy args for each run so worker-local mutations don't interfere
        arg_list = [copy.deepcopy(args) for _ in range(args.numRuns)]
        ctx = multiprocessing.get_context("spawn")
        with ctx.Pool(processes=nproc) as pool:
            runs = pool.map(one_run, arg_list)

    return args, runs

def evosim_plot(args, runs):
    N, numPops, minMu, numClades, s, T, numEpochs, aToB, numRuns = \
    args.N, args.numPops, args.minMu, args.numClades, args.s, args.T, args.numEpochs, args.aToB, args.numRuns

    strN = "Inf" if N < 0 else f"{e_format(N)}"
    colors = [col.ColorConverter.to_rgb(x) for x in ["red", "green", "blue", "magenta", "cyan", "#ff9933", "black"]]
    
    xlog = args.plotLog
    ylog = args.plotLog
    
    lw = 1
    lsA = "-"
    lsB = "--"

    if True: # plot counts over time
        labels = []
        handles = []

        fig, axes = plt.subplots(numPops, layout='constrained', figsize=(4.0, numPops*4.0)) # 6.4x numPops*4.8.
        if numPops == 1: axes = [axes]
        fig.suptitle(f"N={strN}, s={e_format(s)}, T={e_format(T)}"+f" ({numRuns} runs)" if numRuns>1 else "")

        for r in range(numRuns): # each independent run
            grid = runs[r].grid
            args = runs[r].args

            for idx, pop in enumerate(grid.pops): # each deme in the grid
                ax = axes[idx]
                for idx2, clade in enumerate(pop.clades): # each clade in each pop
                    # Fixed colors: M2 is red, M3 is green, M4 is blue, M5 is magenta, M6 is cyan, M7 is orange, 
                    # darker shade for allele A, lighter shade for allele B
                    shades = [scale_lightness(colors[clade.m + minMu - 2], scale) for scale in [0.5, .75, 1., 1.25, 1.5]]
                    if idx==0: # with handles and labels
                        if args.plotAB: # break out A and B
                            # plot A only
                            handles.append(ax.plot(clade.countsA, color=shades[1], linestyle=lsA, linewidth=lw)[0])
                            if r == 0: # only add labels once
                                labels.append(r"$\mathit{M}_{{%d}},A$" % (clade.m + minMu))
                            if not args.hideB:
                                # plot B as well as A
                                handles.append(ax.plot(clade.countsB, color=shades[3], linestyle=lsB, linewidth=lw)[0])
                                if r == 0: # only add labels once
                                    labels.append(r"$\mathit{M}_{{%d}},B$" % (clade.m + minMu))
                        else:
                            handles.append(ax.plot(clade.counts, color=shades[2], linestyle=lsA, linewidth=lw)[0]) # note [0]
                            if r == 0: # only add labels once
                                labels.append(r"$\mathit{M}_{{%d}},(A+B)$" % (clade.m + minMu))
                    else: # just plot, no handles or labels
                        if args.plotAB:
                            ax.plot(clade.countsA, color=shades[1], linestyle=lsA, linewidth=lw)
                            if not args.hideB:
                                # plot B as well as A
                                ax.plot(clade.countsB, color=shades[3], linestyle=lsB, linewidth=lw)
                        else:
                            ax.plot(clade.counts, color=shades[2], linestyle=lsA, linewidth=lw)
                            
                # formats y axis labels nicely
                # thanks copilot!
                from matplotlib.ticker import FuncFormatter
                ax.get_yaxis().set_major_formatter(FuncFormatter(
                    lambda x, pos: "0" if x == 0 else (
                    "{:.0f}".format(x) if (abs(x) >= 1 and abs(x) < 10000) else
                    ("{:.2g}".format(x) if abs(x) < 1 and abs(x) < 10000 else
                    "{:.1e}".format(x).replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")))
                ))

                if idx != numPops-1: # all except bottom plot get no x axis labels
                    axes[idx].tick_params(labelbottom=False)
                    #axes[idx].locator_params(axis='x', nbins=5)  # just put 5 major tics

                if numPops > 1:
                    ax.set_ylabel(f"pop {idx}")
                else:
                    ax.set_ylabel("N" if N>0 else "frequency")
                
                if xlog:
                    ax.set_xscale("log", nonpositive='mask')

                if ylog:
                    ax.set_yscale("log", nonpositive='mask')
                    if N<0: # for infinite pop
                        ax.set_ylim(1e-8, 1.5)
                    else: 
                        ax.set_ylim(1, 1.5*N) #(N*1e-6, 1.5*N)
                else:
                    ax.set_ylim(-0.001*abs(N), 1.005*abs(N))
                    
                # print vertical lines and labels for environment swaps
                # (ylim has to be set before this)
                if r == 0:
                    for swap_gen, env_idx in pop.env.active_env:
                        ax.axvline(x=swap_gen+(1 if xlog else 0), color='gray', linestyle='--', linewidth=1)
                        label = 'A' if env_idx == 'A' else 'B'
                        y_pos = ax.get_ylim()[1] # *0.96
                        ax.text(swap_gen+(1 if xlog else 0), y_pos, label, color='gray', ha='center', va='bottom', fontsize=10) # -ax.get_xlim()[1]*0.015 + 
                
            plt.xlabel("generations")

        #plt.figlegend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        #print(f"pop {idx} env: {pop.env.active_env}")

        #ax = axes[numPops]
        #ax.plot(grid.pops[0].envs)

        #plt.figlegend(handles=handles, labels=labels, loc='outside right center')
        plt.figlegend(handles=handles, labels=labels,
            loc='upper center', bbox_to_anchor=(0.5, -0.05),
            fancybox=True, shadow=True, ncol=3)

        #plt.tight_layout()
        #plt.savefig(f"output/plotx.png", dpi=300)
        plt.show()
    
    ###
    
    # compute and plot fit, surv, surprisal for each clade, at each time, averaged over runs
    # fitness stats are for entire clades, not types. Because we might start with 0 of any type
    if args.plotKFit:
        labels = []
        handles = []

        allp = np.zeros((numPops, numClades, numRuns, T*numEpochs))
        for r in range(numRuns):
            grid = runs[r].grid
            args = runs[r].args
            for idx, pop in enumerate(grid.pops):
                for idx2, clade in enumerate(pop.clades):
                    if N < 0:
                        clade.fit = np.array(clade.counts) / clade.counts[0]  # fitness (i.e., relative to initial count) for clade idx2 in pop idx over t
                        clade.surv = np.array([1.0 if x>0 else 0.0 for x in clade.counts]) # survival for clade idx2 in pop idx over t
                    else:
                        clade.fit = np.array(clade.counts) / clade.counts[0]  # fitness (i.e., relative to initial count) for clade idx2 in pop idx over t
                        clade.surv = np.array([1.0 if x>0 else 0.0 for x in clade.counts]) # survival for clade idx2 in pop idx over t

        # now compute across-run averages
        for idx, pop in enumerate(grid.pops):
            for idx2, clade in enumerate(pop.clades):
                clade.avgFit = np.zeros_like(clade.fit)
                clade.avgSurv = np.zeros_like(clade.surv)
                for r in range(numRuns):
                    grid = runs[r].grid
                    pop = grid.pops[idx]
                    clade_r = pop.clades[idx2]
                    clade.avgFit += clade_r.fit
                    clade.avgSurv += clade_r.surv
                    allp[idx][idx2][r] = clade_r.surv # the prob that clade idx2 survives at time t in run r
                clade.avgFit /= float(numRuns)
                clade.avgSurv /= float(numRuns)
                clade.surprisal = -clade.avgSurv*np.log2(clade.avgSurv + 1e-12) - (1.0 - clade.avgSurv)*np.log2(1.0 - clade.avgSurv + 1e-12)
            
        # compute joint surprisal over all clades for all run
        joint_surprisal = np.zeros((numPops, T*numEpochs))
        numStates = 2**numClades
        for t in range(T*numEpochs):
            for idx, pop in enumerate(grid.pops):
                joint_surv = np.zeros(numStates) # for this whole pop at time t
                for state in range(numStates):
                    prob = 1.0 # will hold the probability of this joint state at time t
                    for idx2, clade in enumerate(pop.clades): # each clade in each pop
                        if (state >> idx2) & 1: # clade idx2 survives in this state
                            prob = prob * allp[idx][idx2][:, t] # the prob that clade idx2 survives at time t in all runs
                        else: # clade idx2 does NOT survive in this state
                            prob = prob * (1.0 - allp[idx][idx2][:, t]) # the prob that clade idx2 does NOT survive at time t in all runs
                    joint_surv[state] = np.mean(prob) # average over runs
                joint_surprisal[idx][t] = -np.sum(joint_surv * np.log2(joint_surv + 1e-12)) 


        # now plot the avg fit (left axis) and surv (right axis) data, with one row for each pop
        fig, axes = plt.subplots(numPops, 3, layout='constrained', figsize=(3*4.0, numPops*4.0)) # Adjust figsize as needed
        fig.suptitle(f"N={strN}, s={e_format(s)}, T={e_format(T)}{f" ({numRuns} runs)" if numRuns>1 else ""}")
        for idx, pop in enumerate(grid.pops): # each deme in the grid
            if numPops == 1:
                ax1, ax2, ax3 = axes[0], axes[1], axes[2]
            else:
                ax1, ax2, ax3 = axes[idx][0], axes[idx][1], axes[idx][2]

            lw = 1
            lsA = "-"
            lsB = "--"
            for idx2, clade in enumerate(pop.clades): # each clade in each pop
                # Fixed colors: M2 is red, M3 is green, M4 is blue, M5 is magenta, M6 is cyan, M7 is orange, 
                # darker shade for allele A, lighter shade for allele B
                shades = [scale_lightness(colors[clade.m + minMu - 2], scale) for scale in [0.5, .75, 1., 1.25, 1.5]]
                handles.append(ax1.plot(clade.avgFit, color=shades[2], linestyle=lsA, linewidth=lw)[0]) # note [0]
                labels.append(r"$\mathit{M}_{{%d}},(A+B)$" % (clade.m + minMu))
                ax2.plot(clade.avgSurv, color=shades[2], linestyle=lsA, linewidth=lw)
                ax3.plot(clade.surprisal, color=shades[2], linestyle=lsA, linewidth=lw)
            # plot joint surprisal
            ax3.plot(joint_surprisal[idx], color='#ff4d00', linestyle=':', linewidth=lw, label='joint')

            ax1.set_ylabel("W_k(0)")
            ax2.set_ylabel("S_k(0)")
            ax3.set_ylabel("H(k)")
            if idx != numPops-1: # all except bottom plot get no x axis labels
                ax1.tick_params(labelbottom=False)
                ax2.tick_params(labelbottom=False)
                ax3.tick_params(labelbottom=False)
            ax1.locator_params(axis='x', nbins=5)  # just put 5 major tics
            ax2.locator_params(axis='x', nbins=5)  # just put 5 major tics
            ax3.locator_params(axis='x', nbins=5)  # just put 5 major tics
            ax1.set_xlabel("generation (k)")
            ax2.set_xlabel("generation (k)")
            ax3.set_xlabel("generation (k)")
        
            if xlog:
                ax1.set_xscale("log", nonpositive='mask')
                ax2.set_xscale("log", nonpositive='mask')
                ax3.set_xscale("log", nonpositive='mask')
            if ylog:
                ax1.set_yscale("log", nonpositive='mask')
                
            ax2.set_ylim(-0.05, 1.05)
            ax3.set_ylim(-0.05, math.log2(numClades + 1))

        plt.figlegend(handles=handles, labels=labels,
                    loc='upper center', bbox_to_anchor=(0.5, -0.05),
                    fancybox=True, shadow=True, ncol=3)

    
def evosim(override_args=DotAccessibleDict()):
    args, runs = evosim_run(override_args)
    evosim_plot(args, runs)

if __name__ == '__main__':
    evosim()
