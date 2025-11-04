from scipy.stats import binom
import matplotlib.pyplot as plt
import numpy as np
from helpers import *

class Environment():
    def __init__(self, sA, sB, idx):
        self.idx = idx
        self.sA = sA
        self.sB = sB
        self.active_env = [] # pairs of (gen, env)

    def __repr__(self):
        return f"idx:{self.idx} sA:{self.sA} sB:{self.sB}"

    def swap(self, gen): # swap which allele is favored
        self.idx = (self.idx+1)%2
        tmp = self.sA
        self.sA = self.sB
        self.sB = tmp

        # log data for plotting
        self.active_env.append((gen, "A" if self.idx == 0 else "B"))

class Clade:
    # a clade contains two genotypes: A and B, with the same fixed mutation rate
    def __init__(self, m, args):
        # initialize one clade
        popSize = args['popSize']
        aToB = args['aToB']
        numClades = args['numClades']
        minMu = args['minMu']

        self.m = m # m is 0, 1, 2, ... for clades labelled M2, M3, M4, ...
        self.mu = 10.0**(-1.0 * (m+minMu)) # mutation rate
        self.label = f"m{m+minMu}" # clade labels M2, M3, M4, ... for printing

        fracA = 1.0 / numClades * aToB  # initial fraction of A genotype
        fracB = 1.0 / numClades * (1.0 - aToB) # initial fraction of B genotype
        if popSize < 0:  # infinite pop - floating point "frequency" of each genotype
            self.nA = fracA
            self.nB = fracB
        else: # finite pop - integer count of each genotype
            self.nA, self.nB = binom.rvs(n=round(popSize), p=[fracA, fracB])

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
        ret = ""
        ret = ret + f"mu:{self.mu:.2e} nA: {self.nA:10.2f} nB:{self.nB:10.2f}, (n:{self.nA+self.nB:10.2f})"
        return ret

    def select(self, env):
        # apply selection, in frequency space
        # nA and nB in each clade are non-integer after this, and sumN != popSize
        self.nA = self.nA*(1.0+env.sA)
        self.nB = self.nB*(1.0+env.sB)

    def mutate(self):
        # apply mutation, in frequency space
        # nA and nB in each clade are non-integer after this, and sumN != popSize
        aToB = self.nA*self.mu
        AtoA = self.nA*(1.0-self.mu)
        BtoA = self.nB*self.mu
        BtoB = self.nB*(1.0-self.mu)
        #print(f"AtoA:{AtoA} aToB:{aToB} BtoB:{BtoB} BtoA:{BtoA}")
        self.nA = AtoA + BtoA
        self.nB = aToB + BtoB

    def normalize(self, sumN, popSize):
        # convert from frequency space back to integer counts
        fracA = float(self.nA) / sumN # fraction of A
        fracB = float(self.nB) / sumN # fraction of B
        if popSize < 0: # this means it's an infinite pop
            # infinite pop - don't convert to integer
            self.nA = fracA # total N of this type
            self.nB = fracB # total N of this type
        else:
            # finite pop - convert to integer
            # save computing binom if frac is 0.0 anyway, just set count to 0
            self.nA = 0 if fracA == 0.0 else binom.rvs(n=round(popSize), p=fracA) # total N of this type
            self.nB = 0 if fracB == 0.0 else binom.rvs(n=round(popSize), p=fracB) # total N of this type


class Pop:
    # a population contains multiple clades
    def __init__(self, args):
        # initialize one population
        self.popSize = args['popSize']
        self.env = args['env']
        numClades = args['numClades']

        # the list of clades in this population
        self.clades = []
        for m in range(0, numClades):
            clade = Clade(m, args)
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
        ret = ""
        # this includes printing each clade
        for idx, clade in enumerate(self.clades):
            ret = ret + f"  clade {clade.label}: {clade}\n"
        ret = ret + f"  sumN: {self.sumN():10.2f}\n"
        ret = ret + f"  env: {self.env}\n"
        return ret

    def sumN(self): # counts individuals in this population
        sumN = 0.0
        for clade in self.clades:
            sumN = sumN + (clade.nA + clade.nB)
        assert sumN > 0.0, f"uh oh, sumN for a population should be > 0. Probably popSize {self.popSize} is too small."
        return sumN

    def select(self):
        # nA and nB in each clade are frequencies after this, and sum != popSize
        for clade in self.clades:
            clade.select(self.env)

    def mutate(self):
        # nA and nB in each clade are frequencies after this, and sum != popSize
        for clade in self.clades:
            clade.mutate()

    def normalize(self):
        # convert from frequency space back to integer counts
        sumN = self.sumN()
        for clade in self.clades:
            clade.normalize(sumN, self.popSize)

    # def one_gen_with_mutation(self, gen):
    #     self.update_gen_data()
    #     self.select() # converts from int to freq
    #     self.mutate() # works on freq
    #     self.normalize() # converts from freq to int


class Grid():
    # A grid contains multiple populations
    # The point is that you can have migration between populations
    def __init__(self, args):
        # initialize the grid of populations (there is only one grid per run)
        self.numPops = args['numPops']
        self.numClades = args['numClades']
        self.mode = args['mode'] # 'mutation' or 'migration
        s = args['s']
        aToB = args['aToB']

        self.pops = []
        for idx in range(self.numPops):
            if idx%2 == 0:
                args['env'] = Environment(idx=1, sA=0.0, sB=s)  # immediately swaps at gen 0
                args['aToB'] = aToB
            else:
                args['env'] = Environment(idx=0, sA=s, sB=0.0) # immediately swaps at gen 0
                args['aToB'] = 1.0-aToB
            pop = Pop(args)
            self.pops.append(pop)

    def __repr__(self):
        # for printing the grid info in pretty format
        ret = ""
        # this includes printing each population
        for idx, pop in enumerate(self.pops):
            ret = ret + f" pop {idx}:\n{pop}"
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
        # swap environments every epochGen generations
        if gen % epochGen == 0:
            print(f"gen:{gen} epoch:{int(gen / epochGen)}")# envt:{envt}")

            print("swap environments")
            for pop in self.pops:
                pop.env.swap(gen)

            print(self) # start of every epoch

        # the mode is fixed for each run - either 'mutation' or 'migration'
        if self.mode == 'mutation':
            for pop in self.pops:
                pop.update_gen_data() # log data
                pop.select()  # converts from int to freq
                pop.mutate()  # works on freq
                pop.normalize()  # converts from freq to int

        elif self.mode == 'migration':
            for pop in self.pops:
                pop.update_gen_data() #log data
                pop.select() # converts from int to freq

            self.migrate() # works on freq

            for pop in self.pops:
                pop.normalize() # converts from freq to int
        else:
            raise Exception(f"invalid mode: {self.mode}")

def evosim(user_args={}):
    # numPops: number of populations
    # popSize: maximum population size (popSize <0 means infinite pop)
    # s: selection coefficient
    # minMu: minimum mutator allele (2 means M2, 3 means M
    # numClades: number of clades (mutator alleles). Mutation rates are 10^-minMu, 10^-(minMu+1), ... 10^-(minMu+numClades-1)
    # aToB: initial fraction of A genotype
    # numEpochs: number of epochs
    # T: number of generations per epoch
    # mode: 'mutation' or 'migration' - whether to simulate mutation within one pop, or migration between multiple pops

    default_args = {'numPops': 1, 'popSize': 1e9, 's':0.1, 'minMu':2, 'numClades': 3, 'aToB': 0.0, 'numEpochs':10, 'T':1e5}
    default_args['mode'] = 'mutation' if default_args['numPops'] == 1 else 'migration'

    print("user_args:", user_args)
    print("default_args:", default_args)
    # merge args with default_args, overriding default_args
    unknown_keys = set(user_args.keys()) - set(default_args.keys())
    if unknown_keys:
        raise KeyError(f"Unknown argument keys: {sorted(unknown_keys)}")
    args = default_args.copy()
    args.update(user_args)
    print("merged_args:", args)

    popSize = args['popSize']
    numPops = args['numPops']
    minMu = args['minMu']
    mode = args['mode']
    numClades = args['numClades']
    s = args['s']
    T = args['T']
    numEpochs = args['numEpochs']

    xlog = False
    ylog = False

    # check constraints on arguments
    assert popSize <0 or popSize >=10, "popSize must be <0 (which means infinite pop) or >=10"
    assert minMu >=2 and minMu <=8, "minMu must be between 2 and 8"
    assert numClades >=2 and numClades <=9, "numClades must be between 2 and 9"
    assert mode in ['mutation', 'migration'], "mode must be 'mutation' or 'migration'"
    assert ((numPops == 1 and mode=='mutation') \
        or (numPops >=1)), "numPops must be 1 for mutation mode, or >=1 for migration mode"
    if s*T<5.0:
        print(f"s {s} is pretty small for given T {T} (sT should be >=5.0)")
    print(f"sT:{s*T:.1f}")
    if popSize>0 and s>0: print(f"Ns:{e_format(popSize*s)}")  

    ## run simulation
    
    grid = Grid(args)

    for gen in range(int(numEpochs*T)):
        grid.one_gen(gen, T)

    ## plot simulation

    strN = "Inf" if popSize < 0 else f"{e_format(popSize)}"
    colors = [col.ColorConverter.to_rgb(x) for x in ["red", "green", "blue", "magenta", "cyan", "yellow", "black"]]

    fig, axes = plt.subplots(numPops, layout='constrained', figsize=(6.4, numPops*4.8)) # 6.4x4.8.
    if numPops == 1: axes = [axes]
    if 1:
        fig.suptitle(f"N={strN}, T={e_format(T)}, s={e_format(s)}")

    labels = []
    handles = []
    plot_total = True
    
    for idx, pop in enumerate(grid.pops):
        ax = axes[idx]
        for idx2, clade in enumerate(pop.clades):
            # Fixed colors: M2 is red, M3 is green, M4 is blue, M5 is magenta, M6 is cyan, M7 is yellow, 
            # darker shade for allele A, lighter shade for allele a
            shades = [scale_lightness(colors[clade.m + minMu - 2], scale) for scale in [0.5, .75, 1., 1.25, 1.5]]
            if idx==0: # with handles and labels
                if plot_total:
                    handles.append(ax.plot(clade.counts, color=shades[3], linestyle="-")[0]) # note [0]
                    labels.append(r"$\mathit{M}_{{%d}},(A+a)$" % (clade.m + minMu))
                else:
                    handles.append(ax.plot(clade.countsA, color=shades[2], linestyle="-")[0])
                    labels.append(r"$\mathit{M}_{{%d}},A$" % (clade.m + minMu))
                    handles.append(ax.plot(clade.countsB, color=shades[4], linestyle=":")[0])
                    labels.append(r"$\mathit{M}_{{%d}},a$" % (clade.m + minMu))
            else:
                if plot_total:
                    ax.plot(clade.counts, color=shades[3], linestyle="-")
                else:
                    ax.plot(clade.countsA, color=shades[2], linestyle="-")
                    ax.plot(clade.countsB, color=shades[4], linestyle=":")
        ax.locator_params(axis='x', nbins=5)  # just put 5 major tics
        
        if xlog:
            ax.set_xscale("log", nonpositive='mask')

        if ylog:
            ax.set_yscale("log", nonpositive='mask')
            if popSize<0: # for infinite pop
                ax.set_ylim(1e-8, 1.0)
            else: 
                ax.set_ylim(popSize*1e-6, 1.5*popSize)
        for swap_gen, env_idx in pop.env.active_env:
            ax.axvline(x=swap_gen+(1 if xlog else 0), color='gray', linestyle='--', linewidth=1)
            label = 'A' if env_idx == 'A' else 'B'
            y_pos = ax.get_ylim()[1]
            ax.text(swap_gen+(1 if xlog else 0), y_pos, label, color='gray', ha='center', va='bottom', fontsize=10)

        plt.ylabel("N" if popSize>0 else "frequency")
        plt.xlabel("generations")
        [axes[idx].tick_params(labelbottom=False) for idx in range(numPops-1)]

        if mode=='migration':
            ax.set_ylabel(f"pop {idx}")

        #plt.figlegend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        
        print(f"pop {idx} env: {pop.env.active_env}")

    #ax = axes[numPops]
    #ax.plot(grid.pops[0].envs)

    plt.figlegend(handles=handles, labels=labels, loc='outside right center')
    #plt.tight_layout()
    plt.savefig(f"output/plotx.png", dpi=300)
    plt.show()


if __name__ == '__main__':
    evosim()
