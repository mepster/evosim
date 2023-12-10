from scipy.stats import binom
import matplotlib.pyplot as plt
import numpy as np
from helpers import *

class Environment():
    def __init__(self, sA, sB):
        self.sA = sA
        self.sB = sB

    def __repr__(self):
        return f"sA:{self.sA} sB:{self.sB}"

    def swap(self):
        tmp = self.sA
        self.sA = self.sB
        self.sB = tmp


class Clade:
    def __init__(self, nA, nB, mu, label):
        self.nA = nA
        self.nB = nB
        self.mu = mu
        self.label = label
        self.init_data()

    def __repr__(self):
        ret = ""
        ret = ret + f"nA: {self.nA:10.2f} nB:{self.nB:10.2f}, (n:{self.nA+self.nB:10.2f}) mu:{self.mu}"
        return ret

    def select(self, env):
        # nA and nB in each clade are non-integer after this, and sumN != maxN
        self.nA = self.nA*(1.0+env.sA)
        self.nB = self.nB*(1.0+env.sB)

    def mutate(self):
        # nA and nB in each clade are non-integer after this, and sumN != maxN
        AtoB = self.nA*self.mu
        AtoA = self.nA*(1.0-self.mu)
        BtoA = self.nB*self.mu
        BtoB = self.nB*(1.0-self.mu)
        #print(f"AtoA:{AtoA} AtoB:{AtoB} BtoB:{BtoB} BtoA:{BtoA}")
        self.nA = AtoA + BtoA
        self.nB = AtoB + BtoB

    def normalize(self, sumN, maxN):
        fracA = float(self.nA)/sumN # fraction of this clade's A out of entire pop
        fracB = float(self.nB)/sumN # fraction of this clade's B out of entire pop
        if maxN < 0: # infinite pop
            self.nA = fracA
            self.nB = fracB
        else:
            # save computing binom if frac is 0.0 anyway
            self.nA = 0 if fracA == 0.0 else int(binom.rvs(n=maxN, p=fracA))
            self.nB = 0 if fracB == 0.0 else int(binom.rvs(n=maxN, p=fracB))

    def init_data(self):
        self.countsA = []
        self.countsB = []
        self.counts = []

    def update_data(self):
        self.countsA.append(self.nA)
        self.countsB.append(self.nB)
        self.counts.append(self.nA+self.nB)


class Pop:
    def __init__(self, args):
        # fracA is initial fraction of A in this clade
        minMu = args['minMu']
        numClades = args['numClades']
        fracA = args['fracA']
        self.maxN = args['maxN']
        self.env = args['env']

        self.clades = []
        for idx, m in enumerate(range(minMu, minMu+numClades)):
            if self.maxN < 0: # infinite pop
                nA = 1.0 / numClades * fracA
                nB = 1.0 / numClades * (1.0 - fracA)
            else:
                nA = binom.rvs(n=self.maxN, p=1.0 / numClades * fracA)
                nB = binom.rvs(n=self.maxN, p=1.0 / numClades * (1.0 - fracA))

            clade = Clade(nA=nA, nB=nB, mu=10**(-1.0 * m), label=str(idx+minMu))
            self.clades.append(clade)

    def __repr__(self):
        sumN = 0
        ret = ""
        for idx, clade in enumerate(self.clades):
            ret = ret + f"  clade {idx}: {clade}\n"
            sumN = sumN + clade.nA + clade.nB
        ret = ret + f"  sumN: {sumN:.2f}\n"
        ret = ret + f"  env: {self.env}\n"
        return ret

    def select(self):
        # nA and nB in each clade are non-integer after this, and sum != maxN
        for clade in self.clades:
            clade.select(self.env)

    def mutate(self):
        # nA and nB in each clade are non-integer after this, and sum != maxN
        for clade in self.clades:
            clade.mutate()

    def normalize(self):
        sumN = 0
        for clade in self.clades:
            sumN = sumN + clade.nA + clade.nB
        for clade in self.clades:
            clade.normalize(sumN, self.maxN)

    def update_data(self):
        for clade in self.clades:
            clade.update_data()

    # def one_gen_with_mutation(self, gen):
    #     self.update_data()
    #     self.select() # converts from int to freq
    #     self.mutate() # works on freq
    #     self.normalize() # converts from freq to int


class Grid():
    def __init__(self, args):
        self.numPops = args['numPops']
        self.numClades = args['numClades']
        self.mode = args['mode'] # 'mutation' or 'migration
        s = args['s']
        fracA = args['fracA']

        self.pops = []
        for idx in range(self.numPops):
            if idx%2 == 0:
                args['env'] = Environment(sA=s, sB=0.0)  # immediately swaps at gen 0
                args['fracA'] = fracA
            else:
                args['env'] = Environment(sA=0.0, sB=s) # immediately swaps at gen 0
                args['fracA'] = 1.0-fracA
            pop = Pop(args)
            self.pops.append(pop)

    def __repr__(self):
        ret = ""
        for idx, pop in enumerate(self.pops):
            ret = ret + f" pop {idx}:\n{pop}"
        return ret

    def migrate(self):
        if self.numPops == 1: # no migration
            return

        # for each pop, accumulate inbound and outbound counts of each type
        # there are numPops * numClades * numTypes total counts
        numTypes = 2
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
        go = np.multiply(current, mig)
        #print("go\n", go)

        # this many leave each pop
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
        if gen % epochGen == 0:
            print(f"gen:{gen} epoch:{int(gen / epochGen)}")# envt:{envt}")

            print("swap environments")
            for pop in self.pops:
                pop.env.swap()

            print(self) # start of every epoch

        if self.mode == 'mutation':
            for pop in self.pops:
                pop.update_data()
                pop.select()  # converts from int to freq
                pop.mutate()  # works on freq
                pop.normalize()  # converts from freq to int
        elif self.mode == 'migration':
            for pop in self.pops:
                pop.update_data()
                pop.select() # converts from int to freq

            self.migrate() # works on freq

            for pop in self.pops:
                pop.normalize() # converts from freq to int
        else:
            raise Exception(f"invalid mode: {self.mode}")

def main():
    maxN = -1#1000000 # maxN = -1 means infinite pop
    numClades = 5

    if 1:
        mode = 'mutation'
        numPops = 1
    else:
        mode = 'migration'
        numPops = 2

    fracA = 0.5
    s = 0.001

    numEpochs = 8
    epochGen = 10000

    args = {'numPops': numPops, 'numClades': numClades, 'maxN': maxN, 's':s, 'fracA': fracA, 'minMu':2, 'mode':mode}
    grid = Grid(args)

    for gen in range(numEpochs*epochGen):
        grid.one_gen(gen, epochGen)
    print()

    # check sT and N/s
    sT = round(s*epochGen, 2)
    if sT<5.0:
        print(f"*** LOW sT:{sT:.1f}")
    else:
        print(f"sT:{sT:.1f}")
    if maxN>1: print(f"N/s:{e_format(maxN/s)}")

    print("mode:", args['mode'])

    ## plots

    strN = "Inf" if maxN==-1 else f"{e_format(maxN)}"
    colors = [col.ColorConverter.to_rgb(x) for x in ["darkred", "orange", "darkgreen", "navy", "darkviolet" ]]

    fig, axes = plt.subplots(numPops, layout='constrained', figsize=(6.4, numPops*4.8)) # 6.4x4.8
    if numPops == 1: axes = [axes]
    fig.suptitle(f"N={strN}, T={e_format(epochGen)}, s={e_format(s)}, {mode} mode")

    labels = []
    handles = []
    for idx, pop in enumerate(grid.pops):
        ax = axes[idx]
        for idx2, clade in enumerate(pop.clades):
            shades = [scale_lightness(colors[idx2], scale) for scale in [0.5, .75, 1., 1.25, 1.5]]
            if idx==0: # with handles and labels
                #handles.append(ax.plot(clade.counts, color=shades[3], linestyle="-")[0]) # note [0]
                #labels.append(f"m{clade.label} all")
                handles.append(ax.plot(clade.countsA, color=shades[2], linestyle="-")[0])
                labels.append(f"m{clade.label} A")
                handles.append(ax.plot(clade.countsB, color=shades[4], linestyle=":")[0])
                labels.append(f"m{clade.label} B")
            else:
                #ax.plot(clade.counts, color=shades[3], linestyle="-")
                ax.plot(clade.countsA, color=shades[2], linestyle="-")
                ax.plot(clade.countsB, color=shades[4], linestyle=":")
        # ax.xscale("log")
        ax.locator_params(axis='x', nbins=5)  # just put 5 major tics
        ax.set_yscale("log", nonpositive='mask')

        plt.xlabel("generations")
        [axes[idx].tick_params(labelbottom=False) for idx in range(numPops-1)]

        if mode=='migration':
            ax.set_ylabel(f"pop {idx}")

        #plt.figlegend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.figlegend(handles=handles, labels=labels, loc='outside right center')
    #plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
