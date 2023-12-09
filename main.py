from scipy.stats import binom
import matplotlib.pyplot as plt

import matplotlib.colors as col
import colorsys
def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

class Clade:
    def __init__(self, nA, nB, mu):
        self.nA = nA
        self.nB = nB
        self.mu = mu
        self.init_data()

    def __repr__(self):
        ret = ""
        ret = ret + f"nA: {self.nA:.2f} nB:{self.nB:.2f}, (n:{self.nA+self.nB:.2f}) mu:{self.mu}"
        return ret

    def select(self, sA, sB):
        # nA and nB in each clade are non-integer after this, and sumN != maxN
        self.nA = self.nA*(1.0+sA)
        self.nB = self.nB*(1.0+sB)

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
        if maxN == -1: # infinite pop
            self.nA = fracA
            self.nB = fracB
        else:
            self.nA = int(binom.rvs(n=maxN, p=fracA))
            self.nB = int(binom.rvs(n=maxN, p=fracB))

    def init_data(self):
        self.countsA = []
        self.countsB = []
        self.counts = []

    def update_data(self):
        self.countsA.append(self.nA)
        self.countsB.append(self.nB)
        self.counts.append(self.nA+self.nB)

class Pop:
    def __init__(self, maxN, fracA=0.5, minMu=1, numClades=5):
        # fracA is initial fraction of A in this clade
        self.maxN = maxN

        self.clades = []
        for idx, m in enumerate(range(minMu, minMu+numClades)):
            if self.maxN == -1: # infinite pop
                nA = 1.0 / numClades * fracA
                nB = 1.0 / numClades * (1 - fracA)
            else:
                nA = binom.rvs(n=self.maxN, p=1.0 / numClades * fracA)
                nB = binom.rvs(n=self.maxN, p=1.0 / numClades * (1 - fracA))

            clade = Clade(nA=nA, nB=nB, mu=10**(-1.0 * m))
            self.clades.append(clade)

    def __repr__(self):
        sumN = 0
        ret = ""
        for idx, clade in enumerate(self.clades):
            ret = ret + f"clade {idx}: {clade}\n"
            sumN = sumN + clade.nA + clade.nB
        ret = ret + f"sumN: {sumN}\n"
        return ret

    def select(self, sA, sB):
        # nA and nB in each clade are non-integer after this, and sum != maxN
        for clade in self.clades:
            clade.select(sA, sB)

    def mutate(self):
        # nA and nB in each clade are non-integer after this, and sum != maxN
        for clade in self.clades:
            clade.mutate()

    def normalize(self, maxN):
        sumN = 0
        for clade in self.clades:
            sumN = sumN + clade.nA + clade.nB
        for clade in self.clades:
            clade.normalize(sumN, maxN)

    def update_data(self):
        for clade in self.clades:
            clade.update_data()

    def one_gen(self, gen, sA, sB):
        #print("gen:", gen)
        #print(self)
        self.update_data()
        self.select(sA, sB)
        self.mutate()
        self.normalize(self.maxN)

def main():
    gen = 0

    # maxN = -1 means infinite pop
    maxN = -1
    fracA = 0.1
    minMu = 2
    numClades = 5
    pop = Pop(maxN=maxN, fracA=fracA, minMu=minMu, numClades=numClades)

    maxEpochs = 2
    epochGen = 500
    s = 0.01

    for epoch in range(maxEpochs):
        print(f"epoch: {epoch}")
        for _ in range(epochGen):
            pop.one_gen(gen, sA=s, sB=0.0)
            gen=gen+1
        for _ in range(epochGen):
            pop.one_gen(gen, sA=0.0, sB=s)
            gen=gen+1

    sT = s*epochGen
    if sT<=5.0:
        print(f"LOW sT:{sT}")
    else:
        print(f"sT:{sT}")

    colors = [col.ColorConverter.to_rgb(x) for x in ["darkred", "darkorange", "darkgreen", "navy", "darkviolet" ]]
    for idx, clade in enumerate(pop.clades):
        shades = [scale_lightness(colors[idx], scale) for scale in [0, .5, 1, 1.5, 2]]
        plt.plot(clade.counts, color=shades[3], linestyle="-", label=f"m{minMu+idx} all")
        plt.plot(clade.countsA, color=shades[3], linestyle=":", label=f"m{minMu+idx} A")
        plt.plot(clade.countsB, color=shades[3], linestyle="--", label=f"m{minMu+idx} B")
    #plt.xscale("log")
    plt.ylim(bottom=0)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
