import matplotlib.colors as col
import colorsys
def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

def e_format(num, m=1, e=1, signFlag=False):
    mantissa = float(f"{num:.{m}e}".split('e')[0])
    exponent = int(f"{num:{e}e}".split('e')[1])
    sign = ('+' if signFlag else '') if exponent>=0 else '-'
    return f"{mantissa}e{sign}{abs(exponent):0>{e}}"

# nums = [1.2345678*10.0**x for x in range(-3,3)]
# for num in nums:
#     print(e_format(num, 1, 1, False))

# def main_one_pop():
#     gen = 0
#
#     numClades = 5
#     maxN = -1#100000 # maxN = -1 means infinite pop
#     fracA = 0.5
#     minMu = 2
#     pop = Pop(numClades=numClades, maxN=maxN, fracA=fracA, minMu=minMu)
#
#     maxEpochs = 4
#     epochGen = 10000
#     s = 0.002
#
#     envt = 1 # immediately switches to 0
#     for gen in range(maxEpochs*epochGen):
#         if gen%epochGen == 0:
#             envt = (envt+1)%2
#             print(f"gen:{gen} epoch:{int(gen/epochGen)} envt:{envt}")
#
#         if envt == 0:
#             pop.one_gen(gen, sA=s, sB=0.0)
#         else:
#             pop.one_gen(gen, sA=0.0, sB=s)
#     print()
#
#     # check sT
#     sT = round(s*epochGen, 2)
#     if sT<5.0:
#         print(f"*** LOW sT:{sT}")
#     else:
#         print(f"sT:{sT}")
#
#     print(f"N/s:{maxN/s:.1e}")
#
#     strN = "Inf" if maxN==-1 else f"{maxN:.1e}"
#     plt.title(f"N={strN} sT={sT:.1e}") #s={s:.1e} T={epochGen:.1e}")#
#     if 1:
#         colors = [col.ColorConverter.to_rgb(x) for x in ["darkred", "darkorange", "darkgreen", "navy", "darkviolet" ]]
#         for idx, clade in enumerate(pop.clades):
#             shades = [scale_lightness(colors[idx], scale) for scale in [0, .5, 1, 1.5, 2]]
#             plt.plot(clade.counts, color=shades[3], linestyle="-", label=f"m{minMu+idx} all")
#             plt.plot(clade.countsA, color=shades[3], linestyle=":", label=f"m{minMu+idx} A")
#             plt.plot(clade.countsB, color=shades[3], linestyle="--", label=f"m{minMu+idx} B")
#     else:
#         list_of_As = [x.countsA for x in pop.clades]
#         allA = [sum(x) for x in zip(*list_of_As)]
#         plt.plot(allA, label="all A")
#         list_of_Bs = [x.countsB for x in pop.clades]
#         allB = [sum(x) for x in zip(*list_of_Bs)]
#         plt.plot(allB, label="all B")
#
#     #plt.xscale("log")
#     plt.locator_params(axis='x', nbins=5)
#     plt.ylim(bottom=0)
#     plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
#     plt.tight_layout()
#     plt.show()
