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

from collections import UserDict
class DotAccessibleDict(UserDict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{key}'")
