import matplotlib.colors as col
import colorsys
def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

import string
class MyFormatter(string.Formatter):
    def format_field(self, value, format_spec):
        ss = string.Formatter.format_field(self,value,format_spec)
        if format_spec.endswith('e'):
            if ( 'e' in ss):
                mantissa, exp = ss.split('e')
                return mantissa + 'e'+ exp[0] + str(int(exp[1:]))
        return ss


# strs = MyFormatter().format('{0:.1e}', s)
# strsT = MyFormatter().format('{0:.1e}',sT)
# strT = MyFormatter().format('{0:.1e}',epochGen)
