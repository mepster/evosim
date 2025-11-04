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

class DotAccessibleDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__()
        if args:
            if len(args) > 1:
                raise TypeError("Expected at most one positional argument")
            self.update(args[0])
        if kwargs:
            self.update(kwargs)

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{key}'")

    def __setattr__(self, key, value):
        if key.startswith('_'):
            super().__setattr__(key, value)
        else:
            self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{key}'")

    def copy(self):
        return DotAccessibleDict(super().copy())
    
import functools
import termcolor
import time
timeit_depth = 0
def timeit(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global timeit_depth
        timeit_depth += 1
        print(termcolor.colored(f"{'=== ' * timeit_depth}Starting {func.__name__}", "green"))
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        print(termcolor.colored(f"{'=== ' * timeit_depth}Elapsed time for {func.__name__}: {elapsed:.2f} seconds", "red"))
        timeit_depth -= 1
        return result
    return wrapper
