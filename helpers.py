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

import sys
def assertx(cond, msg=None, show_stack=False):
    """
    Like `assert`, but the third argument controls whether an unhandled
    AssertionError prints a traceback.
    Usage: assertx(a == b, "a not equal b", False)
    """
    import sys

    if cond:
        return

    if show_stack:
        # behave like usual - raise an AssertionError (will show traceback if uncaught)
        raise AssertionError(msg)
    else:
        # print a single-line error and exit without printing a traceback
        sys.stderr.write(f"\nError: {msg}\n\n")
        sys.exit(1)
        
def plot_binomial_samples(N, p, samples=100000):
    from scipy.stats import binom
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import skew

    assertx(isinstance(N, int) and not isinstance(N, bool) and 1 <= N <= 10000,
            "N must be an integer between 1 and 100000")
    
    # draw draws from the binomial and plot histogram + theoretical pmf
    draws = binom.rvs(n=N, p=p, size=samples)

    mean = np.mean(draws)
    stddev = np.std(draws, ddof=1)  # sample standard deviation
    #var = np.var(draws, ddof=1)
    #skewness = skew(draws, bias=False)
    print(f"mean: {mean:.3f}\t\tshould be ~Np: {N*p:.3f}")
    print(f"stddev/N: {stddev/N:.5f}\tshould be ~sqrt(p(1-p)/N): {np.sqrt(p*(1-p)/N):.5f}")
    #print(f"skewness: {skewness:.5f}")

    plt.figure(figsize=(10, 4))
    bins = np.arange(0, N + 2) - 0.5
    plt.hist(draws, bins=bins, density=True, alpha=0.6, label=f"1 sample" if samples==1 else f"{samples} samples")

    k = np.arange(0, N + 1)
    pmf = binom.pmf(k, N, p)
    plt.plot(k, pmf, 'r.', markersize=3, label='theoretical pmf')

    (xmin, xmax) = (-0.5, N+0.5) if N<100 else (0, N)
    plt.xlim(xmin,#max(0, int(N * p - 4 * np.sqrt(N * p * (1 - p)))), 
             xmax)#min(N, int(N * p + 4 * np.sqrt(N * p * (1 - p)))))
 
    
    # plot vertical lines for + and - stddev around the mean
    #plt.axvline(mean, color='k', linestyle='--', label='mean')
    plt.axvline(mean - stddev, color='gray', linestyle=':', label='mean +/- stddev')
    plt.axvline(mean + stddev, color='gray', linestyle=':')
    
    # xticks should be integers with reasonable spacing
    xmin, xmax = plt.xlim()
    start = max(0, int(np.floor(xmin)))
    end = int(np.ceil(xmax))
    max_ticks = 20
    step = max(1, int(np.ceil((end - start) / max_ticks)))
    ticks = np.arange(start, end, step)
    plt.xticks(ticks, [str(int(t)) for t in ticks])

    # # set y-axis to log scale and ensure a small positive lower bound to avoid zeros
    # ax = plt.gca()
    # ax.set_yscale('log')
    # positive = pmf[pmf > 0]
    # ymin = float(positive.min()) if positive.size else 1e-12
    # ax.set_ylim(bottom=max(ymin * 1e-3, 1e-12))
    
    plt.xlabel('Number of successes (offspring)')
    plt.ylabel('Probability')
    plt.title(f'Binomial(n={N}, p={p}) samples vs theoretical PMF')
    plt.legend()
    plt.show()
    
