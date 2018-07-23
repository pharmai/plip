"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
mp.py - Functions for parallel processing
"""

# Python Standard Library
from __future__ import division
from builtins import zip
import multiprocessing
import itertools
from numpy import asarray
from functools import partial


class SubProcessError(Exception):
    def __init__(self, e, exitcode=1):
        self.exitcode = exitcode
        super(SubProcessError, self).__init__(e)
    pass


def universal_worker(input_pair):
    """This is a wrapper function expecting a tiplet of function, single
       argument, dict of keyword arguments. The provided function is called
       with the appropriate arguments."""
    function, arg, kwargs = input_pair
    return function(arg, **kwargs)


def pool_args(function, sequence, kwargs):
    """Return a single iterator of n elements of lists of length 3, given a sequence of len n."""
    return zip(itertools.repeat(function), sequence, itertools.repeat(kwargs))


def parallel_fn(f):
    """Simple wrapper function, returning a parallel version of the given function f.
       The function f must have one argument and may have an arbitray number of
       keyword arguments. """

    def simple_parallel(func, sequence, **args):
        """ f takes an element of sequence as input and the keyword args in **args"""
        if 'processes' in args:
            processes = args.get('processes')
            del args['processes']
        else:
            processes = multiprocessing.cpu_count()

        pool = multiprocessing.Pool(processes)  # depends on available cores

        result = pool.map_async(universal_worker, pool_args(func, sequence, args))
        pool.close()
        pool.join()
        cleaned = [x for x in result.get() if x is not None]  # getting results
        cleaned = asarray(cleaned)
        return cleaned
    return partial(simple_parallel, f)
