"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
mp.py - Functions for parallel processing
Copyright 2014-2015 Sebastian Salentin, Joachim Haupt

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# Python Standard Library
from __future__ import division
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
    return itertools.izip(itertools.repeat(function), sequence, itertools.repeat(kwargs))


def parallel_fn(f):
    """Simple wrapper function, returning a parallel version of the given function f.
       The function f must have one argument and may have an arbitray number of
       keyword arguments. """
       
    def simple_parallel(func, sequence, **args):
        """ f takes an element of sequence as input and the keyword args in **args"""
        multiprocessing.freeze_support()
        if 'processes' in args:
            processes = args.get('processes')
            del args['processes']
        else:
            processes = multiprocessing.cpu_count()
        
        pool = multiprocessing.Pool(processes=processes)  # depends on available cores
        result = pool.map(universal_worker, pool_args(func, sequence, args))
        cleaned = [x for x in result if x is not None]  # getting results
        cleaned = asarray(cleaned)
        pool.close()
        pool.join()
        return cleaned
    return partial(simple_parallel, f)

