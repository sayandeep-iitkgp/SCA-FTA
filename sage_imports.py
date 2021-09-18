r"""
Sage functions and imports for enabling GF(2) operations.

AUTHOR: Sayandeep Saha		
"""

from random import randint
import fileinput
import collections
import os
import math
import random as pyrandom

import commands, re

from sage.misc.misc import alarm, walltime

from sage.rings.all import FiniteField as GF
from sage.rings.integer_ring import ZZ

from sage.crypto.mq.sbox import SBox
from sage.structure.sequence import Sequence
from sage.misc.prandom import random

from sage.combinat.permutation import Permutation

from multiprocessing import Process, Queue
from Queue import Empty

from sage.misc.misc import *




def bitstringtohexstring(l):
    """
    Return a hex string in PRESENT style for l a list of bits.

    INPUT:
        l -- a list with bit entries of length divisible by 4
    """
    r = []
    for i in xrange(0,len(l),4):
        z = list(reversed(map(int, l[i:i+4])))
        r.append(hex(ZZ(z,2)))

    r = sum([r[i:i+8]+[" "] for i in xrange(0,len(r),8) ],[])

    return "".join(r)[:-1]

def hexstringtobitstring(n, length=64):
    """
    Return a hex string in PRESENT style for l a list of bits.

    INPUT:
        l -- a list with bit entries of length divisible by 4
    """
    n = int(n,16)
    l = []
    for i in xrange(length):
        l.append(1 if 2**i & n else 0)
    l = map(GF(2),l)
    return list(reversed(l))

def intarraytohexstring(l):
	"""
	Return a hex string for l a list of integers in [0,15].

	INPUT:
		l -- a list with bit entries of length divisible by 4
	"""	
	hexstr = ''.join('{:0x}'.format(x) for x in l)
	return hexstr

def bitstringtointarray(l):
	return [int(t, 16) for t in bitstringtohexstring(l).replace(" ", "")]	
