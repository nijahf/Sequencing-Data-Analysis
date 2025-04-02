
import regex
import gzip
from Bio import SeqIO

import scipy
import numpy as np
from scipy.sparse import dok_array, save_npz



def encodeATGC(bcdStr):
    return int(bcdStr.replace("A","0").replace("T","1").replace("G","2").replace("C","3"), 4)


def encodeATGC_revComp(bcdStr):
    bcdStr = bcdStr[::-1]
    return int(bcdStr.replace("A","1").replace("T","0").replace("G","3").replace("C","2"), 4)

def decodeATGC(bcdInt, bcdLength = 12):
    bcdString = np.base_repr(bcdInt, base = 4)
    bcdString = ("0"*(bcdLength-len(bcdString)) ) + bcdString
    return bcdString.replace("0","A").replace("1","T").replace("2","G").replace("3","C")
