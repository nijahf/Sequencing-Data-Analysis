
#### Data files


dataDir = " insert/directory/ " 
dataPrefix = "da-ta-pre-fix"
fastqList = ['fastq','q','file','indenfier']
             
dataPostfix = "da-ta-post-fix.gz"

F2primers = False # Is the data with the F2 + R2-2 primers? If False, assumed to be F1 + R1-2


#### Code

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



bcdlength = 12
umilength_1 = 8

##### Specifying Regular expressions

# F1 + R1-1 Primer Regex
UMI = "([A|T|G|C]{" + str(umilength_1) + "})"
primer = "(GGATCCGATCATGCTT){e<=3}"
barcode = "([A|T|G|C]{4}(TT){e<=1}[A|T|G|C]{4}(TT){e<=1}[A|T|G|C]{4})"
threeprime = "(GGTACCGCTGATTAGT){e<=3}"

read1regex_F1 = regex.compile(UMI + primer + barcode + threeprime)


# F2 + R2-2 Primer regex
UMI = "([A|T|G|C]{" + str(umilength_1) + "})"
primer = "(GGATAAAATGTGATAACTAATCAGCGGTACC){e<=4}"
barcode = "([A|T|G|C]{4}(AA){e<=1}[A|T|G|C]{4}(AA){e<=1}[A|T|G|C]{4})"
threeprime = "(AAGCATGATCGGATCC){e<=3}"

read1regex_F2 = regex.compile(UMI + primer + barcode + threeprime)


if F2primers:
    encodeBcd = encodeATGC_revComp
    read1regex = read1regex_F2
else:
    encodeBcd = encodeATGC
    read1regex = read1regex_F1


for fastq in fastqList:
    print("Starting " + fastq)
    fname = dataPrefix + fastq + dataPostfix

    bcdUMICounts = dok_array((4**bcdlength, 4**(umilength_1)), dtype=np.uint64)
    bcdCounts = dok_array((4**bcdlength, 1), dtype=np.uint64)

    allseq = 0
    bcdseq = 0

    with gzip.open(dataDir + fname, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            allseq = allseq + 1
            m = read1regex.match(str(record.seq))
            if m:
                bcd = m.groups()[2][0:4] + m.groups()[2][6:10] + m.groups()[2][12:16] 
                if "N" not in bcd:
                    bcdseq = bcdseq + 1
                    umi = m.groups()[0]
                
                    bcdUMICounts[encodeBcd(bcd), encodeATGC(umi)] = bcdUMICounts[encodeBcd(bcd), encodeATGC(umi)] + 1
                    bcdCounts[encodeBcd(bcd), 0] = bcdCounts[encodeBcd(bcd), 0] + 1

                
    print(fastq)            
    print("Total reads: " + str(allseq))
    print("Reads passing quality control: " + str(bcdseq))
    print("i.e. " + str(100*bcdseq/allseq) + "%")

    save_npz(fastq.split("S")[0] + "bcdUMI.npz",bcdUMICounts.tocsr())
    save_npz(fastq.split("S")[0] + "bcd.npz",bcdCounts.tocsr())
