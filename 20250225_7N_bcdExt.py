import regex
import gzip
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

import scipy
import numpy as np
from scipy.sparse import dok_array, save_npz
from scipy.sparse import csr_matrix



dataDir = "" 
dataPrefix = "AM-LF01-"
fastqList = ['D769','D770','D771']
             
dataPostfix = "_R1.fastq.gz"

fastqnames = [dataDir + dataPrefix + fastqList[i] + dataPostfix for i in range(len(fastqList))]


file = "Plate32_bcds_canavanine"
bcds = np.loadtxt(file, dtype = 'str')


def encodeATGC(bcdStr):
    'Associates neucleotides into numeric IDs'
    return int(bcdStr.replace("A","0").replace("T","1").replace("G","2").replace("C","3"), 4)


def encodeATGC_revComp(bcdStr):
    'Returns the encoded reverse complement'
    bcdStr = bcdStr[::-1]
    return int(bcdStr.replace("A","1").replace("T","0").replace("G","3").replace("C","2"), 4)

def decodeATGC(bcdInt, bcdLen):
    'When given the numeric ID of the strand, returns the neucleotide names'
    bcdString = np.base_repr(bcdInt, base = 4)
    bcdString = ("0"*(bcdLen-len(bcdString)) ) + bcdString
    return bcdString.replace("0","A").replace("1","T").replace("2","G").replace("3","C")




umilength = 7
bcdLength = 20

# UF Primer Regex
UMI = "([A|T|G|C]{" + str(umilength) + "})"
primer = "(CCACGAGGTCTCT){e<=2}"
barcodes = [f"({bcds[i]})"+'{e<=1}' for i in range(len(bcds))]
threeprime = "(CGTACGCTGCAGGT){e<=2}"


regexObjects = [regex.compile(UMI + primer + barcodes[i] + threeprime) for i in range(len(bcds)) ]


'Extracting barcode counts from the fastq data'

bcdCounts = np.zeros((len(fastqnames),len(bcds)))
encoded_bcds = [encodeATGC(i) for i in bcds]

# Open the FASTQ file
for j in range(len(fastqnames)):
    print("Starting " + fastqList[j])
    counter = 0
    qc = 0 
    with gzip.open(fastqnames[j],"rt") as fq: 
            for record in SeqIO.parse(fq, "fastq"):

                
                for i in range(len(bcds)):
                    match = regexObjects[i].match(str(record.seq))
                    if match:
                        bcdCounts[j,i] += 1
                        if "N" not in match.group():
                            qc += 1
                         
                        
                        
    bcds_with_counts = [encoded_bcds, bcdCounts[j]] #2D list with barcodes in row 0 and their count in row 1\n",

    totalReads = np.sum(bcdCounts[j])                    
    print(fastqList[j])
    print(f'Total reads in fastq: {counter}')
    print(f"Total reads that match regex: {totalReads}")
    print(f'Reads passing quality control: {qc}')
    #print(f'i.e. {100 * qc/totalReads}')
    #print(f'Percentage of reads that match with a regex: {100 * totalReads/counter}')
    
    file_name = str(fastqList[j]) + 'bcd.npz'
    
    
    save_npz(file_name, csr_matrix(bcds_with_counts))


'Plotting barcode abundances in each fastq'
    
fig, axs = plt.subplots(1, 1, figsize=(30,5))

with np.errstate(divide='ignore'):
    plot = axs.imshow(np.log10(bcdCounts))
axs.set_aspect(aspect = 5)
print('Percentage of nonzero barcodes: ')
print(100*len(np.nonzero(bcdCounts[0])[0])/len(bcdCounts[0]))
plt.tight_layout()
fig.colorbar(plot)
plt.title('Barcode Abundances')
plt.show()
plt.savefig("7N-UF" + "bcd.pdf")   
