### The purpose of this script is to give me information about RNA editing sites using SNP data (in the form of an intersect gff)
# with a gene annotation (.gff) file and fasta. Specifically, I want to know is it U-C or` C-U? next is it modifying a start codon or
# a stop codon or an AA? (eventually perhaps identify which amino acids it is changing). In another script I would like to combine a .vcf
# from RNAseq data to get a count file containing proportion of edited transcripts at a given time point.
# /Users/david/PycharmProjects/Isoetes/RNAedit_info.py /Users/david/Documents/Isoetes/RNA_editing/Isoetes_taiwanensis_cp_v1.gff /Users/david/Documents/Isoetes/RNA_editing/Isoetes_taiwanensis_gene.fasta /Users/david/Documents/Isoetes/RNA_editing/CUedit_intersect_gene.out

import sys
from Bio import SeqIO
from collections import Counter
import operator

inFile1 = open(sys.argv[1], 'r')  # gff
inFile2 = open(sys.argv[2], 'r')  # fasta
inFile3 = open(sys.argv[3], 'r')  # intersect

outFile = open(sys.argv[3]+".vcf", 'w')  # new vcf
outFile2 = open(sys.argv[3]+"nonsynonymous.vcf", 'w')  # vcf containing nonsynonymous edits
outFile3 = open(sys.argv[3]+"startCodon.vcf", 'w')  # vcf containing nonsynonymous edits

# First create a simple codon dictionary (  adapted from http://www.petercollingridge.co.uk/tutorials/bioinformatics/codon-table/)
bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

geneCt = 0
with inFile1 as f:
    gffDict = {}
    lenDict = {}
    for line in f:
        if line.startswith('#'):
            pass
        else:
            spLine = line.strip().split('\t')
            if spLine[2] == 'gene':  # this could easily be modified to allow the user to specify a tag e.g. gene CDS mRNA etc...
                geneCt += 1
                ID = spLine[8].split(";")[0].lstrip('ID=')
                Name = spLine[8].split(";")[1].lstrip('Name=')# ditto see above
                lenDict[Name] = int(spLine[4]) - int(spLine[3])
                gffDict[ID] = [Name, spLine[3], spLine[4], spLine[6]]  # this gets gene coordinates and strandedness though it has occured to me I could just take the whole spLine as a list
            else:
                pass
# print(gffDict)

# all gene sequences are printed in the + orientation using the -s flag in:
# bedtools getfasta -fi Isoetes_taiwanensis_cp.fasta -bed Isoetes_taiwanensis_cp_v1-geneOnly.bed -fo Isoetes_taiwanensis_gene.fasta -name -s
with inFile2 as f:
    fastaDict = {}
    seqiter = SeqIO.parse(f, 'fasta')
    for record in seqiter:
        seqID = record.id[:-3]  # cut off last 3 characters
        if seqID in gffDict:
            fastaDict[seqID] = list(str(record.seq))  # I had to use str here to get the sequence to print out on its own otherwise I get the entire record from seqIter

# print(fastaDict)                                                  # Also converted string to a list to allow substitution
# Set up some counters to get a few basic statistics:
pos1Ct = 0  # number of edits at codon position 1
pos2Ct = 0
pos3Ct = 0
startCt = 0  # number of edits in the first three nucleotides of a gene
stopCt = 0
CUedit = 0
UCedit = 0
twoFer = 0  # number of codons with 2 editing sites
silentCt = 0
NameList = []  # list of gene Names to calculate edits per gene
twoFerList = []  # genes containing codons with 2 editing sites
UCList = []  # genes with U-C edits


# bedtools intersect -a Isoetes_taiwanensis_cp_v1.gff -b KBTI_Itai-cpFILT-MQB0.1QUAL300_forIntersect.vcf -wo | grep -P "gene\t" - > CUedit_intersect_gene.out
with inFile3 as f:
    prevPos = 0
    for line in f:
        spLine = line.strip().split("\t")
        cuID = spLine[8].split(";")[0].lstrip('ID=')  # get ID and name as above
        Name = spLine[8].split(";")[1].lstrip('Name=')
        gnmPos = int(spLine[10])
        NameList.append(Name)
        editBySampleList = spLine[18:-1]  # get vcf data from each edit
        if cuID in gffDict:
            outFile.write(line)
            if gffDict[cuID][3] == '+':
                genePos = gnmPos - int(gffDict[cuID][1])
                initPos = spLine[12]
                editPos = spLine[13]
            elif gffDict[cuID][3] == '-':
                genePos = int(gffDict[cuID][2]) - gnmPos
                if spLine [12] == 'G':
                    initPos = 'C'
                    editPos = 'T'
                elif spLine [12] == 'A':
                    initPos = 'T'
                    editPos = 'C'
        else:
            print(cuID)
        if genePos == 1:
            startCt += 1
            print("{} to {} edit in {}: {}_{} creates a START codon".format(initPos, editPos, Name, cuID, gnmPos))
            edit_list = fastaDict[cuID][genePos - 1:genePos + 2]
            init_codon = ''.join(edit_list)
            init_AA = codon_table[init_codon]
            edit_codon = ''.join([edit_list[0], str(editPos), edit_list[2]])
            print("{} : {}".format(init_codon, init_AA))
            print("{} : {}".format(edit_codon, edit_AA))
            outFile3.write(line)
        else:
            if genePos % 3 == 1:
                pos1Ct += 1
                print("{} to {} edit in {}: {}_{} is at the first position".format(initPos, editPos, Name, cuID, gnmPos))
                edit_list = fastaDict[cuID][genePos-1:genePos + 2]
                init_codon = ''.join(edit_list)
                init_AA = codon_table[init_codon]
                edit_codon = ''.join([str(editPos), edit_list[1], edit_list[2]])
                edit_AA = codon_table[edit_codon]
                print("{} : {}".format(init_codon, init_AA))
                print("{} : {}".format(edit_codon, edit_AA))
            elif genePos % 3 == 2:
                pos2Ct += 1
                print("{} to {} edit in {}: {}_{} is at the second position".format(initPos, editPos, Name, cuID, gnmPos))
                edit_list = fastaDict[cuID][genePos - 1:genePos + 2]
                init_codon = ''.join(edit_list)
                init_AA = codon_table[init_codon]
                edit_codon = ''.join([edit_list[0], str(editPos), edit_list[2]])
                edit_AA = codon_table[edit_codon]
                print("{} : {}".format(init_codon, init_AA))
                print("{} : {}".format(edit_codon, edit_AA))
                if gnmPos - prevPos <= 2:
                    twoFer += 1
                    print("TWO IN ONE!!!")
                    twoFerList.append("{} {}_{}\t{} -> {}".format(Name, cuID, gnmPos, init_AA, edit_AA))
            else:
                pos3Ct += 1
                print("{} to {} edit in {}: {}_{} is at the third position".format(initPos, editPos, Name, cuID, gnmPos))
                edit_list = fastaDict[cuID][genePos - 2:genePos + 1]
                init_codon = ''.join(edit_list)
                init_AA = codon_table[init_codon]
                edit_codon = ''.join([edit_list[0], edit_list[1], str(editPos)])
                edit_AA = codon_table[edit_codon]
                print("{} : {}".format(init_codon, init_AA))
                print("{} : {}".format(edit_codon, edit_AA))
                if gnmPos - prevPos <= 2:
                    twoFer += 1
                    print("TWO IN ONE!!!")
                    twoFerList.append("{} {}_{}\t{} -> {}".format(Name, cuID, gnmPos, init_AA, edit_AA))
        if edit_AA == "*" or init_AA == "*":
            stopCt += 1
        if initPos == "C":
            CUedit += 1
        if initPos == "T":
            UCedit += 1
            UCList.append("{}_{}\t{} -> {}".format(cuID, gnmPos, init_AA, edit_AA))
        if init_AA == edit_AA:
            silentCt += 1
        if init_AA != edit_AA:
            outFile2.write(line)
        prevPos = gnmPos
        # print(editBySampleList)
editPerGene = Counter(NameList)
sorted_editPerGene = dict(sorted(editPerGene.items(), key=operator.itemgetter(1), reverse=True))

editPropDict = {}
for k in sorted_editPerGene:
    editProp = sorted_editPerGene[k]/lenDict[k]
    editPropDict[k] = editProp

sorted_editPropDict = dict(sorted(editPropDict.items(), key=lambda item: item[1], reverse=True))

print("{} edits in {} genes".format(UCedit+CUedit, geneCt))

print("{} C to U\n{} U to C\n{} edits at 1st position\n{} edits at 2nd position\n{} edits at 3rd position\n{} edits modify a start codon\n{} edits modify a stop codon\n{} are silent (do not alter amino acid)".format(CUedit, UCedit, pos1Ct, pos2Ct, pos3Ct, startCt, stopCt, silentCt))
for k in sorted_editPerGene:
    print("{} : {},".format(k, sorted_editPerGene[k]))
print("{} codons with 2 sites edited".format(twoFer))
for i in twoFerList:
    print(i)
print("U-C edit locations")
for i in UCList:
    print(i)


for k in sorted_editPropDict:
    print("{} : {},".format(k, sorted_editPropDict[k]))

print(len(set(NameList)))
