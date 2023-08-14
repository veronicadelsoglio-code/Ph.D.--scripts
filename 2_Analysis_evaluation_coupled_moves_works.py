#!/usr/bin/env python
# coding: utf-8

import Bio.PDB
import os
import weblogo
import numpy as np
import pathlib
import glob
import pandas as pd
import seaborn as sns
import matplotlib as plt

# current path
pwd=pathlib.Path().absolute()
print (pwd)


# # converts PDB to FASTA

def PDB2FASTA(pdb_file):
    with open(pdb_file) as input_file:
        output_file= open(pdb_file.split('.')[-2]+'.fasta', 'w')

        aa = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V', 'BCY':' BCY '}

        output_file.write(">"+pdb_file.split('/')[-1]+"\n")

        previous=0
        for line in input_file:
            if line[0:4] == 'ATOM':
                if line [21] == 'A':
                    new=aa[line[17:20]]
                    resnumber=line[23:26]
                    if resnumber != previous:
                        output_file.write(new)
                        previous=resnumber
        output_file.write('\n')
        input_file.close()
        output_file.close()


# creates WEBLOGO 3
# Sets up WebLogo3 using a reference (original) sequence, which gets colored in orange

class RefSeqColor(weblogo.ColorScheme):
    """
    Color the given reference sequence in its own color, so you can easily see 
    which positions match that sequence and which don't.
    """


    def __init__(self, ref_seq, color, description=None):
        self.ref_seq = ref_seq
        self.color = weblogo.Color.from_string(color)
        self.description = description

    def symbol_color(self, seq_index, symbol, rank):
        if symbol == self.ref_seq[seq_index]:
            return self.color

baserules = [
            weblogo.SymbolColor("GSTYC", "green", "polar"),
            weblogo.SymbolColor("NQ", "purple", "neutral"),
            weblogo.SymbolColor("KRH", "blue", "basic"),
            weblogo.SymbolColor("DE", "red", "acidic"),
            weblogo.SymbolColor("PAWFLIMV", "black", "hydrophobic")
        ]

protein_alphabet=weblogo.Alphabet('ACDEFGHIKLMNPQRSTVUWY', zip('acdefghiklmnpqrstvuwy','ACDEFGHIKLMNPQRSTVUWY'))

def plotseqlogo(refseq, fasta_file, name):
    seqs = weblogo.read_seq_data(fasta_file, alphabet=protein_alphabet)

    colorscheme = weblogo.ColorScheme([RefSeqColor(refseq, "orange", "refseq")] + baserules,
                                alphabet = protein_alphabet)

    data = weblogo.LogoData.from_seqs(seqs)
    options = weblogo.LogoOptions()
    options.logo_title = name
    options.unit_name="probability"
    options.show_fineprint = True
    options.yaxis_label = ""
    options.color_scheme = colorscheme
    ###
    #eps_binary = weblogolib.eps_formatter( data, logo_format )
    #eps_str = eps_binary.decode()
    #eps_str = replace_logo_numbers(eps_str, self.design_positions, self.starting_seq )
    
    #with open( eps_logo_filename, 'w') as f:
        #f.write( eps_str )
    #eps_to_pdf(eps_logo_filename)
    #####
    mformat = weblogo.LogoFormat(data, options)

    fname = "%s_energy.pdf" % name
    with open(fname, "wb") as f:
        f.write(weblogo.pdf_formatter(data, mformat))

    fname = "%s_energy.png" % name
    with open(fname, "wb") as f:
        f.write(weblogo.png_print_formatter(data, mformat))
 
# adapted from source: https://baoilleach.blogspot.com/2017/05/using-weblogo3-to-create-sequence-logo.html


#creates fasta from pdb files

filepath = str(pwd)+"/C*_S*_3*_0*_0*_0*pdb"   
pdbfiles = glob.glob(filepath)
for file in pdbfiles:
    PDB2FASTA(file)
    


# concenates all fasta files to all.fasta'

filepath = str(pwd)+"/C*_S*_3*_0*_0*_0*fasta"   
fasta = glob.glob(filepath)
with open('all.fasta', 'w') as outfile:
    for file in fasta:
        with open(file) as infile:
                for line in infile:
                    outfile.write(line)



# get positions which were designed (resid)
# creates original.fasta containing non design sequence
#filepath = str(pwd)+"/*out"
#for outfile in glob.glob(filepath):
#    with open(outfile) as log_file:
#        design_resid=set()
#        found="false" #to stop the search after the first found instance
#        for line in log_file:
# this line contains design positions as numbers and original sequenceat the end
#            if found=="false" and line.startswith("protocols.CoupledMovesProtocol: Design Positions:"):
#                found="true"
                #with open('original.fasta', 'w') as outfile:
                #    outfile.write(">original"+"\n"+str(line.split(' ')[-1]))
                #    outfile.close()
                #original_seq= list(line.split(' ')[-1])[0:-1] # to remove "\n" at the end of sequence
#                for i in line.split(' '):
#                    if i.isdigit() == True:
#                        design_resid.add(int(i))

#design_resid=sorted(list(design_resid))
#print (design_resid)

from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)




## name of native pdb ahs to be here:
structure = parser.get_structure("native", ''.join(glob.glob("C*_S*_3*_0*_0*_N*pdb")))
lenght_protein=(len(structure[0]["A"]))
print (len(structure[0]["A"]))
#lenght_protein=1
#for model in structure[0]["B"]:
#    for residue in model.get_residues():
#        if Bio.PDB.is_aa(residue):
#            lenght_protein += 1

print (lenght_protein)
            
original_seq3=[]
for res in range(1,lenght_protein+1):
    residue = structure[0]["A"][res]
    original_seq3.append(residue.get_resname())

    


print (original_seq3)

from Bio.SeqUtils import seq1
original_seq=[]
for i in original_seq3:
    original_seq.append(seq1(i).upper())
print (original_seq)


# # calculates PSSM

import Bio.AlignIO
import Bio.SubsMat
from Bio.Align import AlignInfo

align = Bio.AlignIO.read("all.fasta", "fasta")
summary_align = AlignInfo.SummaryInfo(align)
my_pssm = summary_align.pos_specific_score_matrix('')

#print (my_pssm[2])


def substitutions(mutations, cutoff=0.0001):
    allowed_substitutions={}
    min_pssm={}
    for i in range(0,len(mutations)):
        total = sum(my_pssm[i].values(), 0.0)
        norm_pssm = {k: v / total for k, v in my_pssm[i].items()}
        #print (norm_pssm)

        min_pssm[mutations[i]]={}
        for (k, v) in norm_pssm.items():
            if v >= cutoff:
                min_pssm[mutations[i]][k]= v
                
    return (min_pssm)
            
pssm_table=substitutions(range(1,lenght_protein+1))                
print (pssm_table)


all_fasta = open("all.fasta", "r")
name="allsequences"
print (len(original_seq))
print (len(pssm_table))

plotseqlogo(original_seq, all_fasta, name)


# <img src="allsequences_energy.png" width=1000 height=80 />

## Gets energyterm (default=total) of a designed residue for a residue and a pdb file

def get_energy(residuenumber, pdb_file, energy_term="total", ligand="none"):
    with open(pdb_file) as pdb_file:
        energy_parts=[]
        weights=[]
        for line in pdb_file:
            #if ligand != "none":
            #if line.startswith("RNI"):
            #    return ("Ligand found")
                    #term=energy_parts.index(energy_term)
                    #return (str(line.rsplit()[0])+" "+str(line.rsplit()[term]))
            #print (energy_parts)
            if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                #print (next(pdb_file))
                for i in next(pdb_file).rsplit():
                    energy_parts.append(i)
                for i in next(pdb_file).rsplit():
                    weights.append(i)
            if "_"+str(residuenumber)+" " in line:
                term=energy_parts.index(energy_term)
                return (str((line.rsplit()[0])).split("_")[0]+" "+str((line.rsplit()[0])).split("_")[-1]+" "+str(line.rsplit()[term]))

def get_total_energy(pdb_file, energy_term="total"):
    with open(pdb_file) as pdb_file:
        energy_parts=[]
        weights=[]
        for line in pdb_file:
 
            if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                #print (next(pdb_file))
                for i in next(pdb_file).rsplit():
                    energy_parts.append(i)
                for i in next(pdb_file).rsplit():
                    weights.append(i)
            if "pose"+" " in line:
                term=energy_parts.index(energy_term)
                return ((pdb_file.name).split("/")[-1]+" "+str(str(line.rsplit()[term])))
            
def get_total_link_res(pdb_file):
    with open(pdb_file) as pdb_file:
        for line in pdb_file:
 
            if line.startswith("LINK"):
                return (line.split()[4])



# concenates all total energies

filepath_pdb = str(pwd)+"/C*_S*_3*_0*_0*_0*pdb"
pdbfiles = glob.glob(filepath_pdb)
outfile = open("total_energy.energy", "a+")
for file in pdbfiles: 
    outfile.write(get_total_energy(file)+"\n")
outfile.close()


# concenates all designed residue energies


outfile = open("design_res.energy", "a+")
for file in pdbfiles: 
    for res in range(1,lenght_protein+1):    
        outfile.write(get_energy(res, file)+"\n")
outfile.close()

# energies of every single designed residue 



for file in pdbfiles:    
    for i in range(1,lenght_protein+1):    
        outfile = open(str(i)+".energy", "a+")
        outfile.write(str(get_energy(i, file, ligand="RNI"))+"\n")
        outfile.close()
    #else:
    #    continue




energy_res=[] # lists all residues were more than one amino acids is possible after the first design round

for i in range(1,lenght_protein+1):
    
    df=pd.read_csv(str(i)+".energy", sep=' ', header=None, names=["filename","res","res_number","energy"])
    grouped = (df.groupby(df.res)).nunique()
    data=df.groupby('res').agg([np.mean, np.std])
    if len(grouped)>0:
        data=df.groupby('res').agg([np.mean, np.std])
        #print (data)
        #print (df.quantile(q=0.25, axis=0, numeric_only=True, interpolation='linear'))



df=pd.read_csv("design_res.energy", sep=' ', header=None, names=["res","res_number","energy"])
#print(df)
sns.set_theme(style="ticks", color_codes=True)
g=sns.catplot(x="res", y="energy", col="res_number",col_wrap=5, data=df, kind="strip", sharey=False, sharex=False)
(g.set_axis_labels("", "res")).set(ylim=(-6, 2))
g = g.map(lambda y, **kw: plt.pyplot.axhline(y.quantile(0.25), color="k"), 'energy')
  
g.savefig("energies_per_res.png")



df=pd.read_csv("design_res.energy", sep=' ', header=None, names=["res","res_number","energy"])
#print(df)

