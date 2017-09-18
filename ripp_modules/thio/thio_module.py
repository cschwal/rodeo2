#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 17:09:04 2017

@author: bryce
"""
import csv
import os
import re
import numpy as np
from ripp_modules.thio.svm import svm_classify as svm
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.Virtual_Ripp import Virtual_Ripp
import hmmer_utils

index = 0

def write_csv_headers(output_filename):
    dir_prefix = 'output/thio/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification,Contains TOMM YcaO PF02624,Contains LanB Nterm PF04738,Contains LanB Cterm PF14028,Contains TOMM dehy PF00881,Contains rSAM MTase PF04055,Contains P450 PF00067,Contains ABC trans1 PF00005,Contains ABC trans2 PF01061,Contains ABC trans3 PF12698,Contains abhydrolase1 PF12697,Contains abhydrolase2 PF00561,CSS motif,CTT motif,SS motif,SSS motif,SSSS motif,CC motif,CCC motif,CCCC motif,TT motif,TTT motif,TTTT motif,No Cys core residues,No Ser core residues,No Thr core residues,Core mass < 2100,Sum of repeating Cys/Ser/Thr > 4,Avg heterocycle block length > 3,Leader net charge < 5,Leader net charge > 0,Leader contains a Cys?,Peptide terminates Cys/Ser/Thr,Core contains >= 2 positive residues,Heterocycle ratio > 0.4,Number of core repeating blocks,Number of core repeating Cys,Number of core repeating Ser,Number of core repeating Thr,Number of core heterocycle blocks,avg core heterocycle block length,Precursor peptide mass (unmodified),Leader peptide mass (unmodified),Core peptide mass (unmodified),Length of Precursor,Length of Leader,Length of Core,Leader / core ratio,Heterocycle residues/lenth of core,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl'
    svm_headers = svm_headers.split(',')
    features_headers = ["Accession_id", "Genus/Species/Code", "Leader", "Core", "Start", "End", "Total Score", "Valid Precursor" ] + svm_headers
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("ripp_modules/thio/svm/fitting_set.csv", 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(features_headers)
    svm_writer.writerow(svm_headers)

def ripp_write_rows(output_filename, accession_id, genus_species, list_of_rows):
    dir_prefix = 'output/thio/'
    global index
    features_csv_file = open(dir_prefix + "temp_features.csv", 'a')
    svm_csv_file = open("ripp_modules/thio/svm/fitting_set.csv", 'a')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    for row in list_of_rows:
        features_writer.writerow([accession_id, genus_species] + row[0:5] + ["valid_precursor_placeholder", index, ''] + row[5:])
        svm_writer.writerow([index, ''] + row[5:]) #Don't include accession_id, leader, core sequence, start, end, or score
        index += 1
        
def run_svm():
    svm.run_svm()
    svm_output_reader = csv.reader(open("ripp_modules/thio/svm/fitting_results.csv"))
    final_output_writer = csv.writer(open("output/thio/thio_features.csv", 'w'))
    features_reader = csv.reader(open("output/thio/temp_features.csv"))
    header_row = features_reader.next() #skip header
    final_output_writer.writerow(header_row)
    for row in features_reader:
        svm_output = svm_output_reader.next()[1]
        row[9] = svm_output
        if int(svm_output) == 1:
            row[6] = int(row[6]) + 10
        if int(row[6]) > 20: #CUTOFF
            row[7] = 'Y'
        else:
            row[7] = 'N'
        final_output_writer.writerow(row) 

class Ripp(Virtual_Ripp):
    def __init__(self, 
                 start, 
                 end, 
                 sequence,
                 upstream_sequence,
                 pfam_2_coords):
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords)
        self.peptide_type = 'thio'
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        
    def set_split(self):
        """Try to identify cleavage site using regular expressions"""
        #Regular expressions; try 1 first, then 2, etc.
        rex1 = re.compile('([I|V]AS)')
        rex2 = re.compile('([G|A|S]AS)')
    
        #For each regular expression, check if there is a match that is <10 AA from the end
        if re.search(rex1,self.sequence) and len(re.split(rex1,self.sequence)[-1]) > 10:
            start, end = [m.span() for m in rex1.finditer(self.sequence)][-1]
            end -= 5
        elif re.search(rex2,self.sequence) and len(re.split(rex1,self.sequence)[-1]) > 10:
            start, end = [m.span() for m in rex2.finditer(self.sequence)][-1]
            end -= 5 #TODO +- 5? why
        else:
            self.split_index =  -1
            self.core = self.sequence
            self.leader = ''
            return
        self.split_index = end
        self.leader = self.sequence[:end]
        self.core = self.sequence[end:]
        
    def set_score(self):
        tabs = []
        score = 0
        pfams = []
        for pfam_dot in self.pfam_2_coords.keys():
            pfams.append(pfam_dot.split('.')[0])
        
        #Contains TOMM YcaO (PF02624)
        if "PF02624" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #TODO ask Chris why not just use pfam id...
        #Contains LanB N-terminal domain (PF04738)
        if "Lant_dehyd_N" in pfams or "PF04738" in pfams or "tsrC" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains LanB C-terminal domain (PF14028)
        if "PF14028" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains TOMM dehydrogenase (PF00881)
        if "PF00881" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains rSAM methyltransferase (PF04055)
        if "PF04055" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains P450 (PF00067)
        if "PF00067" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains ABC transporter or abhydrolase
        abc_transp_abhydrolases = ["PF00005", "PF01061", "PF12698", "PF12697", "PF00561"]
        for dom in abc_transp_abhydrolases:
            if dom in pfams:
                score += 1
                tabs.append(1)
            else:
                tabs.append(0)
        #CSS/CTT, SS/SSS/SSS, CC/CCC/CCCC, TT/TT/TTTT motifs
        motifs = (('[C][S]{2,}', 1), ('[C][T]{2,}', 1), ('[S]{2,}', 1), ('[S]{3,}', 1), \
           ('[S]{4,}', 2), ('[C]{2,}', 1), ('[C]{3,}', 1), ('[C]{4,}', 2), \
           ('[T]{2,}', 1), ('[T]{3,}', 1), ('[T]{4,}', 2))
        for motif in motifs:
            if re.search(motif[0], self.core) != None:
                score += motif[1]
                tabs.append(1)
            else:
                tabs.append(0)
        #No Cys/Ser/Thr core residues
        for aa in "CST":
            if aa not in self.core:
                score -= 2
                tabs.append(1)
            else:
                tabs.append(0)    
        #Mass of core peptide (unmodified) < 2100
        if "X" in self.core:
            Xs = self.core.count("X") #TODO why?>??
            noXcore = self.core.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            core_analysis = ProteinAnalysis(noXcore, monoisotopic=True)
        else:
            core_analysis = ProteinAnalysis(self.core, monoisotopic=True)
        if core_analysis.molecular_weight() < 2100:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Sum of repeating Cys/Ser/Thr > 4
        number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks = thioscout(self.core)
        if sum([int(nr) for nr in number_of_repeating_CST.split(", ")]) > 4:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Avg heterocycle block length > 3
        if avg_heteroblock_length != "nan" and avg_heteroblock_length > 3:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Leader net charge < 5
        charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
        leader_charge = sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])
        if leader_charge < 5:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Leader net charge > 0
        if leader_charge > 0:
            score -= 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Leader contains a Cys
        if "C" in self.leader:
            score -= 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Peptide terminates Cys/Ser/Thr
        if self.core[-1] in ["C", "S", "T"]:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Core contains >= 2 positive residues
        if sum([self.core.count(aa) for aa in "RK"]) >= 2:
            score -= 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Number of heterocyclizable residues to core ratio > 0.4
        if float(sum([self.core.count(aa) for aa in "CST"])) / len(self.core) >= 0.4:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
    
        #Now svm heuristics
        columns = []
        #append classification and index
        columns += tabs
        #Number repeating blocks of heterocyclizable residues in core
        number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks = thioscout(self.core)
        columns.append(number_of_repeat_blocks)
        #Number of core repeating Cys
        columns.append(int(number_of_repeating_CST.split(", ")[0]))
        #Number of core repeating Ser
        columns.append(int(number_of_repeating_CST.split(", ")[1]))
        #Number of core repeating Thr
        columns.append(int(number_of_repeating_CST.split(", ")[2]))
        #Number of blocks of heterocyclizable residues in core
        columns.append(number_of_heteroblocks)
        #Average core heterocycle block length
        if avg_heteroblock_length == "nan":
            columns.append(0)
        else:
            columns.append(avg_heteroblock_length)
        #Precursor peptide mass (unmodified)
        if "X" in self.sequence:
            Xs = self.sequence.count("X")
            noXprecursor = self.sequence.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            precursor_analysis = ProteinAnalysis(noXprecursor, monoisotopic=True)
            columns.append(float(precursor_analysis.molecular_weight()) + 110 * Xs)
        else:
            precursor_analysis = ProteinAnalysis(self.sequence, monoisotopic=True)
            columns.append(float(precursor_analysis.molecular_weight()))
        #Unmodified leader peptide mass
        if "X" in self.leader:
            Xs = self.leader.count("X")
            noXleader = self.leader.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            leader_analysis = ProteinAnalysis(noXleader, monoisotopic=True)
            columns.append(float(leader_analysis.molecular_weight()) + 110 * Xs)
        else:
            leader_analysis = ProteinAnalysis(self.leader, monoisotopic=True)
            columns.append(float(leader_analysis.molecular_weight()))
        #Unmodified core peptide mass
        if "X" in self.core:
            Xs = self.core.count("X")
            noXcore = self.core.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            core_analysis = ProteinAnalysis(noXcore, monoisotopic=True)
            columns.append(float(core_analysis.molecular_weight()) + 110 * Xs)
        else:
            core_analysis = ProteinAnalysis(self.core, monoisotopic=True)
            columns.append(float(core_analysis.molecular_weight()))
        #Length of Precursor
        columns.append(len(self.sequence))
        #Length of Leader
        columns.append(len(self.leader))
        #Length of Core
        columns.append(len(self.core))
        #Ratio of length of leader / length of core
        columns.append(float(len(self.core)) / float(len(self.leader)))
        #Ratio of heterocyclizable  residues / length of core
        columns.append(float(sum([self.core.count(aa) for aa in "CST"])) / len(self.core))
        #Estimated core charge at neutral pH
        #charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        #columns.append(sum([charge_dict[aa] for aa in core if aa in charge_dict]))
        #Estimated leader charge at neutral pH
        #columns.append(sum([charge_dict[aa] for aa in leader if aa in charge_dict]))
        #Estimated precursor charge at neutral pH
        #columns.append(sum([charge_dict[aa] for aa in leader+core if aa in charge_dict]))
        #Absolute value of core charge at neutral pH
        #columns.append(abs(sum([charge_dict[aa] for aa in core if aa in charge_dict])))
        #Absolute value of leader charge at neutral pH
        #columns.append(abs(sum([charge_dict[aa] for aa in leader if aa in charge_dict])))
        #Absolute value of precursor charge at neutral pH
        #columns.append(abs(sum([charge_dict[aa] for aa in leader+core if aa in charge_dict])))
        #Number in leader of each amino acid
        columns += [self.leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in leader
        columns.append(sum([self.leader.count(aa) for aa in "FWY"]))
        #Neg charged in leader
        columns.append(sum([self.leader.count(aa) for aa in "DE"]))
        #Pos charged in leader
        columns.append(sum([self.leader.count(aa) for aa in "RK"]))
        #Charged in leader
        columns.append(sum([self.leader.count(aa) for aa in "RKDE"]))
        #Aliphatic in leader
        columns.append(sum([self.leader.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in leader
        columns.append(sum([self.leader.count(aa) for aa in "ST"]))
        #Counts of AAs in core
        columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in core
        columns.append(sum([self.core.count(aa) for aa in "FWY"]))
        #Neg charged in core
        columns.append(sum([self.core.count(aa) for aa in "DE"]))
        #Pos charged in core
        columns.append(sum([self.core.count(aa) for aa in "RK"]))
        #Charged in core
        columns.append(sum([self.core.count(aa) for aa in "RKDE"]))
        #Aliphatic in core
        columns.append(sum([self.core.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in core
        columns.append(sum([self.core.count(aa) for aa in "ST"]))
        #Counts of AAs in entire precursor (leader+core)
        columns += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "FWY"]))
        #Neg charged in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "DE"]))
        #Pos charged in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "RK"]))
        #Charged in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "RKDE"]))
        #Aliphatic in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "ST"]))
        self.score = score
        self.csv_columns += [score] + columns
        return
        
        
def thioscout(core):
    """ThioScout function from Chris Schwalen to count repeat blocks"""
    #rex1 repeating Cys Ser Thr residues
    rex1 = re.compile('[C]{2,}|[S]{2,}|[T]{2,}')
    #rex2 contiguous cyclizable residues
    rex2 = re.compile('[C|S|T]{2,}')
    rexout1 = re.findall(rex1,core)
    number_of_repeat_blocks = len(rexout1)
    temp = "".join(rexout1)
    number_of_repeating_CST = str([temp.count("C"), temp.count("S"), temp.count("T")]).strip("[]")
    rexout2 = re.findall(rex2,core)
    number_of_heteroblocks = len(rexout2)
    rexout2 = np.mean([len(x) for x in rexout2])
    avg_heteroblock_length = str(rexout2).strip("[]")
    return number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks
    
