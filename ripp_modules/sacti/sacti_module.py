#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 17:41:36 2017

@author: bryce
"""
import csv
import os
import re
import numpy as np
from ripp_modules.lanthi.svm import svm_classify as svm
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.Virtual_Ripp import Virtual_Ripp
import hmmer_utils

peptide_type = "sacti"
index = 0

def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/sacti/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'Precursor Index,Classification,Precursor peptide mass,Leader peptide mass,Core peptide mass,Distance,Within 500 nt?,Within 150 nt?,Further than 1000 nt?,Ratio of N-term to 1st C 0.25<x<0.60,Ratio of N-term to 1st C <0.25 or >0.60,Three or more C,Less than 3 C,CXXXXXC,CXXXC,CXXC,CXC,CC,CCC,No cys in last 1/4th?,"2 Cys in first 2/3rds and 1 Cys in last 1/3rd",Peptide matches SboA hmm,Peptide matches SkfA hmm,Peptide matches SCIFF hmm,Cluster has PF05402 (PqqD/RRE),Cluster has PF13186 (SPASM),PF04055 (rSAM) domain start,BOOL peptidase PF05193,BOOL S8 peptidase PF00082,BOOL S41 peptidase PF03572,BOOL M16 peptidase PF00675,BOOL ABC trans PF00005,BOOL ABC mem PF00664,BOOL response reg PF00072,BOOL maj facilit PF07690,BOOL ATPase PF13304,BOOL Fer4_12 PF13353,BOOL rSAM PF04055,no  recognized peptidase,C-terminal portion is < 0.35 or > 0.65,C-terminal portion is > 0.35 and < 0.65,SS profile sum > 1,Leader length (aa),Precursor length (aa),Core length (aa),Core/precursor ratio,Core/leader ratio,Ratio of N-terminus to first Cys,number of CxnC,avg dist between CxnC,Ratio from last CxnC to C-terminus,SS profile 1 aa,SS profile 2 aa,SS profile 3 aa,SS profile 4 aa,SS profile 5 aa,SS profile 6 aa,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,peptidase PF05193,CAAX immunity PF02517,S8 peptidase PF00082,S41 peptidase PF03572,M16 peptidase PF00675,M50 peptidase PF02163,S9 peptidase PF00326,ABC trans PF00005,ABC mem PF00664,response reg PF00072,maj facilit PF07690,HATPase PF02518,ATPase PF13304,Fer4_12 PF13353,rSAM PF04055,Length of rSAM maturase'
    svm_headers = svm_headers.split(',')
    features_headers = ["Accession_id", "Genus/Species/Code", "Leader", "Core", "Start", "End", "Total Score", "Valid Precursor" ] + svm_headers
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("ripp_modules/sacti/svm/fitting_set.csv", 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(features_headers)
    svm_writer.writerow(svm_headers)
    
def ripp_write_rows(output_dir, accession_id, genus_species, list_of_rows):
    dir_prefix = output_dir + '/sacti/'
    global index
    features_csv_file = open(dir_prefix + "temp_features.csv", 'a')
    svm_csv_file = open("ripp_modules/sacti/svm/fitting_set.csv", 'a')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    for row in list_of_rows:
        features_writer.writerow([accession_id, genus_species] + row[0:5] + ["valid_precursor_placeholder", index, ''] + row[5:])
        svm_writer.writerow([index, ''] + row[5:]) #Don't include accession_id, leader, core sequence, start, end, or score
        index += 1


def run_svm(output_dir):
    svm.run_svm()
    svm_output_reader = csv.reader(open("ripp_modules/sacti/svm/fitting_results.csv"))
    final_output_writer = csv.writer(open(output_dir + "/sacti/sacti_features.csv", 'w'))
    features_reader = csv.reader(open(output_dir + "/sacti/temp_features.csv"))
    header_row = features_reader.next() #skip header
    final_output_writer.writerow(header_row)
    for row in features_reader:
        svm_output = svm_output_reader.next()[1]
        row[9] = svm_output
        if int(svm_output) == 1:
            row[6] = int(row[6]) + 10
        if int(row[6]) > 17: #CUTOFF
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
        self.peptide_type = 'sacti'
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        
    #TODO better split for sactis?    
    def set_split(self):
        self.split_index = int(len(self.sequence)*.25)
        self.leader = self.sequence[:self.split_index]
        self.core = self.sequence[self.split_index:]
        return
    
    #TODO no FIMO file?
    
    
    def set_score(self):
#        tabs = [-666] #spaceholder for index. index set later
        tabs = []
        score = 0
        if "X" in self.sequence:
            Xs = self.sequence.count("X")
            noXprecursor = self.sequence.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            precursor_analysis = ProteinAnalysis(noXprecursor, monoisotopic=True)
            tabs.append(float(precursor_analysis.molecular_weight()) + 110 * Xs)
        else:
            precursor_analysis = ProteinAnalysis(self.sequence, monoisotopic=True)
            tabs.append(float(precursor_analysis.molecular_weight()))
        #Calcd. leader peptide mass (Da)
        if "X" in self.leader:
            Xs = self.leader.count("X")
            noXleader = self.leader.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            leader_analysis = ProteinAnalysis(noXleader, monoisotopic=True)
            tabs.append(float(leader_analysis.molecular_weight()) + 110 * Xs)
        else:
            leader_analysis = ProteinAnalysis(self.leader, monoisotopic=True)
            tabs.append(float(leader_analysis.molecular_weight()))
        #Calcd. core peptide mass (Da)
        if "X" in self.core:
            Xs = self.core.count("X")
            noXcore = self.core.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            core_analysis = ProteinAnalysis(noXcore, monoisotopic=True)
            tabs.append(float(core_analysis.molecular_weight()) + 110 * Xs)
        else:
            core_analysis = ProteinAnalysis(self.core, monoisotopic=True)
            tabs.append(float(core_analysis.molecular_weight()))
            
        rSAM = "PF04055"
        if rSAM in self.pfam_2_coords.keys():
            dist = self.get_min_dist(self.pfam_2_coords[rSAM])
        else:
            dist = 999999
        tabs.append(dist)
        #rSAM within 500 nt?
        if dist < 500:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #rSAM within 150 nt?
        if dist < 150:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #rSAM further than 1000 nt?
        if dist > 1000:
            score -= 2
            tabs.append(1)
        else:
            tabs.append(0)
        
        #Ratio of N-term to 1st Cys 0.25<x<0.60; Ratio of N-term to 1st Cys <0.25 or >0.60
        if "C" not in self.sequence:
            score -= 2
            tabs += [0,1]
        elif 0.25 <= self.sequence.find("C") / float(len(self.sequence)) <= 0.60:
            score += 2
            tabs += [1,0]
        else:
            score -= 2
            tabs += [0,1]
        #Three or more Cys; Less than 3 Cys
        if self.sequence.count("C") >= 3:
            score += 4
            tabs += [1,0]
        else:
            score -= 4
            tabs += [0,1]
        #CxC/CxxC/CxxxC/CxxxxxC; #CC/CCC
        motifs = (('C[ARNDBCEQZGHILKMFPSTWYV]{5}C', 2), ('C[ARNDBCEQZGHILKMFPSTWYV]{3}C', 1), \
           ('C[ARNDBCEQZGHILKMFPSTWYV]{2}C', 1),('C[ARNDBCEQZGHILKMFPSTWYV]{1}C', 2), \
           ('CC', -2), ('CCC', -2))
        for motif in motifs:
            if re.search(motif[0], self.core) != None:
                score += motif[1]
                tabs.append(1)
            else:
                tabs.append(0)
        #No Cys in last 1/4th?
        quarter_length = (len(self.sequence) / 4) * -1
        if not "C" in self.sequence[quarter_length:]:
            score += 1
            tabs.append(1)
        else:
            score -= 1
            tabs.append(0)
        #2 Cys in first 2/3rds of precursor, 1 Cys in last 1/3rd of precursor
        two_thirds = (len(self.sequence) / 3) * 2
        if self.sequence[:two_thirds].count("C") == 2 and self.sequence[two_thirds:].count("C") == 1:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
            
        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence,n=5,e_cutoff=1, query_is_accession=False)
        pfams = []
        for pfam_dot, _, _, in precursor_hmm_info:
            pfams.append(pfam_dot.split('.')[0])
        
        #Peptide matches SboA hmm (PF11420)
        if "PF11420" in pfams:
            score += 3
            tabs.append(1)
        else:
            tabs.append(0)
        #Peptide matches SkfA hmm
        if  "TIGR04404" in pfams:  #TODO check for TIGRfams
            score += 3
            tabs.append(1)
        else:
            tabs.append(0)
        #Peptide matches SCIFF hmm
        if "TIGR03973" in pfams: #TODO check for TIGRfams
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has PqqD/RRE (PF05402)
        if "PF05402" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has SPASM domain (PF13186)
        if "PF13186" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
             
        #TODO HMMSEARCH for rSAM
        tabs.append(-666)
        #####
        
        pfams = []
        for pfam_dot in self.pfam_2_coords.keys():
            pfams.append(pfam_dot.split('.')[0])
            
        peptidase_domains = ["PF05193","PF00082","PF03572","PF00675"]
        no_peptidase = True
        for pepdom in peptidase_domains:
            if pepdom in pfams:
                score += 1
                tabs.append(1)
                no_peptidase = False
            else:
                tabs.append(0)
                
        transport_domains = ["PF00005", "PF00664"]
        for transpdom in transport_domains:
            if transpdom in pfams:
                score += 1
                tabs.append(1)
            else:
                tabs.append(0)
        #Cluster has response regulator (PF00072)
        if "PF00072" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has major facilitator (PF07690)
        if "PF07690" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has ATPase (PF13304)
        if "PF13304" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has Fer4_12 (PF13353)
        if "PF13353" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has rSAM (PF04055)
        if "PF04055" in pfams or "TIGR03975" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Cluster has no recognized peptidase
        if no_peptidase:
            score -= 2
            tabs.append(1)
        else:
            tabs.append(0)
            
        #C-terminal portion is < 0.35 or > 0.65; C-terminal portion is defined as the part from the last cysteine in the last identified Cx(n)C motif to the C-terminus
        #And criterion "C-terminal portion is > 0.35 and < 0.65"
        if "C" not in self.sequence:
            score -= 2
            tabs += [1,0]
        else:
            last_motif_C = 0
            index = -1
            for aa in reversed(self.sequence):
                if aa == "C" and "C" in self.sequence[index-6:index]:
                    last_motif_C = len(self.sequence[:index]) + 1
                index -= 1
            if not (0.35 <= last_motif_C / float(len(self.sequence)) <= 0.65):
                score -= 2
                tabs += [1,0]
            else:
                score += 3
                tabs += [1,0]
        #TODO no need for Sactiscout.py? SS profile var? ;anti-smash
        lanlower = 1
        lanupper = 6
        cysrex = '(?=(C.{%d,%d}C))' % (lanlower, lanupper)
        rex4 = re.compile(cysrex)
        totrings = rex4.findall(self.core)
        if len(totrings) > 1:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
            
        #Length of leader peptide
        tabs.append(len(self.leader))
        #Length of precursor peptide
        tabs.append(len(self.sequence))
        #Length of core peptide
        tabs.append(len(self.core))
        #Length of core / length of precursor ratio
        tabs.append(float(len(self.core)) / float(len(self.sequence)))
        #Length of core / length of leader ratio
        tabs.append(float(len(self.core)) / float(len(self.leader)))
        #Ratio of length of N-terminus to first Cys / length of precursor
        tabs.append(self.sequence.count("C") / float(len(self.sequence)))
        #Number of occurrences of CxNC motifs
        numrings, strings, avgs, cterm, profile = lanscout([self.sequence])
        if np.isnan(avgs[0]):
            avgs = [0]
        tabs.append(numrings[0])
        #Average distance between CxNC motifs
        tabs.append(avgs[0])
        #Ratio of length from last CxNC to C-terminus / length of core
        tabs.append(cterm[0])
        #Number of instances of CxNC where N = 1
        tabs.append(profile[0].split(",")[0])
        #Number of instances of CxNC where N = 2
        tabs.append(profile[0].split(",")[1])
        #Number of instances of CxNC where N = 3
        tabs.append(profile[0].split(",")[2])
        #Number of instances of CxNC where N = 4
        tabs.append(profile[0].split(",")[3])
        #Number of instances of CxNC where N = 5
        tabs.append(profile[0].split(",")[4])
        #Number of instances of CxNC where N = 6
        tabs.append(profile[0].split(",")[5])
        
        #Number in entire precursor of each amino acid
        tabs += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Number in entire precursor of each amino acid type (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
        tabs.append(sum([self.sequence.count(aa) for aa in "FWY"]))
        tabs.append(sum([self.sequence.count(aa) for aa in "DE"]))
        tabs.append(sum([self.sequence.count(aa) for aa in "RK"]))
        tabs.append(sum([self.sequence.count(aa) for aa in "RKDE"]))
        tabs.append(sum([self.sequence.count(aa) for aa in "GAVLMI"]))
        tabs.append(sum([self.sequence.count(aa) for aa in "ST"]))
        #Number in leader of each amino acid
        tabs += [self.leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Number in leader of each amino acid type (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
        tabs.append(sum([self.leader.count(aa) for aa in "FWY"]))
        tabs.append(sum([self.leader.count(aa) for aa in "DE"]))
        tabs.append(sum([self.leader.count(aa) for aa in "RK"]))
        tabs.append(sum([self.leader.count(aa) for aa in "RKDE"]))
        tabs.append(sum([self.leader.count(aa) for aa in "GAVLMI"]))
        tabs.append(sum([self.leader.count(aa) for aa in "ST"]))
        #Number in core of each amino acid
        tabs += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Number in core of each amino acid type (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
        tabs.append(sum([self.core.count(aa) for aa in "FWY"]))
        tabs.append(sum([self.core.count(aa) for aa in "DE"]))
        tabs.append(sum([self.core.count(aa) for aa in "RK"]))
        tabs.append(sum([self.core.count(aa) for aa in "RKDE"]))
        tabs.append(sum([self.core.count(aa) for aa in "GAVLMI"]))
        tabs.append(sum([self.core.count(aa) for aa in "ST"]))
        #Number of each peptidase Pfam hit (PF05193/PF00082/PF03572/PF00675/PF02517/PF02163/PF00326)
        pfams = []
        for pfam_dot in self.pfam_2_coords.keys():
            pfams.append(pfam_dot.split('.')[0])
        peptidase_domains = ["PF05193", "PF00082", "PF03572", "PF00675", "PF02517", "PF02163", "PF00326"]
        for pepdom in peptidase_domains:
            tabs.append(pfams.count(pepdom))
        #Number of each ABC transporter Pfam hit (PF00005/PF00664)    
        transp_domains = ["PF00005", "PF00664"]
        for transp_dom in transp_domains:
            tabs.append(pfams.count(transp_dom))
        #Number of each response regulator Pfam hit (PF00072)
        tabs.append(pfams.count("PF00072"))
        #Number of each major facilitator Pfam hit (PF07690)
        tabs.append(pfams.count("PF07690"))
        #Number of each ATPase Pfam hit (PF02518/PF13304)
        atpase_domains = ["PF02518", "PF13304"]
        for atpase_dom in atpase_domains:
            tabs.append(pfams.count(atpase_dom))
        #Number of each Fer4_12 Pfam hit (PF13353)
        tabs.append(pfams.count("PF13353"))
        #Number of each rSAM Pfam hit (PF04055)
        tabs.append(pfams.count("PF04055"))
        #Length of rSAM including PqqD domain
        
        #TODO part of rSAM hmmscan todo
        tabs.append(-666)
#        if len(hitends) == 0:
#            tabs.append(0)
#        else:
#            tabs.append(max(hitends))
        self.csv_columns += [score] + tabs
        self.score = score
        return

def lanscout(seq):
    lanlower = 1
    lanupper = 6
    #nis = 'ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK'
    #define lanthionine ring with a c-terminal Cys
    #lan = '(?=([T|S].{%d,%d}C))' % (lanlower, lanupper)
    cysrex = '(?=(C.{%d,%d}C))' % (lanlower, lanupper)
    
    #db = []
    strings = []
    #numringlist = []

    db = []
    avgs = []
    numrings = []
    cterm = []
    for n in range(0, len(seq)):
        core = str(seq[n][:])
        strings.append(core)

        rex1 = re.compile(r'([G|A]{1,2}C)')
        rex2 = re.compile(r'(C[G|A]{1,2})')
        rex3 = re.compile(r'(C.{2,5})(?=C)')
        rex4 = re.compile(cysrex)
        
        loc1 = []
        match1 = []
        for county in rex1.finditer(core):
            loc1.append(county.start())
            match1.append(county.group())
            
        for county in rex2.finditer(core):
            loc1.append(county.start())
            match1.append(county.group())
        tempy = []
        for m in range(len(loc1[:-1])):
            if (loc1[m+1]-loc1[m]) > 0:
                tempy.append(loc1[m+1]-loc1[m])
            else: 
                tempy.append(1)
        if len(tempy) == 0:
            avgs.append(np.nan)
        else:
            avgs.append(np.mean(tempy))
        numrings.append(len(match1))
        #print loc1
        #print match1
        cterm.append(len(rex3.split(core)[-1])/float(len(core)))
        #print cterm

        numringlist = []
        totrings = rex4.findall(core)
        size = []
        for i in range(0, len(totrings)):
            size.append(len(totrings[i]))
        db.append(size)
        rings = len(totrings)
        numringlist.append(numrings)
    
    profile = []
    for i in range(0,len(db)):
        #print db[i]
        temp = []
        for j in range(lanlower+2,lanupper+3):
            temp.append(db[i].count(j))
        profile.append(temp)
    #print profile

    for i in range(0,len(profile)):
        #profile[i] = str(profile[i])
        profile[i]=str(profile[i]).strip('[]')
    #print profile
    
    return numrings, strings, avgs, cterm, profile