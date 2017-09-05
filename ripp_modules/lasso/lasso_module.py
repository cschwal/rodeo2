#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:40:56 2017

@author: bryce
"""
import csv
import os
import re
import tempfile
import subprocess
from ripp_modules.lasso.svm import svm_classify as svm
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.Virtual_Ripp import Virtual_Ripp
index = 0

def write_csv_headers(output_filename):
    dir_prefix = 'output/lasso/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    headers = ['Accession_id', 'Genus/Species', 'Leader', 'Core', 'Start', 'End' , 'Primary Key', 'VALID PRECURSOR', 'Calcd. Lasso Mass (Da)', 'Distance', 'Within 500 nt?', 'Within 150 nt?', 'Further than 1000 nt?', 'Core has 2 or 4 Cys?', 
     'Leader longer than core?', 'Plausible lasso ring?', 'Leader has GxxxxxT motif?', 'Core starts with G?', 'Core and BGC in same direction?',
	 'Ratio leader/core < 2 and > 0.5'	, 'Core starts with Cys and even number of Cys?', 'No Gly in core?', 'Core has at least 1 aromatic aa?',
     'Core has at least 2 aromatic aa?', 'Core has odd number of Cys?', 'Leader has Trp?', 'Leader has Lys?', 'Leader has Cys?',
     'Cluster has PF00733?', 'Cluster has PF05402?', 'Cluster has PF13471?', 'Leader has LxxxxxT motif?', 'Core has adjacent identical aas (doubles)?',
     'Core length (aa)', 'Leader length (aa)', 'Precursor length (aa)', 'Leader/core ratio', 'Number of Pro in first 9 aa of core?',
     'Estimated core charge', 'Estimated leader charge', 'Estimated precursor charge', 'Absolute value of core charge',
     'Absolute value of leader charge', 'Absolute value of precursor charge',
     'LEADER A', 'R',	'D',	'N',	'C',	'Q',	'E',	'G',	'H',	'I',	'L',	'K',	'M',	'F',	'P',	'S',	'T',	'W',	'Y',	'V',	
     'Aromatics', 'Neg charged', 'Pos charged', 'Charged', 'Aliphatic', 'Hydroxyl',
     'CORE A', 'R',	'D',	'N',	'C',	'Q',	'E',	'G',	'H',	'I',	'L',	'K',	'M',	'F',	'P',	'S',	'T',	'W',	'Y',	'V',	
     'Aromatics', 'Neg charged', 'Pos charged', 'Charged', 'Aliphatic', 'Hydroxyl',
     'FIRST CORE RESIDUE A', 'R',	'D',	'N',	'C',	'Q',	'E',	'G',	'H',	'I',	'L',	'K',	'M',	'F',	'P',	'S',	'T',	'W',	'Y',	'V',	
     'PRECURSOR A', 'R',	'D',	'N',	'C',	'Q',	'E',	'G',	'H',	'I',	'L',	'K',	'M',	'F',	'P',	'S',	'T',	'W',	'Y',	'V',	
     'Aromatics', 'Neg charged', 'Pos charged', 'Charged', 'Aliphatic', 'Hydroxyl',
     'Motif1?', 'Motif2?', 'Motif3?', 'Motif4?', 'Motif5?', 'Motif6?', 'Motif7?', 'Motif8?', 'Motif9?',
     'Motif10?', 'Motif11?', 'Motif12?', 'Motif13?', 'Motif14?', 'Motif15?', 'Motif16?',	
     'Total motifs hit',	'Score, Motif1',	'Score, Motif2',	'Score, Motif3',	'Score, Motif4',	
     'Score, Motif5',	'Score, Motif6',	'Score, Motif7',	'Score, Motif8',	'Score, Motif9',	
     'Score, Motif10', 'Score, Motif11', 'Score, Motif12', 'Score, Motif13', 'Score, Motif14',	
     'Score, Motif15', 'Score, Motif16', 'Sum of MEME scores',
     'No Motifs?', 'Alternate Start Codon?', 'Total Score', 'SVM Classification']
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("ripp_modules/lasso/svm/fitting_set.csv", 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(headers)
    svm_writer.writerow(headers[6:-2])#Don't include accession_id, genus/species,
                                        #leader, core sequence, score, or svm classification


def ripp_write_rows(output_filename, accession_id, genus_species, list_of_rows):
    dir_prefix = 'output/lasso/'
    global index
    features_csv_file = open(dir_prefix + "temp_features.csv", 'a')
    svm_csv_file = open("ripp_modules/lasso/svm/fitting_set.csv", 'a')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    for row in list_of_rows:
        features_writer.writerow([accession_id, genus_species] + row[0:4] + [index, ''] + row[4:])
        svm_writer.writerow([index, ''] + row[4:-1]) #Don't include accession_id, leader, core sequence, start, end, or score
        index += 1
        
def run_svm():
    svm.run_svm()
    svm_output_reader = csv.reader(open("ripp_modules/lasso/svm/fitting_results.csv"))
    final_output_writer = csv.writer(open("output/lasso/lasso_features.csv", 'w'))
    features_reader = csv.reader(open("output/lasso/temp_features.csv"))
    header_row = features_reader.next() #skip header
    final_output_writer.writerow(header_row)
    for row in features_reader:
        svm_output = svm_output_reader.next()[1]
        row.append(svm_output)
        if int(svm_output) == 1:
            row[-2] = int(row[-2]) + 10
        if int(row[-2]) > 17: #CUTOFF
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
        self.peptide_type = 'lasso'
        self.set_split()
        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        
    #Used for scoring to determine minimum distance from this ripp to a set
    #of coordinates.
    #For example, if coords_list was a list of the coordinates of some pfam,
    #this function would return the distance to the closest coordinate of that pfam
#    def get_min_dist(self, coords_list):
#        if coords_list == []:
#            return None
#        min_dist = abs(self.start-coords_list[0][0])
#        for coord in coords_list:
#            min_dist = min(abs(self.start-coord[0]), abs(self.end-coord[0]),
#                           abs(self.start-coord[1]), abs(self.end-coord[1]),
#                           min_dist)
#        return min_dist
#    
#    #TODO error catching
#    
#    
#        
#    def run_fimo_simple(self):
#        "Run FIMO"
#        #TODO change to temp file
#        with open("TEMP.seq", 'w+') as tfile:
#            tfile.write(">query\n%s" % (self.sequence))
#    
##       command = ["$HOME/meme/bin/fimo", "--text", "--verbosity", "1", self.query_motif_file, tfile.name]
#        command = ["$HOME/meme/bin/fimo --text --verbosity 1 " + self.query_motif_file + " " + "TEMP.seq"]
#        try:
#            out, err, retcode = execute(command)
#        except OSError:
#            print("ERROR:\tCould not run FIMO")
#            return ""
#        if retcode != 0:
#            print('FIMO returned %d: %r while searching %r', retcode,
#                            err, self.query_motif_file)
#            return []
#        return out   

    
####################################################################################
################################ Lasso functions ###################################
####################################################################################

        
#    def set_leader_core(self):
#        self.split()
#        if self.split_index == -1:
#            self.core = self.sequence
#            self.leader = ''
#        else:
#            self.leader = self.sequence[0:self.split_index]
#            self.core = self.sequence[self.split_index:]
        
        
    def set_split(self):
        #TODO add more regexes
        match = re.search('(T[A-Z]{7,10}(D|E)[A-Z]{5,20}\*)', self.sequence + '*')
        if match is None:
            self.split_index = -1
        else:
            self.split_index = match.start() + 2
        
        if self.split_index == -1:
            self.core = self.sequence
            self.leader = ''
        else:
            self.leader = self.sequence[0:self.split_index]
            self.core = self.sequence[self.split_index:]
                
    def get_fimo_score(self):
        #TODO better way for this?
        fimo_output = self.run_fimo_simple()
        fimo_motifs = [int(line.partition("\t")[0]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()]
        fimo_scores = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        #Calculate score
        motif_score = 0
        if 2 in fimo_motifs:
            motif_score += 4
        elif len(fimo_motifs) > 0:
            motif_score += 2
        else:
            motif_score += -1
        return fimo_motifs, motif_score, fimo_scores
    
    def set_monoisotopic_mass(self):
        self._set_number_bridges()
        CC_mass = 2*self._num_bridges
        # dehydration indicative of cyclization     
        bond = 18.02
        monoisotopic_mass = ProteinAnalysis(self.core.replace('X', ''), monoisotopic=True).molecular_weight()
        self._monoisotopic_weight = monoisotopic_mass + CC_mass - bond
        
    def _set_number_bridges(self):
        '''
        Predict the lassopeptide number of disulfide bridges
        '''
        self._num_bridges = 0
        if self.core.count("C") == 2:
            self._num_bridges = 1
        if self.core.count("C") >= 4:
            self._num_bridges = 2
        return self._num_bridges
    
    def set_score(self):
        scoring_csv_columns = []
        cPFam = "PF00733"
        ePFam = "PF05402"
        bPFam = "PF13471"
        self.score = 0
        scoring_csv_columns.append(self._monoisotopic_weight)
        
        bsp_coords = []
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in [cPFam, ePFam, bPFam]):
                bsp_coords += self.pfam_2_coords[pfam]
        min_distance = self.get_min_dist(bsp_coords)
        if min_distance is None:
            scoring_csv_columns.append(999999)
        else:
            scoring_csv_columns.append(min_distance)
        
        has_bPFam = False
        has_cPFam = False
        has_ePFam = False
        within_500 = False
        within_150 = False
        within_1000 = False
        cyclase_same_strand= False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in [cPFam, ePFam, bPFam]):
                #TODO 'in' or 'is'/==
                if bPFam in pfam:
                    has_bPFam = True
                elif ePFam in pfam:
                    has_ePFam = True
                elif cPFam in pfam:
                    if self.start < self.end and self.pfam_2_coords[pfam][0] < self.pfam_2_coords[pfam][0] or\
                       self.start > self.end and self.pfam_2_coords[pfam][0] > self.pfam_2_coords[pfam][0]:
                       cyclase_same_strand = True
                        
                    has_cPFam = True
                dist = self.get_min_dist(self.pfam_2_coords[pfam])
                if dist < 1000:
                    within_1000 = True
                    if dist < 500:
                        within_500 = True
                        if dist < 150:
                            within_150 = True
                            self.score += 2
                        else:
                            self.score += 1
                        break
        if within_500:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_150:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if not within_1000:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if self.core.count("C") == 2 or self.core.count("C") == 4:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if len(self.leader) > len(self.core):
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if any(aa in self.core[6:9] for aa in ["D", "E"]): 
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        #Check for GxxxxxT motif
        match = re.search('G.{5}T', self.leader)
        if match is not None:
            self.score += 3
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if self.core[0] == "G":
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        #check if peptide and lasso cyclase are on same strand +1
        if cyclase_same_strand:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if 0.5 < len(self.leader)/float(len(self.core)) < 2:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if self.core[0] == 'C' and self.core.count('C') % 2 == 0:
            self.score += 0
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if "G" not in self.core:
            self.score += -4
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Check for aromatic residues
        if any(aa in self.core for aa in ["H", "F", "Y", "W"]):
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if sum(aa in self.core for aa in ["H", "F", "Y", "W"]) >= 2:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if self.core.count("C") % 2 == 1:
            self.score += -2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
        if "W" in self.leader:
            self.score += -1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if "K" in self.leader:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if "C" in self.leader:
            self.score += -2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if has_cPFam:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if has_ePFam:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if has_bPFam:
            scoring_csv_columns.append(1)
        else:
            self.score += -2
            scoring_csv_columns.append(0)

        match = re.search('L.{5}T', self.leader)
        if match is not None:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        #'Core has adjacent identical aas'
        prev_aa = ''
        adjacent_aas= False
        for aa in self.core:
            if aa == prev_aa:
                adjacent_aas = True
                break
            else:
                prev_aa = aa
        if adjacent_aas:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        scoring_csv_columns.append(len(self.core))
        scoring_csv_columns.append(len(self.leader))
        scoring_csv_columns.append(len(self.sequence))
        scoring_csv_columns.append(float(len(self.leader))/len(self.core))
        scoring_csv_columns.append(self.core[:9].count('P'))
        
        charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        scoring_csv_columns.append(sum([charge_dict[aa] for aa in self.core if aa in charge_dict]))
        #Estimated leader charge
        scoring_csv_columns.append(sum([charge_dict[aa] for aa in self.leader if aa in charge_dict]))
        #Estimated precursor charge
        scoring_csv_columns.append(sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict]))
        #Absolute value of core charge
        scoring_csv_columns.append(abs(sum([charge_dict[aa] for aa in self.core if aa in charge_dict])))
        #Absolute value of leader charge
        scoring_csv_columns.append(abs(sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])))
        #Absolute value of precursor charge
        scoring_csv_columns.append(abs(sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict])))
        #Counts of AAs in leader
        scoring_csv_columns += [self.leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "FWY"]))
        #Neg charged in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "DE"]))
        #Pos charged in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "RK"]))
        #Charged in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "RKDE"]))
        #Aliphatic in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "ST"]))
        #Counts of AAs in core
        scoring_csv_columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "FWY"]))
        #Neg charged in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "DE"]))
        #Pos charged in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "RK"]))
        #Charged in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "RKDE"]))
        #Aliphatic in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "ST"]))
        #Counts (0 or 1) of amino acids within first AA position of core sequence
        scoring_csv_columns += [self.core[0].count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Counts of AAs in leader+core
        scoring_csv_columns += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"] #Temp to work with current training CSV
        #Aromatics in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "FWY"]))
        #Neg charged in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "DE"]))
        #Pos charged in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "RK"]))
        #Charged in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "RKDE"]))
        #Aliphatic in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "ST"]))
        
        
        #TODO MEME stuff
        fimo_motifs, motif_score, fimo_scores = self.get_fimo_score()
        self.fimo_motifs = fimo_motifs
        self.fimo_scores = fimo_scores
        self.score += motif_score
        #Motifs
        scoring_csv_columns += [1 if motif in fimo_motifs else 0 for motif in range(1, 17)]
        #Total motifs hit
        scoring_csv_columns.append(len(fimo_motifs))
        #Motif scores
        scoring_csv_columns += [fimo_scores[motif] if motif in fimo_motifs else 0 for motif in range(1, 17)]
        #Sum of MEME scores
        scoring_csv_columns.append(sum([fimo_scores[motif] if motif in fimo_motifs else 0 for motif in range(1, 17)]))
        #No Motifs?
        if len(fimo_motifs) == 0:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if self.leader[0] != 'M':
            scoring_csv_columns.append(1)
            self.score += -1
        else:
            scoring_csv_columns.append(0)
        scoring_csv_columns.append(self.score)
        self.csv_columns += scoring_csv_columns

def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, 
                                shell=True)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        print(e)
        raise