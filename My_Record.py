#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 21:09:19 2017

@author: bryce
"""
import hmmer_utils

class My_Record(object):
    """ """
    #TODO get genus and species frecom gb file
    def __init__(self, query_accession_id):
        self.query_accession_id = query_accession_id
        self.cluster_accession = ""
        self.cluster_sequence = ""
        self.cluster_length = ""
        self.genus = ""
        self.species = ""
        self.CDSs = []
        self.intergenic_seqs = []
        self.intergenic_orfs = []
        self.window_start = 0
        self.window_end = 0
        self.start_codons = ['ATG','GTG', 'TTG']
        self.stop_codons = ['TAA','TAG','TGA']
        self.ripps = []
        self.ripp_type_to_ripps = {}
        
    class Sub_Seq(object):
        """Useful for storing subsequences and their coordinates"""
        def __init__(self, seq_type, seq, start, end, accession_id=None):
            self.start = start
            self.end = end
            self.sequence = seq
            self.accession_id = accession_id
            self.type = seq_type ##aa, nt etc.
    
    def _get_query_index(self):
        """Get the index of the query CDS in the list of CDSs"""
        query_index = 0
        for cds in self.CDSs:
            if cds.accession_id == self.query_accession_id:
                return query_index
            else:
                query_index += 1
                
    def _clean_CDSs(self):
        """Called after trimming to get rid of CDSs not in the window"""
        i = 0
        while i < len(self.CDSs):
            #TODO entire CDS outside of window
            if self.CDSs[i].end < self.window_start or self.CDSs[i].start < self.window_start \
                or self.CDSs[i].start > self.window_end or self.CDSs[i].end > self.window_end:
                del self.CDSs[i]
            else:
                i += 1
    
    #TODO cutoff or keep if in middle of gene?
    def trim_to_n_nucleotides(self, n):
        """Trim the window down to -n nucleotides of the start of the 
        query CDS and +n nucleotides of the end of the CDS"""
        query_index = self._get_query_index()
        self.window_start = max(0, self.CDSs[query_index].start - n)
        self.window_end = min(len(self.cluster_sequence), 
                              self.CDSs[query_index].end + n)
        self._clean_CDSs()
        
    def trim_to_n_orfs(self, n, fetch_distance):
        """Trims the window to +- n CDSs from the query CDS, adding on "fetch_distance"
        nucleotides to each end of the window past the end CDSs"""
        query_index = self._get_query_index()
        first_cds = max(0, query_index - n)
        last_cds = min(len(self.CDSs)-1, query_index + n)
        self.window_start = min(self.CDSs[first_cds].start, self.CDSs[first_cds].end)
        self.window_end = max(self.CDSs[first_cds].start, 
                              self.CDSs[last_cds].end)
        self.window_start = max(0, self.window_start - fetch_distance)
        self.window_end = min(len(self.cluster_sequence), 
                              self.window_end + fetch_distance)
        self._clean_CDSs()
        
    def annotate_w_hmmer(self):
        self.pfam_2_coords = {}
        for CDS in self.CDSs:
            CDS.pfam_descr_list = hmmer_utils.get_hmmer_info(CDS.accession_id) #Possible input for n and e_cutoff here
            for annot in CDS.pfam_descr_list:
                if annot[0] not in self.pfam_2_coords.keys():
                    self.pfam_2_coords[annot[0]] = []
                self.pfam_2_coords[annot[0]].append((CDS.start, CDS.end))

    def set_intergenic_seqs(self):
        """Sets the sequences between called CDSs"""
        
        #First need to check if we have trimmed our sequence yet
        if self.window_end == 0:
            self.window_end = len(self.cluster_sequence)
        start = self.window_start
        for cds in self.CDSs:
            if len(cds.pfam_descr_list) == 0:
                continue
            end = cds.start
            if end > start:
                #end == start could happen if the first cds starts at 0
                nt_seq = self.cluster_sequence[start:end]
                intergenic_sequence = self.Sub_Seq(seq_type='IGS',
                                                   seq=nt_seq, 
                                                   start=start,
                                                   end=end)
                self.intergenic_seqs.append(intergenic_sequence)
            start = cds.end
        nt_seq = self.cluster_sequence[start:]
        end = self.window_end
        if end > start:
            intergenic_sequence = self.Sub_Seq(seq_type='IGS',
                                                   seq=nt_seq, 
                                                   start=start,
                                                   end=end)
            self.intergenic_seqs.append(intergenic_sequence)
        return
        
    def set_intergenic_orfs(self, min_aa_seq_length, max_aa_seq_length, overlap):
        """Examines intergenic sequences to determine whether or not there are ORFs
        that code for valid aa sequences"""
        for strand, sequence in [(1, self.cluster_sequence),
                                         (-1, self.cluster_sequence.reverse_complement())]:
            for intergenic_seq in self.intergenic_seqs:
                #Do start codons iteratively to ensure you don't skip any
                #In theory, should only be a 3x slowdown in runtime MAX.
                for start_codon in self.start_codons:
                    if strand == -1:
                        next_start = max(len(sequence) - intergenic_seq.end - overlap, 
                                         len(sequence)-self.window_end)
                        intergenic_seq_end = min(len(sequence) - intergenic_seq.start + overlap,
                                                 len(sequence) - self.window_start)
                    else:
                        next_start = max(intergenic_seq.start - overlap, self.window_start)
                        intergenic_seq_end = min(self.window_end, intergenic_seq.end + overlap)
                        
                    start = 0
                    #Stay in the loop until we can't find a stop codon
                    while start < intergenic_seq_end:
                        start = sequence.find(start_codon, next_start)
                        if start == -1:
                            #Can't find any more of this codon
                            break
                        next_start = start + 1
                        end = start 
                        found_stop = False
                        while end < len(sequence) and end < intergenic_seq_end:
                            codon = sequence[end:end+3]
                            if str(codon) in self.stop_codons:
                                found_stop = True
                                break
                            else:
                                end = end + 3
                        #If we didn't find a stop codon, do what?
                        if not found_stop:
                            continue
                        if not (min_aa_seq_length < (end-start)/3 < max_aa_seq_length):
                            continue
                        nt_subsequence = sequence[start:end] #Add 3 if you want to inclue stop codon
                        aa_sequence = nt_subsequence.translate(11)
                        #Get nucleotide coords for original strand
                        if strand == -1:
                            old_end = end
                            end = len(sequence) - start
                            old_start = len(sequence) - old_end - 3
                            potential_orf = self.Sub_Seq('ORF',
                                                         aa_sequence,
                                                         end,
                                                         old_start)
                        else:
                            end = end + 3
                            potential_orf = self.Sub_Seq('ORF',
                                                         aa_sequence,
                                                         start,
                                                         end)
                        if potential_orf.start > potential_orf.end:
                            upstream_sequence = str(self.cluster_sequence[potential_orf.start+4:potential_orf.start+13].reverse_complement())
                        else:
                            upstream_sequence = str(self.cluster_sequence[potential_orf.start-13:potential_orf.start+4])
                        potential_orf.upstream_sequence = upstream_sequence
                        self.intergenic_orfs.append(potential_orf)
        self.intergenic_orfs.sort(key=lambda seq: seq.start)
        #Get rid of duplicates. Duplicate ORFs will appear when the overlap is
        #set such that two intergenic sequences are expanded to a point where 
        #they share nucleotides with eachother
        i = 1
        while i < len(self.intergenic_orfs):
            if self.intergenic_orfs[i].start == self.intergenic_orfs[i-1].start:
                del self.intergenic_orfs[i]
            i += 1
        return
    
    def set_ripps(self, module):
        self.ripps = []
        for orf in self.intergenic_orfs:
            ripp = module.Ripp(orf.start, orf.end, str(orf.sequence), orf.upstream_sequence, self.pfam_2_coords)
            if ripp.split_index != -1:
                self.ripps.append(ripp)
                
    def score_ripps(self, module):
        for ripp in self.ripps:
            ripp.score_leader_core()
                
                    
    def print_info(self):
        print("="*50)
        counter = 0
        print("CDSs")
        for sub_seq in self.CDSs:
            print(counter)
            counter+=1
            if sub_seq.accession_id == self.query_accession_id:
                print("QUERY Accession:  " + sub_seq.accession_id)
            else:
                print("Accession:  " + sub_seq.accession_id)
            print("Coords:  " + str(sub_seq.start) + " to " +  str(sub_seq.end))
        print("="*50)
        print("IGSs")
        for sub_seq in self.intergenic_seqs:
            print(counter)
            counter+=1
            print("Coords:  " + str(sub_seq.start) + " to " +  str(sub_seq.end))
        print("="*50)
        print("Intergenic ORFs")
        for sub_seq in self.intergenic_orfs:
            print(counter)
            counter+=1
            if sub_seq.end < sub_seq.start:
                strand = -1
            else:
                strand = 1
            print("Potential ORF of length " + str(len(sub_seq.sequence)-1) + 
                  " found at " + str(sub_seq.start) + ":" + str(sub_seq.end) +
                  " on strand " + str(strand))
            print(sub_seq.sequence + '\n')
        print("="*50)