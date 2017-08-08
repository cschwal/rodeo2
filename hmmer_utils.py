#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 21:38:20 2017

@author: bryce
"""

import csv
import subprocess
       
#TODO try/except blocks for error checking processes
def _generate_fasta(accession):
    out_file = open("fasta_file.tmp.fasta", 'w+')
    esearch_process = subprocess.Popen(["/home/bryce/edirect/esearch", "-db", "protein", "-query", accession],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    efetch_process = subprocess.Popen(["/home/bryce/edirect/efetch", "-format", "fasta"], 
                                      stdin=esearch_process.stdout, stdout=out_file)
    efetch_process.communicate()
    out_file.close()

def _generate_hmmer():
    hmmer_process = subprocess.call(["hmmscan", "-o", "hmm_out.tmp.tab", "--noali",
                                      "--domtblout", "pFamInfo.tmp.tab", "Pfam-A.hmm",
                                      "fasta_file.tmp.fasta"])    

#TODO try/except blocks for parsing to double check format
def _parse(n, e_cutoff):
    ret_list = []
    NUM_COLUMNS = 23
    pfam_accession = ''
    with open('pFamInfo.tmp.tab') as handle:
        for line in csv.reader(handle, delimiter=' '):
            if line[0][0] == '#':
                continue
            tmp_list = []
            for item in line:
                if item != '':
                    tmp_list.append(item)
            if tmp_list[1] != pfam_accession: #Check to make sure its not the same as the last one
                pfam_accession = tmp_list[1]
                description = ' '.join(tmp_list[NUM_COLUMNS-1:])
                e_val = tmp_list[6]
                if float(e_val) < e_cutoff and len(ret_list) < n: #Will go through one extra line but oh well #TODO
                    ret_list.append((pfam_accession, description, float(e_val)))
                else:
                    return ret_list
    return ret_list

def get_hmmer_info(accession, n=3, e_cutoff=.001): #TODO handle lists of accessions
    """Returns top n hmmscan hits with e_values lower than e_cutoff"""
    _generate_fasta(accession)
    _generate_hmmer()
    pfam_desc_list = _parse(n, e_cutoff)
    return pfam_desc_list