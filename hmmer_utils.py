#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 21:38:20 2017

@author: bryce
"""

import os
import csv
import subprocess
   
pid = None    
#TODO try/except blocks for error checking processes
def _generate_fasta(accession):
    out_file = open("fasta_file.tmp.fasta", 'w+')
    esearch_process = subprocess.Popen(["/home/bryce/edirect/esearch", "-db", "protein", "-query", accession],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    efetch_process = subprocess.Popen(["/home/bryce/edirect/efetch", "-format", "fasta"], 
                                      stdin=esearch_process.stdout, stdout=out_file)
    efetch_process.communicate()
    out_file.close()

def _generate_fasta_from_string(query_string):
    handle = open(pid+"fasta_file.tmp.fasta", "w+")
    handle.write(">temp_string\n" + query_string)

def _generate_hmmer():
    hmmer_process = subprocess.call(["hmmscan", "-o", pid+"hmm_out.tmp.tab", "--noali",
                                      "--domtblout", pid+"pFamInfo.tmp.tab", "Pfam-A.hmm",
                                      pid+"fasta_file.tmp.fasta"])    

#TODO try/except blocks for parsing to double check format
def _parse(n, e_cutoff):
    ret_list = []
    NUM_COLUMNS = 23
    pfam_accession = ''
    with open(pid+'pFamInfo.tmp.tab') as handle:
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

def get_hmmer_info(query, n=3, e_cutoff=.001, query_is_accession=False): #TODO handle lists of accessions
    """Returns top n hmmscan hits with e_values lower than e_cutoff"""
    global pid
    pid = str(os.getpid())
    if query_is_accession:
        _generate_fasta(query)
    else:
        _generate_fasta_from_string(query)
    _generate_hmmer()
    pfam_desc_list = _parse(n, e_cutoff)
    os.remove(pid+'pFamInfo.tmp.tab')
    os.remove(pid+"fasta_file.tmp.fasta")
    os.remove(pid+"hmm_out.tmp.tab")
    return pfam_desc_list