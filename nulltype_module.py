#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:45:26 2017

@author: bryce
"""

import csv

def main_write_headers(output_filename):
    headers = ['Query', 'Genus/Species', 'Nucleotide_acc', 'start', 'end', 'dir', 'AA_seq']
    features_writer = csv.writer(open('output/' + output_filename + "_main.csv", 'w'))
    features_writer.writerow(headers)
    
def main_write_row(output_filename, row):
    features_writer = csv.writer(open('output/' + output_filename + "_main.csv", 'a'))
    features_writer.writerow(row)
    
def co_occur_write_headers(output_filename):
    headers = ['Query', 'Genus/Species', 'Nucleotide_acc', 'Protein_acc', 'start', 'end', 'dir',\
               'PfamID1', 'Description', 'E-value1',
               'PfamID2', 'Description', 'E-value2',
               'PfamID3', 'Description', 'E-value3']
    features_writer = csv.writer(open('output/' + output_filename + "_co_occur.csv", 'w'))
    features_writer.writerow(headers)
    
def co_occur_write_row(output_filename, row):
    features_writer = csv.writer(open('output/' + output_filename + "_co_occur.csv", 'a'))
    features_writer.writerow(row)