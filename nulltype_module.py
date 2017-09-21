#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:45:26 2017

@author: bryce
"""

import csv

def main_write_headers(output_folder):
    headers = ['Query', 'Genus/Species', 'Nucleotide_acc', 'start', 'end', 'dir', 'AA_seq']
    features_writer = csv.writer(open(output_folder + "/main_results.csv", 'w'))
    features_writer.writerow(headers)
    
def main_write_row(output_folder, row):
    features_writer = csv.writer(open(output_folder + "/main_results.csv", 'a'))
    features_writer.writerow(row)
    
def co_occur_write_headers(output_folder):
    headers = ['Query', 'Genus/Species', 'Nucleotide_acc', 'Protein_acc', 'start', 'end', 'dir',\
               'PfamID1', 'Description', 'E-value1',
               'PfamID2', 'Description', 'E-value2',
               'PfamID3', 'Description', 'E-value3']
    features_writer = csv.writer(open(output_folder + "/main_co_occur.csv", 'w'))
    features_writer.writerow(headers)
    
def co_occur_write_row(output_folder, row):
    features_writer = csv.writer(open(output_folder +  "/main_co_occur.csv", 'a'))
    features_writer.writerow(row)