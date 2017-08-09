#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:34:38 2017

@author: bryce
"""

import argparse
import csv
import nulltype_module
import entrez_utils

def __main__():
#==============================================================================
#     First we will handle the input, whether it be an accession, a list of acc
#   or a file which itself contains a list of accessions.
#==============================================================================
    parser = argparse.ArgumentParser("Main RODEO app.")
    parser.add_argument('acc', type=str,
                        help='Accession number or or .txt file of accession numbers') #accession # or gi
    parser.add_argument('-u', '--upper_limit', type=int, default=100, 
                        help='Maximum size of potential ORF') #better word for potential?
    parser.add_argument('-l', '--lower_limit', type=int, default=20, 
                        help='Minimum size of potential ORF') 
    parser.add_argument('-o', '--overlap', type=int, default=50, 
                        help='Maximum overlap of search with existing CDSs')
    parser.add_argument('-ft', '--fetch_type', type=str, default='orfs', 
                        help='Type of window specification.\n' +
                        '\'orfs\' will make the window +/- n orfs from the query.\n' +
                        '\'nucs\' will make the window +/- n nucleotides from the query')
    parser.add_argument('-fn', '--fetch_n', type=int, default=6, 
                        help='The \'n\' variable for the -ft=orfs')
    parser.add_argument('--window_size', type=int, default=5000, 
                        help='The \'n\' variable for the -ft=nucs')
    parser.add_argument('-fd', '--fetch_distance', type=int, default=300, 
                        help='Number of nucleotides to fetch outside of window')
    parser.add_argument('-pt', '--peptide_types', nargs='*', default = [],
                        help='Type(s) of peptides to score.')
    parser.add_argument('-ea', '--evaluate_all', type=bool, default = False,
                        help='Evaluate all duplicates if accession id corresponds to duplicate entries')
    parser.add_argument('-out', '--output_filename', type=str, default = "rodeo_out",
                        help='Output filename')
    args, _ = parser.parse_known_args()
    if not any(args.fetch_type==ft for ft in ['orfs', 'nucs']):
        print("ERROR:\t Invalid argument for -ft/-fetch_type")
        return None
    query = args.acc
    queries = []
    #WARNING:
        #POSSIBLE INFINITE LOOP if users have circular file references
    if '.txt' == query[-4:]:
        try:
            input_handle = open(query)
            for line in input_handle:
                queries.append(line.rstrip())
        except:
            print("ERROR:\tCould not find %s" % (query)) 
    else:
        queries.append(query)
    min_aa_seq_length = args.lower_limit
    max_aa_seq_length = args.upper_limit
    overlap = args.overlap
    fetch_type = args.fetch_type
    fetch_n = args.fetch_n
    fetch_distance = args.fetch_distance
    window_size = args.window_size
    peptide_types = args.peptide_types
    evaluate_all = args.evaluate_all
    output_filename = args.output_filename
    module = nulltype_module
    
    module.main_write_headers(output_filename)
    module.co_occur_write_headers(output_filename)
    for peptide_type in peptide_types:
            if peptide_type == "lasso":
                from ripp_modules.lasso import lasso_module
                lasso_module.write_csv_headers(output_filename)
                
    for query in queries:
        if not any(query[-4:] == extension for extension in ['.gb', '.gbk']):#accession_id
            print("LOG:\tRunning RODEO for %s" % ( query))
            gb_handles = entrez_utils.get_gb_handles(query)
            if gb_handles < 0:
                print("ERROR:\tentrez_utils.get_gb_handles(%s) returned with value %d"\
                      % (query, gb_handles))
                continue
        else:#gbk file
            try:
                gb_handles = [open(query)]
            except Exception as e:
                print(e)
                continue
        for handle in gb_handles:
            record = entrez_utils.get_record_from_gb_handle(handle, query)
            if record < 0:
                print("ERROR:\tentrez_utils.get_record_from_gb_handle(%s) returned with value %d"\
                  % (query, record))
                continue
            if fetch_type == 'orfs':
                record.trim_to_n_orfs(fetch_n, fetch_distance)
            elif fetch_type == 'nucs':
                record.trim_to_n_nucleotides(fetch_n, window_size)
            record.annotate_w_hmmer()
            record.set_intergenic_seqs()
            record.set_intergenic_orfs(min_aa_seq_length=min_aa_seq_length, 
                                       max_aa_seq_length=max_aa_seq_length,
                                       overlap=overlap)
            module = nulltype_module
            for orf in record.intergenic_orfs:
                if orf.start < orf.end:
                    direction = "+"
                else:
                    direction = "-"
                row = [query, record.cluster_genus_species, record.cluster_accession, 
                       orf.start, orf.end, direction, orf.sequence]
                module.main_write_row(output_filename, row)
            for cds in record.CDSs:
                if cds.start < cds.end:
                    direction = "+"
                else:
                    direction = "-"
                row = [query, record.cluster_genus_species, record.cluster_accession,
                       cds.accession_id, cds.start, cds.end, direction]
                for pfam_acc, desc, e_val in cds.pfam_descr_list:
                    row += [pfam_acc, desc, e_val]
                module.co_occur_write_row(output_filename, row)
            
    #==============================================================================
    #     All generic orf processing is done now. Moving on to RiPP class specifics
    #==============================================================================
            for peptide_type in peptide_types:
                list_of_rows = []
                if peptide_type == "lasso":
                    module = lasso_module
                record.set_ripps(module)
                record.score_ripps(module)
                for ripp in record.ripps:
                    list_of_rows.append(ripp.csv_columns)
                module.ripp_write_rows(output_filename, record.cluster_accession, #cluster acc or query acc?
                                       record.cluster_genus_species, list_of_rows)
            if not evaluate_all:
                #TODO users may want to see what the other entries could be?
                break
    for peptide_type in peptide_types:
        if peptide_type == "lasso":
                    module = lasso_module
        module.run_svm()
__main__()