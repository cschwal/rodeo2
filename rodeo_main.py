#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:34:38 2017

@author: bryce
"""
import logging
import shutil
import os
import argparse
import csv
import nulltype_module
import entrez_utils
import main_html_generator
import ripp_html_generator
import My_Record
import time

VERBOSITY = logging.DEBUG


def __main__():
    logger = logging.getLogger("rodeo_main")
    logger.setLevel(VERBOSITY)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(VERBOSITY)
    
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # add formatter to ch
    ch.setFormatter(formatter)
    
    # add ch to logger
    logger.addHandler(ch)
#==============================================================================
#     First we will handle the input, whether it be an accession, a list of acc
#   or a file which itself contains a list of accessions.
#==============================================================================
    parser = argparse.ArgumentParser("Main RODEO app.")
    parser.add_argument('query', type=str,
                        help='Accession number, genbank file or .txt file with an accession or .gbk query on each line') #accession # or gi
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
    parser.add_argument('-out', '--output_dir', type=str,
                        help='Name of output folder')
    parser.add_argument('-ex', '--exhaustive', type=bool, default = False,
                        help="Score RiPPs even if they don't have a valid split site")
    args, _ = parser.parse_known_args()
    
    if args.output_dir == None:
        args.output_dir = args.query.split('.')[0] + "_rodeo_out"
    overwriting_folder = False
    try:
        os.mkdir(args.output_dir)
    except:
        overwriting_folder = True
        shutil.rmtree(args.output_dir)
        os.mkdir(args.output_dir)
    
    if overwriting_folder:
        logger.warning("Overwriting %s folder." % (args.output_dir))
        
    if not any(args.fetch_type==ft for ft in ['orfs', 'nucs']):
        logger.critical("Invalid argument for -ft/-fetch_type")
        return None
    query = args.query
    queries = []
    
    if '.txt' == query[-4:]:
        try:
            input_handle = open(query)
            for line in input_handle:
                if line.rstrip() == "\n" or line.rstrip() == "": 
                    continue
                queries.append(line.rstrip())
        except:
            logger.critical("Could not find %s" % (query)) 
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
    output_dir = args.output_dir
    exaustive = args.exhaustive
    module = nulltype_module
    
    
        
    module.main_write_headers(output_dir)
    module.co_occur_write_headers(output_dir)
    main_html = open(output_dir + "/main_results.html", 'w')
    main_html_generator.write_header(main_html, args)
    main_html_generator.write_table_of_contents(main_html, queries)
    ripp_modules = {}
    ripp_htmls = {}
    records = []
    
    HMM_SCAN_TIME = 0
    SCORE_TIME = 0
    WRITE_READ_TIME = 0
    REQUEST_TIME = 0
    SVM_TIME = 0
    for peptide_type in peptide_types:
        list_of_rows = []
        if peptide_type == "lasso":
            import ripp_modules.lasso.lasso_module as module
        elif peptide_type == "lanthi":
            import ripp_modules.lanthi.lanthi_module as module
        elif peptide_type == "sacti":
            import ripp_modules.sacti.sacti_module as module
        elif peptide_type == "thio":
            import ripp_modules.thio.thio_module as module
        else:
            logger.error("%s not in supported RiPP types" % (peptide_type))
            continue
        os.mkdir(output_dir + "/" + peptide_type)
        module.write_csv_headers(output_dir)
        ripp_modules[peptide_type] = module
        ripp_htmls[peptide_type] = open(output_dir + "/" + peptide_type + "/" + peptide_type + "_results.html", 'w')
        ripp_html_generator.write_header(ripp_htmls[peptide_type], args, peptide_type)
        ripp_html_generator.write_table_of_contents(ripp_htmls[peptide_type], queries)
    
    query_no = 1
    for query in queries:
        t0 = time.time()
        if '.gbk' != query[-4:] and '.gb' != query[-3:]: #accession_id
            logger.info("Running RODEO for %d.\t%s" % ( query_no, query))
            query_no += 1
            gb_handles = entrez_utils.get_gb_handles(query)
            if gb_handles < 0:
                error_message = "ERROR:\tentrez_utils.get_gb_handles(%s) returned with value %d"\
                      % (query, gb_handles)
                main_html_generator.write_failed_query(main_html, query, error_message)
                continue
        else:#gbk file
            try:
                gb_handles = [open(query)]
            except Exception as e:
                logger.error(e)
                continue
        REQUEST_TIME += time.time()-t0
        for handle in gb_handles:
            t0 = time.time()
            record = entrez_utils.get_record_from_gb_handle(handle, query)
            REQUEST_TIME += time.time()-t0
            if record < 0:
                error_message = "ERROR:\tentrez_utils.get_record_from_gb_handle(%s) returned with value %d"\
                  % (query, record)
                main_html_generator.write_failed_query(main_html, query, error_message)
                continue
            if fetch_type == 'orfs':
                record.trim_to_n_orfs(fetch_n, fetch_distance)
            elif fetch_type == 'nucs':
                record.trim_to_n_nucleotides(fetch_n, window_size)
            t0 = time.time()   
            record.annotate_w_hmmer()
            HMM_SCAN_TIME += time.time()-t0
            record.set_intergenic_seqs()
            record.set_intergenic_orfs(min_aa_seq_length=min_aa_seq_length, 
                                       max_aa_seq_length=max_aa_seq_length,
                                       overlap=overlap)
            module = nulltype_module
            t0 = time.time()
            for orf in record.intergenic_orfs:
                if orf.start < orf.end:
                    direction = "+"
                else:
                    direction = "-"
                row = [query, record.cluster_genus_species, record.cluster_accession, 
                       orf.start, orf.end, direction, orf.sequence]
                module.main_write_row(output_dir, row)
            for cds in record.CDSs:
                if cds.start < cds.end:
                    direction = "+"
                else:
                    direction = "-"
                row = [query, record.cluster_genus_species, record.cluster_accession,
                       cds.accession_id, cds.start, cds.end, direction]
                for pfam_acc, desc, e_val in cds.pfam_descr_list:
                    row += [pfam_acc, desc, e_val]
                module.co_occur_write_row(output_dir, row)
                WRITE_READ_TIME += time.time() - t0
            t0 = time.time()
            main_html_generator.write_record(main_html, record)
            WRITE_READ_TIME += time.time() - t0
    #==============================================================================
    #     All generic orf processing is done now. Moving on to RiPP class specifics
    #==============================================================================
            for peptide_type in peptide_types:
                list_of_rows = []
                module = ripp_modules[peptide_type]
#                module.write_csv_headers(output_dir)
                t0 = time.time()
                record.set_ripps(module, exaustive)
                record.score_ripps(module)
                SCORE_TIME += time.time() - t0
                for ripp in record.ripps[peptide_type]:
                    list_of_rows.append(ripp.csv_columns)
                t0 = time.time()
                module.ripp_write_rows(output_dir, record.cluster_accession, #cluster acc or query acc?
                                       record.cluster_genus_species, list_of_rows)
                WRITE_READ_TIME += time.time() - t0
            records.append(record)
            
            if not evaluate_all:
                #TODO users may want to see what the other entries could be?
                break
    main_html.write("</html>")   
    t0 = time.time()
    for peptide_type in peptide_types:
        module = ripp_modules[peptide_type]
        module.run_svm(output_dir)
    SVM_TIME += time.time() - t0
    t0 = time.time()
    My_Record.update_score_w_svm(output_dir, records)
    
    for peptide_type in peptide_types: 
        for record in records:
            ripp_html_generator.write_record(ripp_htmls[peptide_type], record, peptide_type)
    WRITE_READ_TIME += time.time() - t0
    logger.debug("Time for requests %f" % (REQUEST_TIME))
    logger.debug("Time for I/O %f" % (WRITE_READ_TIME))
    logger.debug("Time for HMMSCAN %f" % (HMM_SCAN_TIME))
    logger.debug("Time for SCORE/SPLIT %f" % (SCORE_TIME))
    logger.debug("Time for SVM %f" % (SVM_TIME))
    
    
if __name__ == '__main__':
    __main__()