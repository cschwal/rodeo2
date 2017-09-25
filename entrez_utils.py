#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:58:10 2017

@author: bryce
"""
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from My_Record import My_Record
import logging
from rodeo_main import VERBOSITY

logger = logging.getLogger(__name__)
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

Entrez.email = 'kille2@illinois.edu'

def get_gb_handles(prot_accession_id):
    """Returns a list of .gb/.gbk filestreams from protein db accession.
    
    ERROR CODES:
        -1 = No results in protein db for Esearch on prot_accession_id
        -2 = No results in nuccore db for value obtained from protein db
        -3 = Any response failure from Entrez database (error on database side)
    """
    tries = 3 #Max number of times to try the database before it gives up
    
    for i in range(tries):
        try:
            record = Entrez.read(Entrez.esearch("protein",term=prot_accession_id))
            total_count = record["Count"]
            if int(total_count) < 1:
                logger.error("Esearch returns no results for query " + prot_accession_id)
                return -1
                
            IdList = record["IdList"]
            link_records = Entrez.read(Entrez.elink(dbfrom="protein",db="nuccore",id=IdList))
            nuccore_ids=[]
            if len(link_records[0]['LinkSetDb']) == 0:
                logger.error("%s has no nuccore entries..." % (prot_accession_id))
                return -2
            for record in link_records[0]['LinkSetDb'][0]['Link']:
                nuccore_ids.append(record['Id']) 
                
            search_handle     = Entrez.epost(db="nuccore", id=",".join(nuccore_ids))
            search_results    = Entrez.read(search_handle)
            webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
            
            batchSize = 1
            
            handles = []
            for start in range(len(nuccore_ids)):
                orig_handle = Entrez.efetch(db="nuccore", dbfrom="protein", rettype="gbwithparts", 
                                               retmode="text", retstart=start, retmax=batchSize, 
                                               webenv=webenv, query_key=query_key)
                handles.append(orig_handle)
            return handles
        except Exception as e:
            logger.error("Entrez fails with error message  \"%s\"" % (e))
    logger.error("Failed to reach Entrez database or Entrez failed to respond")
    return -3

#gb_handle should only be a handle to ONE query 
def get_record_from_gb_handle(gb_handle, nuccore_accession_id):
    """Takes an input gb_filestream and query accession_id.
    Returns a record containing basic information about the query.
    i.e. CDSs, genus/species, sequence.
    
    ERROR CODES:
        -1 = Couldn't process .gb filestream for some reason...
    """
    gb_record = SeqIO.parse(gb_handle, "genbank", generic_dna)
    for record in gb_record: #Should only be one record in gb_record
        ret_record = My_Record(nuccore_accession_id)
        ret_record.cluster_accession = record.id
        ret_record.cluster_sequence = record.seq
        ret_record.cluster_genus_species = record.annotations['organism'] #TODO verify
        for feature in record.features:
            if feature.type == 'CDS':
                start = int(feature.location.start)
                end = int(feature.location.end)
                if 'protein_id' in feature.qualifiers.keys():
                    accession_id = feature.qualifiers['protein_id'][0]
                    direction = feature.strand
                    if direction == -1:
                        tmp = start
                        start = end
                        end = tmp
                    if nuccore_accession_id == accession_id:
                        ret_record.query_index = len(ret_record.CDSs)
                else:
                    continue
                if 'translation' in feature.qualifiers.keys():
                    seq = feature.qualifiers['translation'][0]
                else:
                    continue
                cds = My_Record.Sub_Seq(seq_type='CDS', seq=seq, start=start, end=end, direction=direction, accession_id=accession_id)
                ret_record.CDSs.append(cds)
    if len(ret_record.cluster_sequence) > 0 and len(ret_record.CDSs) > 0:
        logger.info("Record made for %s" % (ret_record.cluster_accession))
        return ret_record
    else:
        logger.error("Could not process gb file for %s" % (ret_record.cluster_accession))
        return -1
