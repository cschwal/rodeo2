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
                print("ERROR:\tEsearch returns no results for query " + prot_accession_id)
                return -1
                
            IdList = record["IdList"]
            link_records = Entrez.read(Entrez.elink(dbfrom="protein",db="nuccore",id=IdList))
            nuccore_ids=[]
            if len(link_records[0]['LinkSetDb']) == 0:
                print("ERROR:\t%s has no nuccore entries..." % (prot_accession_id))
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
            print(e)
    return -3

#gb_handle should only be a handle to ONE query 
def get_record_from_gb_handle(gb_handle, nuccore_accession_id):
    """Takes an input gb_filestream and nuccore accession_id.
    Returns a record containing basic information about the accession_id.
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
                else:
                    #print("ERROR:\tCouldn't get accession ID for CDS") 
                    continue
                if 'translation' in feature.qualifiers.keys():
                    seq = feature.qualifiers['translation'][0]
                else:
                    #print("ERROR:\tCouldn't get sequence for CDS") 
                    continue
                cds = My_Record.Sub_Seq(seq_type='CDS', seq=seq, start=start, end=end, accession_id=accession_id)
                ret_record.CDSs.append(cds)
    if len(ret_record.cluster_sequence) > 0 and len(ret_record.CDSs) > 0:
        print("LOG:\tRecord made for %s" % (ret_record.cluster_accession))
        return ret_record
    else:
        print("ERROR:\tCould not process gb file for %s" % (ret_record.cluster_accession))
        return -1