#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 14:38:46 2017

@author: bryce
"""

import csv
import datetime

def draw_generic_orf_arrow_co_occur(html_file, line, query, bsc_start, bsc_end, sub_by, scale_factor):
    if line[0] != query: #new query
        html_file.write('</svg>')
        return False
    #line length should be 16. some are missing entries. i.e. only 1 pfam
    for i in range(16-len(line)):
        line.append("")
    query, genus, nuc_accession, prot_acc, start, end, direction, pfamid1, desc1, e_val1, pfamid2, desc2, e_val2, pfamid3, desc3, e_val3 = line
    start = int(start)
    end = int(end)
    #skipping fill color. all white for now
    fill_color = "white"
    html_file.write('<polygon points=\"')
    if direction == "+":
        #parth's arrowWid and arrowWid3
        arrow_wid = int((start - sub_by) * scale_factor)
        arrow_wid3 = int((end - sub_by) * scale_factor)
        if arrow_wid3 - arrow_wid < 40:
            arrow_wid2 = (arrow_wid + arrow_wid3) / 2 #middle?
        else:
            arrow_wid2 = arrow_wid3 - 20
            
        str_arrow_wid = str(arrow_wid)
        str_arrow_wid2 = str(arrow_wid2)
        str_arrow_wid3 = str(arrow_wid3)
    else:
        arrow_wid = int((end - sub_by) * scale_factor)
        arrow_wid3 = int((start - sub_by) * scale_factor)
        if arrow_wid - arrow_wid3 < 40:
            arrow_wid2 = (arrow_wid + arrow_wid3) / 2 #middle?
        else:
            arrow_wid2 = arrow_wid3 + 20
    html_file.write(str_arrow_wid + ",10 " + str_arrow_wid + ",40 "\
                        + str_arrow_wid2 + ",40 " + str_arrow_wid2 + ",50 "+ str_arrow_wid3 \
                        + ",25 " + str_arrow_wid2 + ",0 " + str_arrow_wid2 + ",10 " + str_arrow_wid + ",10")
    html_file.write('" style="fill:' + fill_color)
    html_file.write(';stroke:black;stroke-width:.5" onMouseOver="return overlib(')
    html_file.write("'" + prot_acc + " - " + pfamid1 + " : " + desc1 + "'")
    html_file.write(')" onMouseOut="return nd()"/>')
        
#    return True

def draw_orf(html_file, main_reader, main_reader_current_line,
             co_occur_reader, co_occur_current_line, start_end_dictionary):
    return_current_lines = []
    query = co_occur_current_line[0]
    #TODO direction
    html_file.write('<h2 id="' + query +'"> results for ' + query + " [" + co_occur_current_line[1] + "]</h2><p></p>")#TODO why paragraph
    html_file.write("<h3>Architecture</h3>")
    html_file.write('<svg width="1060" height="53">')
    
    bsc_start, bsc_end = start_end_dictionary[query]
    sub_by = bsc_start - 500
    scale_factor = (660./(bsc_end - bsc_start))
    draw_generic_orf_arrow_co_occur(html_file, co_occur_current_line, query, bsc_start, bsc_end, sub_by, scale_factor)
    index = 0
    for co_occur_current_line in co_occur_reader:
        index += 1
        draw_generic_orf_arrow_co_occur(html_file, co_occur_current_line, query, bsc_start, bsc_end, sub_by, scale_factor)
        if co_occur_current_line[0] == query: #If its not the current query, then stop reading it
            continue
        else:
            return_current_lines.append(co_occur_current_line)
            break #Found a new query but not at end of file
    if len(return_current_lines) == 0: #End of co_occur_file
        return_current_lines.append(False)
    
    #draw precursor html
    if main_reader_current_line[0] != query:
        return_current_lines.append(False)
    
    
    html_file.write('</svg>')
    
    return return_current_lines #End of file


#####Peptide type independent########
def write_header(html_file, record):
    html_file.write("""
    <html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css">
    </head>
    <style media="screen" type="text/css">
     
    .square {
      width: 54px;
      height: 14px;
      background-color: white;
      outline: #ffffff solid 1px;
      text-align: center;
      line-height: 14px;
      font-size: 12px;
    }
     
    table {
       font-size: 11px;
    }
     
    </style>
     
    <script src='https://img.jgi.doe.gov//js/overlib.js'></script>
    <div class="container">
    <h1 align="center" id="header">RODEO</h1>
    <div class="row">
         <div class="col-md-5">
            <h3>Parameters</h3>
                 <table class="table table-condensed" style="width:100%;">
                        <tr><th scope="row">Run Time</th><td>""" + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + """</td></tr>
                        """)
    if record.query_index != -1:
        html_file.write('<tr><th scope="row">Gene Window</th><td>+/-' + str(record.fetch_n)+ "  " + record.fetch_type + "</td></tr>\n")
    html_file.write("""
                        <tr><th scope="row">Peptide Range</th><td>""" + str(record.min_aa_seq_length) + "-" + str(record.max_aa_seq_length) + """ aa</td></tr>
                        <tr><th scope="row">Fetch Distance</th><td>""" + str(record.overlap) + """bp</td></tr>
                  </table>
         </div>
         <div class="col-md-6">
                    <div class="panel panel-default">
        </table>
        </div></div>
        </div>""")

def set_fill_color(color_string):
    if color_string ==  "C":
        return "#377eb8"
    if color_string ==  "E":
        return "#ffff33"
    if color_string ==  "B":
        return "#ff7f00"
    if color_string ==  "D":
        return "#984ea3"
    if color_string ==  "G1":
        return "gray"
    if color_string ==  "G2":
        return "gray"
    return "white"

def __main__():
    html_file = open("output/output.html", 'w') #TODO rename this
    min_aa_seq_length = 20
    max_aa_seq_length = 100
    overlap = 200
    fetch_type = "orfs"
    fetch_n = 6
    fetch_distance = "10000"
    window_size = "100000"
    peptide_types = "lasso"
    evaluate_all = False
    output_prefix = "output"
    
    write_header(html_file, fetch_n, fetch_type, min_aa_seq_length, max_aa_seq_length, overlap)
    
    
    data_file = open("output/rodeo_out_co_occur.csv")
    csv_reader = csv.reader(data_file)
    header = csv_reader.next() #Skip header
    
    start_end_dictionary = {}
    current_query = ""
    bsc_start = 0
    bsc_end = 0
    first_line = csv_reader.next()
    current_query = first_line[0]
    start_end_dictionary[current_query] = [int(first_line[4])] #start
    for line in csv_reader:
        if current_query != line[0]:
            start_end_dictionary[current_query].append([int(bsc_end)])
            current_query = line[0]
            start_end_dictionary[current_query] = int(line[4]) #start
        bsc_end = line[5] #end
    start_end_dictionary[current_query].append(int(bsc_end))
    
    main_reader = csv.reader(open(output_prefix + "/rodeo_out_main.csv")) 
    co_occur_reader = csv.reader(open(output_prefix + "/rodeo_out_co_occur.csv")) 
    main_reader.next() #Skip header
    co_occur_reader.next() #Skip header
    
    main_reader_current_line = main_reader.next()
    co_occur_current_line = co_occur_reader.next()
    
    current_lines = draw_orf(html_file=html_file, main_reader=csv_reader, 
                         main_reader_current_line=main_reader_current_line, 
                         co_occur_reader=co_occur_reader,
                         co_occur_current_line=co_occur_current_line,
                         start_end_dictionary=start_end_dictionary)
    
    while current_lines[0] and current_lines[1]:
        current_lines = draw_orf(html_file=html_file, main_reader=csv_reader, 
                         main_reader_current_line=main_reader_current_line, 
                         co_occur_reader=co_occur_reader,
                         co_occur_current_line=co_occur_current_line,
                         start_end_dictionary=start_end_dictionary)
    
    print("done")
    html_file.write("</html>")
    html_file.close()
    
__main__()