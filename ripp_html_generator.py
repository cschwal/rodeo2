#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 23:19:54 2017

@author: bryce
"""

import csv
import datetime
import My_Record
from decimal import Decimal

def write_header(html_file, args, peptide_type):
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
    html_file.write('<tr><th scope="row">Gene Window</th><td>+/-' + str(args.fetch_n)+ "  " + args.fetch_type + "</td></tr>\n")
    html_file.write("""
                        <tr><th scope="row">Peptide Range</th><td>""" + str(args.lower_limit) + "-" + str(args.upper_limit) + """ aa</td></tr>
                        <tr><th scope="row">Fetch Distance</th><td>""" + str(args.overlap) + """bp</td></tr>
                        <tr><th scope="row">Peptide Type</th><td>""" + peptide_type + """</td></tr>
                  </table>
         </div>
         <div class="col-md-6">
                    <div class="panel panel-default">
        </table>
        </div></div>
        </div>""")

def draw_CDS_arrow(main_html, cds, sub_by, scale_factor):
    fill_color = "white" #TODO make changes
    if cds.direction == "+":
        start = cds.start
        end = cds.end
    else:
        end = cds.start
        start = cds.end
    #HMM info
    if len(cds.pfam_descr_list) == 0:
        pfamID = "No Pfam match"
        pfam_desc = ""
    else:
        pfamID = cds.pfam_descr_list[0][0]
        pfam_desc = cds.pfam_descr_list[0][1]
    main_html.write('<polygon points=\"')
    arrow_wid = int((start - sub_by) * scale_factor)
    arrow_wid3 = int((end - sub_by) * scale_factor)
    if arrow_wid3 - arrow_wid < 40:
        arrow_wid2 = (arrow_wid + arrow_wid3) / 2 #middle?
    else:
        arrow_wid2 = arrow_wid3 - 20
        
    str_arrow_wid = str(arrow_wid)
    str_arrow_wid2 = str(arrow_wid2)
    str_arrow_wid3 = str(arrow_wid3)
    main_html.write(str_arrow_wid + ",10 " + str_arrow_wid + ",40 "\
                        + str_arrow_wid2 + ",40 " + str_arrow_wid2 + ",50 "+ str_arrow_wid3 \
                        + ",25 " + str_arrow_wid2 + ",0 " + str_arrow_wid2 + ",10 " + str_arrow_wid + ",10")
    main_html.write('" style="fill:' + fill_color)
    main_html.write(';stroke:black;stroke-width:.5" onMouseOver="return overlib(')
    main_html.write("'" + cds.accession_id + " - " + pfamID + " : " + pfam_desc + "'")
    main_html.write(')" onMouseOut="return nd()"/>')

def draw_orf_arrow(main_html, orf, sub_by, scale_factor, index):
    fill_color = "white"
    start = orf.start
    end = orf.end
    arrow_wid = int((start - sub_by) * scale_factor)
    arrow_wid3 = int((end - sub_by) * scale_factor)
    if arrow_wid3 - arrow_wid < 40:
        arrow_wid2 = (arrow_wid + arrow_wid3) /2
    else:
        arrow_wid2 = arrow_wid3 - 20
    letter_x = (arrow_wid + arrow_wid3) / 2
    main_html.write('<polygon points="')
    main_html.write("%d,10 %d,40 %d,40 %d,50 %d,25 %d,0 %d,10 %d,10" % (arrow_wid, arrow_wid, arrow_wid2, arrow_wid2, arrow_wid3, arrow_wid2, arrow_wid2, arrow_wid))
    main_html.write('" style="fill:' + fill_color + ';stroke:black;stroke-width:.5" />' )
    main_html.write('<text x="%d"y="32" font-family="sans-serif" font-size="12px" text-anchor="middle" fill="grey">' % (letter_x))
    main_html.write(str(index))
    main_html.write("</text>")

def draw_orf_diagram(main_html, record):
    main_html.write('<h3>Architecture</h3>\n')
    main_html.write('<svg width="1060" height="53">')
    bsc_start = record.CDSs[0].start
    bsc_end = record.CDSs[-1].end
    sub_by = bsc_start - 500
    scale_factor = (660./(bsc_end - bsc_start))
    for cds in record.CDSs:
        draw_CDS_arrow(main_html, cds, sub_by, scale_factor)
    main_html.write('</svg>')
    main_html.write('<svg width="1060" height="53">')
    index = 0
    for orf in record.intergenic_orfs:
        index += 1
        draw_orf_arrow(main_html, orf, sub_by, scale_factor, index)
    main_html.write('</svg>')
    bar_length = scale_factor * 1000
    bar_legx = bar_length + 5
    main_html.write('<svg width="500" height="23">')
    main_html.write('<polygon points="')
    main_html.write("0,10 %f, 10" % (bar_length))
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('<text x="%f"y="12"' % (bar_legx))
    main_html.write("""font-famil="sans-serif"
                    font-size="10px"
                    text_anchor="right"
                    fill="black">1000 nucleotides</text>""")
    main_html.write('<polygon points="')
    main_html.write("0,5 0,15")
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('<polygon points="%f,5 %f,15' % (bar_length, bar_length))
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('</svg>')
    
def draw_cds_table(main_html, record):
    main_html.write("""<br><br><table class="table table-condensed">
  <tbody>
    <tr>
      <th scope="col">Accession</th>
      <th scope="col">start</th>
      <th scope="col">end</th>
      <th scope="col">direction</th>
      <th scope="col">length (aa)</th>
      <th scope="col">Pfam/HMM</th>
      <th scope="col">E-value</th>    
      <th scope="col">description</th>
    </tr>""")
    for cds in record.CDSs:
        main_html.write("<tr>\n")
        main_html.write("""\t<td><a href='https://www.ncbi.nlm.nih.gov/protein/%s'>%s</a></td>
            <td>%s</td> 
            <td>%d</td>
            <td>%s</td>
            <td>%d</td>""" % (cds.accession_id, cds.accession_id, cds.start, cds.end, cds.direction, len(cds.sequence)))
        if len(cds.pfam_descr_list) == 0:
            main_html.write("<td>NO PFAM MATCH</td>")
        else:
            main_html.write("<td><a href='http://pfam.xfam.org/family/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            for pfamid, _, _, in cds.pfam_descr_list[1:]:
                main_html.write("<br><a href='http://pfam.xfam.org/family/%s'>%s</a>" % (pfamid, pfamid))
            e_val = cds.pfam_descr_list[0][2]
            main_html.write("</td><td>%.2E" % Decimal(e_val))
            for _, _, e_val, in cds.pfam_descr_list[1:]:
                main_html.write("<br>%.2E" % Decimal(e_val))
            descr = cds.pfam_descr_list[0][1]
            main_html.write("</td><td>%s" % (descr))
            for _, descr, _, in cds.pfam_descr_list[1:]:
                main_html.write("<br>%s" % (descr))
            main_html.write("</td>")
        main_html.write("</tr>") 
    main_html.write("</tbody></table><p></p>")
       

def draw_orf_table(main_html, record, peptide_type):
    if peptide_type in ["lasso", "lanthi", "sacti", "thio"]:
        main_html.write("""<table class="table table-bordered">
          <tbody>
            <tr>
              <th scope="col">leader</th>
              <th scope="col">core</th>""")
    else:
        main_html.write("""<table class="table table-bordered">
      <tbody>
        <tr>
          <th scope="col">peptide</th>""")
    main_html.write("""
      <th scope="col">start</th>
      <th scope="col">end</th>
      <th scope="col">dir</th>
      <th scope="col">score</th>
    </tr>""")
    for ripp in record.ripps[peptide_type]:
        main_html.write("<tr>\n")
        if peptide_type in ["lasso", "lanthi", "sacti", "thio"]:
            main_html.write("<td>%s</td>" % (ripp.leader))
            main_html.write("<td>%s</td>" % (ripp.core))
        else:
            main_html.write("<td>%s</td>" % (ripp.sequence))
        main_html.write("<td>%d</td>" % (ripp.start))
        main_html.write("<td>%d</td>" % (ripp.end))
        main_html.write("<td>%s</td>" % (ripp.direction))
        main_html.write("<td>%d</td>" % (ripp.score))
        main_html.write("</tr>\n")
    main_html.write("</tbody></table>")
        

def write_table_of_contents(main_html, queries):
    main_html.write("<h3> Input Queries (click to navigate)</h3>")
    main_html.write('<ul style="list-style-type:none">')
    for query in queries:
        main_html.write('<li><a href="#%s">%s</a></li>' % (query, query))
    main_html.write("</ul>\n")

def write_failed_query(main_html, query, message):
    main_html.write('<h2 id="%s"> Results for %s\n' % (query, query))
    main_html.write('<a href="#header"><small><small>back to top</small></small></a></h2>') #TODO keep for single?
    main_html.write('<p></p>') # TODO why
    main_html.write(message)
    main_html.write('<p></p>')
    
def write_record(main_html, record, peptide_type):
    #RESULTS FOR xxxx
    #DRAW ORF
    #DRAW ORF scale
    #PUT LINK TO nuc SEQUENCE
    #TABLE of CDS
    #TABLE of ORFs
    
    main_html.write('<h2 id="%s"> Results for %s [%s]\n' % (record.query_accession_id, record.query_accession_id, record.cluster_genus_species))
    main_html.write('<a href="#header"><small><small>back to top</small></small></a></h2>') #TODO keep for single?
    main_html.write('<p></p>') # TODO why
    draw_orf_diagram(main_html, record)
    main_html.write('<p></p>') 
    main_html.write('<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s">Link to nucleotide sequence</a>' % (record.cluster_accession))
    draw_cds_table(main_html, record)
    draw_orf_table(main_html, record, peptide_type)
    