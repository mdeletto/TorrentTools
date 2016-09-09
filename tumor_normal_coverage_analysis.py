#!/usr/bin/env python

import optparse
import re
import sys
import subprocess
import os
import math
from collections import defaultdict
from biomart import BiomartServer

# sys.path.append("/home/michael/YNHH/Code/Bioinformatics/Globals")
# from global_variables import *

global BEDTOOLS_EXE
global opts

BEDTOOLS_EXE = "/home/michael/bin/bedtools2/bin/bedtools"

#-------------------------------------------------------------------------------------------
#----------------------------Command line parser and arguments------------------------------
#-------------------------------------------------------------------------------------------

desc="""Calculates base depth of coverage for each base across a supplied set of target regions, and also calculates mean depth across each region."""

parser = optparse.OptionParser(description=desc)

parser.add_option("-b", action="store", help="Base output",dest="base_output")
parser.add_option("-d", action="store", help="Output directory name (local only).  Default is to use -b option.",dest="dir_output")
parser.add_option('-n','--normal_bam', help="Path to normal bam", dest='normal_bam', action='store', default=None)
parser.add_option('-t','--tumor_bam', help="Path to tumor bam", dest='tumor_bam', action='store')
parser.add_option('-r','--target_region', help="Target regions in BED format", dest='target_bed', action='store',default="/home/michael/YNHH/Reference_Files/OCP/exons/OCP_exons.coding.sorted.ensembl.bed")
parser.add_option('--target_type', help="Select <Exon> or <Amplicon> for determination of target region type.  Please consult README for acceptable inputs.)", dest='target_type', action='store',default="Exon")
parser.add_option('--min_tumor_depth', help="Minimum tumor depth for calculations", dest='min_tumor_depth', action='store',default=str(20))
parser.add_option('--min_normal_depth', help="Minimum normal depth for calculations", dest='min_normal_depth', action='store',default=str(5))
#parser.add_option("--merged_bed", action="store", help="Flag indicating ",dest="base_output")


(opts, args) = parser.parse_args()

if opts.normal_bam is None:
    parser.error("ERROR: Option is unsupported currently.")
    parser.print_help()

def bedtools_mean_depth_command(target_bed,bam_file,base_output,sample_type):
    subprocess.call("%s coverage -sorted -g /home/michael/YNHH/Reference_Files/FASTA/hg19.genome -mean -a %s -b %s > %s-%s.coverage_statistics.mean_depth_per_region.tsv" % (BEDTOOLS_EXE,target_bed,bam_file,base_output,sample_type),shell=True)

def bedtools_per_base_depth_command(target_bed,bam_file,base_output,sample_type):
    subprocess.call("%s coverage -sorted -g /home/michael/YNHH/Reference_Files/FASTA/hg19.genome -d -a %s -b %s > %s-%s.coverage_statistics.coverage_per_base.tsv" % (BEDTOOLS_EXE,target_bed,bam_file,base_output,sample_type),shell=True)

def awk_add_last_column(base_output):
    subprocess.call("""awk '{FS="\t"; getline f1 <"%s-normal.coverage_statistics.coverage_per_base.tsv" ; print f1 "\t" $(NF)}' < %s-tumor.coverage_statistics.coverage_per_base.tsv > %s-normal-tumor.coverage_statistics.coverage_per_base.tsv""" % (base_output,base_output,base_output) ,shell=True)
    subprocess.call("""awk '{FS="\t"; getline f1 <"%s-normal.coverage_statistics.mean_depth_per_region.tsv" ; print f1 "\t" $(NF)}' < %s-tumor.coverage_statistics.mean_depth_per_region.tsv > %s-normal-tumor.coverage_statistics.mean_depth_per_region.tsv""" % (base_output,base_output,base_output) ,shell=True)
    coverage_files = [
                      "%s-normal-tumor.coverage_statistics.coverage_per_base.tsv" % base_output,
                      "%s-normal-tumor.coverage_statistics.mean_depth_per_region.tsv" % base_output
                      ]
    return coverage_files

def split_based_on_gene_category(coverage_files):
    
    class BED_fields:
        def __init__(self,line):

            strip_line = line.strip()
            split_line = strip_line.split("\t")            
            
            if opts.target_type == "Exon":
                self.gene = split_line[3]
            elif opts.target_type == "Amplicon":
                match = re.search("GENE_ID=(.\w+)",line)
                self.gene = match.group(1)
            else:
                pass

    for file in coverage_files:
        tmp_name = str(file).replace(".tsv","")
        hotspot_output = open("%s.hotspots.tsv" % tmp_name,"w")
        tumor_suppressor_output = open("%s.tumor_suppressors.tsv" % tmp_name,"w")
        CNV_output = open("%s.amplified.tsv" % tmp_name,"w")

        with open(file,"rw") as f:
            counter = 0
            for line in f.readlines():
                counter += 1
#                 if counter == 1:
#                         hotspot_output.write(line)
#                         tumor_suppressor_output.write(line)
#                         CNV_output.write(line)
#                 else:
                strip_line = line.strip()
                split_line = strip_line.split("\t")
                bed_fields = BED_fields(line)
                if bed_fields.gene in OCP_HOTSPOT_GENES:
                    hotspot_output.write(line)
                elif bed_fields.gene in OCP_TUMOR_SUPPRESSOR_GENES:
                    tumor_suppressor_output.write(line)
                elif bed_fields.gene in OCP_AMPLIFIED_GENES:
                    CNV_output.write(line)
                else:
                    pass
            hotspot_output.close()
            tumor_suppressor_output.close()
            CNV_output.close()
            
        if re.search("coverage_per_base",file):
            calculate_uniformity(file)

def calculate_uniformity(input_coverage_per_base):

    #---------UNIFORMITY CALCULATION----------
    
    def perform_calculations(chrom_loc_depth_dict,sample_type):
    
        uniformity_counter = 0
        read_depth_sum = 0    
    
        for key,value in chrom_loc_depth_dict.iteritems():
            read_depth_sum += value
           
        mean_read_depth = read_depth_sum/len(chrom_loc_depth_dict) # mean read depth = total depth / total bp
           
        for chrom_loc in chrom_loc_depth_dict:
            if chrom_loc_depth_dict[chrom_loc] >= (0.2 * mean_read_depth): # calculate number of bases that meet uniformity threshold (0.2x mean read depth)
                uniformity_counter += 1                            # count bases      
            else:
                pass
               
        uniformity = (float(uniformity_counter)/float(len(chrom_loc_depth_dict)))*100   # calculate uniformity = bases at 0.2x mean read depth / total evaluated bases
        print "%s: Uniformity = %s" % (sample_type,float(uniformity))

    chrom_loc_depth_dict_normal = {}
    chrom_loc_depth_dict_tumor = {}
    chrom_loc_depth_per_pool = []
    pools = []
    pool_read_depth_dict={}
    
    f = open(input_coverage_per_base,"r")
    
    for line in f.readlines():
        line = line.rstrip().split("\t")
    
        class fields:
            chrom_loc = "%s\t%s" % (str(line[0]),int(line[1])+int(line[-3])-1)
            normal_depth = int(line[-2])
            tumor_depth = int(line[-1])
            pool = line[4]
            
        pools.append(fields.pool) # add all pool numbers to a list (which will be reduced to a set later)
        chrom_loc_depth_dict_normal[fields.chrom_loc] = fields.normal_depth  # unique locations matched to depth --- used to calculate total uniformity
        chrom_loc_depth_dict_tumor[fields.chrom_loc] = fields.tumor_depth  # unique locations matched to depth --- used to calculate total uniformity
#         chrom_loc_depth_per_pool.append([fields.chrom_loc,fields.depth,fields.pool])  # a list of all locations with depth and pool number --- used to calculate uniformity per pool
    
    perform_calculations(chrom_loc_depth_dict_normal,"normal")
    perform_calculations(chrom_loc_depth_dict_tumor,"tumor")
    

    
    
    #----------POOL-SPECIFIC UNIFORMITY CALCULATION-----------
    
#     for pool in set(pools):
#         read_depth_sum_per_pool = 0
#         base_number_sum = 0
#         pool_read_depth_dict[pool] = [read_depth_sum_per_pool,base_number_sum]
#      
#     for chrom_loc in chrom_loc_depth_per_pool:
#         pool = chrom_loc[2]
#         depth = chrom_loc[1]
#         pool_read_depth_dict[pool][0] += depth
#         pool_read_depth_dict[pool][1] += 1

#     for key,value in pool_read_depth_dict.iteritems():
#         uniformity_counter = 0
#         for chrom_loc in chrom_loc_depth_per_pool:
#             pool_number = key
#             read_depth_sum = value[0]
#             base_number_sum = value[1]
#             depth_at_pos = chrom_loc[1]
#             pool = chrom_loc[2]
#             if pool_number == pool and depth_at_pos >= (0.2*(read_depth_sum/base_number_sum)):
#                 uniformity_counter += 1
#         uniformity_per_pool = (float(uniformity_counter)/float(base_number_sum))*100


def plot_low_coverage_regions(coverage_files):
    
    def parse_minimum_requested_depth_thresholds(values):
        values = values.split(",")
        return values

    tumor_depth_values = parse_minimum_requested_depth_thresholds(opts.min_tumor_depth)
    normal_depth_values = parse_minimum_requested_depth_thresholds(opts.min_normal_depth)

    if len(tumor_depth_values) != len(normal_depth_values):
        print "ERROR: Length of provided min_tumor_depths for graphing does not match length of min_normal_depths"
        print parser.print_help()
        sys.exit(1)
    else:

        
        mapped_tumor_normal_values = zip(normal_depth_values,tumor_depth_values)

        read_depth_filter_criteria = {
                                      "YCGA_Depth_Cutoff" : [4,6],
                                      "TPL_Depth_Cutoff" : [5,20]
                                      }
        
        for file in coverage_files:
            if re.search("mean_depth_per_region",file):
                for value_pair in mapped_tumor_normal_values:
                    
                    normal_coverage_threshold, tumor_coverage_threshold = (i for i in value_pair)
                    if normal_coverage_threshold == str(5) and tumor_coverage_threshold == str(20):
                        filter_name = "TPL_Depth_Cutoff"
                    elif normal_coverage_threshold == str(4) and tumor_coverage_threshold == str(6):
                        filter_name = "YCGA_Depth_Cutoff"
                    else:
                        filter_name = "Custom_Depth_Cutoff"
                    
    #             for filter_name,depth_values in read_depth_filter_criteria.iteritems():
    #                 #tmp_name = str(file).replace(".tsv","")
    #                 normal_coverage_threshold, tumor_coverage_threshold = (i for i in depth_values)  
                    subprocess.call("""/usr/bin/Rscript /home/michael/YNHH/Code/Bioinformatics/R/coverage_plot_from_bedtools_coverage.R %s %s %s %s %s %s""" % (file,base_output,filter_name,normal_coverage_threshold,tumor_coverage_threshold, target_type),shell=True)

def create_empty_normal_column_in_file(base_output):
    subprocess.call("""awk '{FS="\t"; printf $1FS$2FS$3FS$4FS$5FS$6FS$7FS0FS$8; print NL}' %s-tumor.coverage_statistics.coverage_per_base.tsv > %s-normal-tumor.coverage_statistics.coverage_per_base.tsv""" % (base_output, base_output), shell=True)         
    subprocess.call("""awk '{FS="\t"; printf $1FS$2FS$3FS$4FS$5FS$6FS0FS$7; print NL}' %s-tumor.coverage_statistics.mean_depth_per_region.tsv > %s-normal-tumor.coverage_statistics.mean_depth_per_region.tsv""" % (base_output, base_output), shell=True)         
    coverage_files = [
                      "%s-normal-tumor.coverage_statistics.coverage_per_base.tsv" % base_output,
                      "%s-normal-tumor.coverage_statistics.mean_depth_per_region.tsv" % base_output
                      ]
    return coverage_files

def convert_coordinate_to_protein(genomic_coding_start, genomic_coding_end, aa_start, strand, coordinate_overlap_float, split_query, ensembl_gene_response):
    coordinate_to_protein = defaultdict(lambda: None)
    aa_start = int(aa_start)
    aa_frame_counter = 1
    if strand == "1" or strand == 1:
        counter = int(genomic_coding_start)
        while counter <= int(genomic_coding_end):
            if aa_frame_counter <= 3:
                coordinate_to_protein[counter] = int(aa_start)
                aa_frame_counter += 1
            else:
                aa_frame_counter = 1
                aa_start += 1
                coordinate_to_protein[counter] = int(aa_start)
                aa_frame_counter += 1
            counter += 1
    elif strand == "-1" or strand == -1:
        counter = int(genomic_coding_end)
        while counter >= int(genomic_coding_start):
            if aa_frame_counter <= 3:
                coordinate_to_protein[counter] = int(aa_start)
                aa_frame_counter += 1
            else:
                aa_frame_counter = 1
                aa_start += 1
                coordinate_to_protein[counter] = int(aa_start)
                aa_frame_counter += 1
            counter -= 1
            
        
    if coordinate_overlap_float == 1:
        aa_start = ensembl_gene_response.aa_start
        aa_end = ensembl_gene_response.aa_end
    else:
        if int(split_query.start) in coordinate_to_protein.keys():
            if strand == 1 or strand == "1":
                aa_start = coordinate_to_protein[int(split_query.start)] 
                #aa_start = coordinate_to_protein[split_query.start]
            elif strand == -1 or strand == "-1":
                aa_start = coordinate_to_protein[int(split_query.start)]
                #aa_start = coordinate_to_protein(split_query.end)
        else:
            if strand == 1 or strand == "1":
                aa_start = ensembl_gene_response.aa_start
            elif strand == -1 or strand == "-1":
                aa_start = ensembl_gene_response.aa_end
                
        if int(split_query.end) in coordinate_to_protein.keys():
            if strand == 1 or strand == "1":
                aa_end = coordinate_to_protein[int(split_query.end)]
                #aa_end = coordinate_to_protein[split_query.end]
            elif strand == -1 or strand == "-1":
                aa_end = coordinate_to_protein[int(split_query.end)]
                #aa_end = coordinate_to_protein[split_query.start]
        else:
            if strand == 1 or strand == "1":
                aa_end = ensembl_gene_response.aa_end
            elif strand == 1 or strand == "-1":
                aa_end = ensembl_gene_response.aa_start
    
    if int(aa_start) > int(aa_end):
        tmp_list = [aa_end, aa_start]
    else:
        tmp_list = [aa_start, aa_end]
    
    
    return tmp_list

def open_low_target_dropout_file(base_output):

    headers = ['#chr:start-end',
               'Normal_Mean_Depth',
               'Tumor_Mean_Depth',
               'HGNC_Gene_Symbol',
               'ENSEMBL_Gene_ID',
               'ENSEMBL_Transcript_ID',
               'Transcript::CDS_Length',
               'Transcript::CDS_Start',
               'Transcript::CDS_End',
               'Transcript::AA_Start',
               'Transcript::AA_End',
               'ENSEMBL_Exon_ID',
               'Exon::Rank',
               'Exon::Start_Coordinate',
               'Exon::End_Coordinate',
               'Exon::Genomic_Coding_Start_Coordinate',
               'Exon::Genomic_Coding_End_Coordinate',
               'Query::AA_Start',
               'Query:AA_End',
               'Query::Overlap_with_Exon_(bp)',
               'Query::Percent_of_Exon_Overlap']
    
    with open("%s/%s-normal-tumor.coverage_statistics.low_coverage_regions.tsv" % (os.getcwd(), base_output), "w") as f:
        f.write("\t".join(headers) + "\n") 
        f.close()
    #print 

def query_biomart_and_collect_overlapping_gene_information(query):

    class Query:
        def __init__(self,query):
            split_query = query.split(":")
            self.chr = split_query[0]
            self.start = split_query[1]
            self.end = split_query[2]
            
    split_query = Query(query)

    response = ensembl_genes.search({
        'filters' : {
            'chromosomal_region' : '%s' % query
        },
        'attributes' : [
            'ensembl_gene_id',
            'external_gene_name',
            'ensembl_transcript_id',
            'cds_start',
            'cds_end',
            'cds_length',
            'rank',
            'cdna_coding_start',
            'cdna_coding_end'
        ]
    })                            
    
    
    for line in response.iter_lines():
        line = line.decode('utf-8')
        split_line = line.strip().split("\t")
        
        class EnsemblGenesResponse():
            def __init__(self, split_line):
                self.ensembl_gene_id = split_line[0]
                self.gene_name = split_line[1]
                self.ensembl_transcript_id = split_line[2]
                self.cds_start = split_line[3]
                self.cds_end = split_line[4]
                self.cds_length = split_line[5]
                self.exon_rank = split_line[6]
        
        
        ensembl_gene_response = EnsemblGenesResponse(split_line)
        if not ensembl_gene_response.cds_length == "":
            if ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] is None:
                ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] = {
                                                                                              'transcript' : ensembl_gene_response.ensembl_transcript_id,
                                                                                              'cds_length' : ensembl_gene_response.cds_length
                                                                                              }
            else:
                
                if ensembl_gene_response.gene_name in ensembl_gene_name_to_canonical_transcript.keys():
                    try:
                        if ensembl_gene_response.ensembl_transcript_id == ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name]['transcript']:
                            ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] = {
                                                                                                         'transcript' : ensembl_gene_response.ensembl_transcript_id,
                                                                                                          'cds_length' : ensembl_gene_response.cds_length
                                                                                                          }
                        else:
                            if ensembl_gene_response.ensembl_transcript_id != ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name]['transcript']:
                                if int(ensembl_gene_response.cds_length) > int(ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name]['cds_length']):
                                    ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] = {
                                                                                                                  'transcript' : ensembl_gene_response.ensembl_transcript_id,
                                                                                                                  'cds_length' : ensembl_gene_response.cds_length
                                                                                                                  }   
                    except KeyError:
                        if ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] == ensembl_gene_response.ensembl_transcript_id:
                            ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] = {
                                                                                                         'transcript' : ensembl_gene_response.ensembl_transcript_id,
                                                                                                          'cds_length' : ensembl_gene_response.cds_length
                                                                                                          }
                        else:
                            pass
        else:
            pass
            # Remove genes with canonical transcripts that may have been added in the first step
            # ONLY THOSE GENES THAT HAVE A CDS LENGTH THAT IS BLANK.  This corresponds to pseudogenes
#             if ensembl_gene_response.gene_name in ensembl_gene_name_to_canonical_transcript.keys():
#                 ensembl_gene_name_to_canonical_transcript.pop(ensembl_gene_response.gene_name)

    
    response = ensembl_genes.search({
        'filters' : {
            'chromosomal_region' : '%s' % query
        },
        'attributes' : [
            'ensembl_gene_id',
            'external_gene_name',
            'ensembl_transcript_id',
            'cds_start',
            'cds_end',
            'cds_length',
            'rank',
            'cdna_coding_start',
            'cdna_coding_end',
            'exon_chrom_start',
            'exon_chrom_end',
            'ensembl_exon_id',
            'genomic_coding_start',
            'genomic_coding_end',
            'strand'
        ]
    })       

    printed_bool = False
    
    low_coverage_annotated_output = open("%s/%s-normal-tumor.coverage_statistics.low_coverage_regions.tsv" % (os.getcwd(), base_output), "a")

    for line in response.iter_lines():
        line = line.decode('utf-8')
        split_line = line.strip().split("\t")
        class EnsemblGenesResponse():
            def __init__(self, split_line):
                self.ensembl_gene_id = split_line[0]
                self.gene_name = split_line[1]
                self.ensembl_transcript_id = split_line[2]
                self.cds_start = split_line[3]
                self.cds_end = split_line[4]
                self.cds_length = split_line[5]
                self.exon_rank = split_line[6]
                try:
                    self.cdna_coding_start = split_line[7]
                except:
                    self.cdna_coding_start = ''
                try:
                    self.cdna_coding_end = split_line[8]
                except:
                    self.cdna_coding_end = ""
                if self.cds_start == '':
                    self.aa_start = ''
                else:
                    self.aa_start = str( int( math.ceil(float(self.cds_start) / 3) )  )
                if self.cds_end == '':
                    self.aa_end = ''
                else:
                    self.aa_end = str( int( math.ceil(float(self.cds_end) / 3) )  )                   
                self.exon_start = split_line[9]
                self.exon_end = split_line[10]
                self.exon_id = split_line[11]
                try:
                    self.genomic_coding_start = split_line[12]
                except:
                    self.genomic_coding_start = ""
                try:
                    self.genomic_coding_end = split_line[13]
                except:
                    self.genomic_coding_end = ""
                self.strand = split_line[14]
                
                
        ensembl_gene_response = EnsemblGenesResponse(split_line)


        if ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name] is None:
            pass
        else:
            if ensembl_gene_name_to_canonical_transcript[ensembl_gene_response.gene_name]['transcript'] == ensembl_gene_response.ensembl_transcript_id and printed_bool is False:
                
                
                gene_attributes_to_print = [
                                            '%s:%s-%s' % (split_query.chr, split_query.start, split_query.end),
                                            round(intervals_with_depths[query]['normal_mean_read_depth'], 1),
                                            round(intervals_with_depths[query]['tumor_mean_read_depth'], 1),
                                            ensembl_gene_response.gene_name,
                                            ensembl_gene_response.ensembl_gene_id,
                                            ensembl_gene_response.ensembl_transcript_id,
                                            ensembl_gene_response.cds_length,
                                            ensembl_gene_response.cds_start,
                                            ensembl_gene_response.cds_end,
                                            ensembl_gene_response.aa_start,
                                            ensembl_gene_response.aa_end,
                                            ensembl_gene_response.exon_id,
                                            ensembl_gene_response.exon_rank,
                                            ensembl_gene_response.exon_start,
                                            ensembl_gene_response.exon_end,
                                            ensembl_gene_response.genomic_coding_start,
                                            ensembl_gene_response.genomic_coding_end
                                            ]
                try:
                    #if int(ensembl_gene_response.genomic_coding_start) <= int(split_query.end) and int(split_query.start) <= int(ensembl_gene_response.genomic_coding_end):
                    if int(ensembl_gene_response.exon_start) <= int(split_query.end) and int(split_query.start) <= int(ensembl_gene_response.exon_end):    
                        def getOverlap(a, b):
                            return max(0, min(a[1], b[1]) - max(a[0], b[0]))
                            
                        # set boundaries
                        
                        exon_boundary = [int(ensembl_gene_response.exon_start), int(ensembl_gene_response.exon_end)]
                        query_boundary = [int(split_query.start), int(split_query.end)]
                        try:
                            genomic_coding_boundary = [int(ensembl_gene_response.genomic_coding_start), int(ensembl_gene_response.genomic_coding_end)]
                        except:
                            genomic_coding_boundary = None
                        
                        # calculate length
                        
                        query_length = int(split_query.end) - int(split_query.start) + 1
                        exon_length = int(ensembl_gene_response.exon_end) - int(ensembl_gene_response.exon_start) + 1
                        try:
                            genomic_coding_length = int(ensembl_gene_response.genomic_coding_end) - int(ensembl_gene_response.genomic_coding_start) + 1
                        except:
                            genomic_coding_length = None
                        
                        # calculate ranges
                        
                        exon_range = range(int(ensembl_gene_response.exon_start), int(ensembl_gene_response.exon_end))
                        query_range = range(int(split_query.start), int(split_query.end))
                        try:
                            genomic_coding_range = range(int(ensembl_gene_response.genomic_coding_start), int(ensembl_gene_response.genomic_coding_end))
                        except:
                            genomic_coding_range = None

                        xs = set(query_range)
                        
                        # calculate overlap between query and other boundary
                        
                        coordinate_overlap_int = getOverlap(exon_boundary, query_boundary) + 1
                        coordinate_overlap_float = float(coordinate_overlap_int) / float (exon_length) 
                        
                        #print coordinate_overlap_float
                        #print min_overlap_coordinate, max_overlap_coordinate + 1

                        # calculate aa_start and aa_end for query based on genomic coding coordinates and aa positions
                        try:
                            if int(ensembl_gene_response.genomic_coding_start) <= int(split_query.end) and int(split_query.start) <= int(ensembl_gene_response.genomic_coding_end):
                                
                                # find min and max overlap coordinate for each genomic coding range
                                
                                min_overlap_coordinate = min(sorted(xs.intersection(genomic_coding_range)))
                                max_overlap_coordinate = max(sorted(xs.intersection(genomic_coding_range)))

                                
                                aa_start, aa_end = (i for i in convert_coordinate_to_protein(ensembl_gene_response.genomic_coding_start, ensembl_gene_response.genomic_coding_end, ensembl_gene_response.aa_start, ensembl_gene_response.strand, coordinate_overlap_float, split_query, ensembl_gene_response))
                                gene_attributes_to_print.append(aa_start)
                                gene_attributes_to_print.append(aa_end)

                        except:
                            pass
                        
                        # append additional values to print
                        
                        gene_attributes_to_print.append(str(coordinate_overlap_int))
                        gene_attributes_to_print.append(str(round(coordinate_overlap_float, 2)))
                        
                        # print out values
                        
                        low_coverage_annotated_output.write("\t".join(str(i) for i in gene_attributes_to_print) + "\n")
                        printed_bool = True
                
                    #elif int(ensembl_gene_response.chromosome_start) <= int(split_query.end) and int(split_query.start) <= int(ensembl_gene_response.chromosome_end):
                    else:
                        gene_attributes_to_print_with_no_exon_overlap = [
                                                                         '%s:%s-%s' % (split_query.chr, split_query.start, split_query.end),
                                                                        round(intervals_with_depths[query]['normal_mean_read_depth'], 1),
                                                                        round(intervals_with_depths[query]['tumor_mean_read_depth'], 1),
                                                                        ensembl_gene_response.gene_name,
                                                                        ensembl_gene_response.ensembl_gene_id,
                                                                        ensembl_gene_response.ensembl_transcript_id,
                                                                        ensembl_gene_response.cds_length]
                        
                        #low_coverage_annotated_output.write("\t".join(str(i) for i in gene_attributes_to_print) + "\n")
                        #printed_bool = True
                
                except Exception, e:

                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)
                    print str(e)


    if printed_bool is False:
        low_coverage_annotated_output.write("\t".join(str(i) for i in gene_attributes_to_print_with_no_exon_overlap) + "\n")

def parse_canonical_ensembl_transcripts():

    ensembl_gene_name_to_canonical_transcript = defaultdict(lambda: None)
    
    with open("/home/michael/YNHH/Reference_Files/transcripts/canonical_ensembl_transcripts.txt") as canonical_ensembl_transcripts:
        for line in canonical_ensembl_transcripts:
            split_line = line.strip().split("\t")
            ensembl_transcript = split_line[5]
            gene_name = split_line[6]
              
            if re.search("ENST", ensembl_transcript):
                ensembl_gene_name_to_canonical_transcript[gene_name] = {'transcript' : ensembl_transcript}
                
    return ensembl_gene_name_to_canonical_transcript

def initialize_biomart():

    server = BiomartServer( "http://grch37.ensembl.org/biomart/martservice" )
    ensembl_genes = server.datasets['hsapiens_gene_ensembl']
    
    #server.verbose = True
    #server.show_databases()
    #server.show_datasets()
    #ensembl_genes.show_filters()
    #ensembl_genes.show_attributes()
    
    return ensembl_genes

def calculate_dropouts():

    def calculate_neighbors(iterable):
        iterator = iter(iterable)
        prev = None
        item = iterator.next()  # throws StopIteration if empty.
        for next in iterator:
            yield (prev,item,next)
            prev = item
            item = next
        yield (prev,item,None)

    # This is the value we are returning.
    # Calculate intervals with low or no coverage, annotate with read depths, and return.
    intervals_with_depths = defaultdict(dict)

    
    with open("%s/%s-normal-tumor.coverage_statistics.coverage_per_base.tsv" % (os.getcwd(), base_output)) as f:
        
        positions_below_depth_threshold = defaultdict(dict)
        positions_above_depth_threshold = defaultdict(dict)
        
        
        for line in f.readlines():
            split_line = line.strip().split("\t")
            
            class NormalTumorCoverageOutput:
                def __init__(self, split_line):
                    self.chr = split_line[0].lstrip("chr")
                    self.start_pos = split_line[1]
                    self.end_pos = split_line[2]
                    
                    # Commented out because we don't need these.  If let in, we would need strict BED format compliance
                    #self.amplicon_id = split_line[3]
                    #self.score = split_line[4]
                    #self.description = split_line[5]
                    
                    self.base_no = split_line[-3]
                    self.normal_depth = int(split_line[-2])
                    self.tumor_depth = int(split_line[-1])
                    self.corrected_start_pos = int(self.start_pos) + (int(self.base_no) - 1)
                    self.corrected_end_pos = int(self.start_pos) + int(self.base_no)
            
            
            split_line = NormalTumorCoverageOutput(split_line)
            
            if split_line.tumor_depth < 20 or split_line.normal_depth < 5:
                if "%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos)) in positions_above_depth_threshold:
                    pass
                else:
                    #positions_below_depth_threshold.append("%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos)))
                    positions_below_depth_threshold["%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos))]['normal_depth'] = split_line.normal_depth
                    positions_below_depth_threshold["%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos))]['tumor_depth'] = split_line.tumor_depth
            else:
                positions_above_depth_threshold["%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos))]['normal_depth'] = split_line.normal_depth
                positions_above_depth_threshold["%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos))]['tumor_depth'] = split_line.tumor_depth
                #positions_above_depth_threshold.append("%s_%s_%s" % (split_line.chr, str(split_line.corrected_start_pos), str(split_line.corrected_end_pos)))


    for prev, item, next in calculate_neighbors(sorted(positions_below_depth_threshold.keys())):
        
        class ChrStartEnd:
            def __init__(self, item, position_below_depth_threshold):
                split_item = item.split("_")
                self.chr = split_item[0]
                self.start = str(split_item[1])
                self.end = str(split_item[2])
                try:
                    self.normal_depth = position_below_depth_threshold[item]['normal_depth']
                except KeyError:
                    self.normal_depth = 0
                try:
                    self.tumor_depth = position_below_depth_threshold[item]['tumor_depth']
                except KeyError:
                    self.tumor_depth = 0
        
        item = ChrStartEnd(item, positions_below_depth_threshold)
        
        # If prev is None, this is the first base we encounter EVER
        # Initialize all counters wit hthe first entry
        if prev is None:
            chr = item.chr
            start = item.start
            end = item.end
            entry = "%s:%s:%s" % (chr, start, end)
            base_counter = 1
            normal_depth_counter = item.normal_depth
            tumor_depth_counter = item.tumor_depth
            
        # If next is not None, we can use it to see if the next item is in the interval
        if next is not None:
            
            next = ChrStartEnd(next, positions_above_depth_threshold)
            
            # If next.start == item.end, then we can continue the interval
            if next.start == item.end:
                end = next.end
                entry = "%s:%s:%s" % (chr, start, end)
                base_counter += 1
                normal_depth_counter += item.normal_depth
                tumor_depth_counter += item.tumor_depth

                
            # If next.start != item.end, we print the entry and reinitialize counters with the values of the next entry
            else:
                base_counter += 1
                normal_depth_counter += item.normal_depth
                tumor_depth_counter += item.tumor_depth

                intervals_with_depths[entry]['normal_mean_read_depth'] = float(float(normal_depth_counter) / float(base_counter))
                intervals_with_depths[entry]['tumor_mean_read_depth'] = float(float(tumor_depth_counter) / float(base_counter))
                
                # Reset all counters if the next starting base is on the same interval
                
                chr = next.chr
                start = next.start
                end = next.end
                entry = "%s:%s:%s" % (chr, start, end)
                base_counter = 0
                normal_depth_counter = next.normal_depth
                tumor_depth_counter = next.tumor_depth
                
        # If next base is None, we have reached the end of the list.     
        else:
            
            prev = ChrStartEnd(prev, positions_above_depth_threshold)
            
            # If prev.end == item.start, then this current base is part of the interval still
            if prev.end == item.start:
                end = item.end
                entry = "%s:%s:%s" % (chr, start, end)
                base_counter += 1
                normal_depth_counter += item.normal_depth
                intervals_with_depths[entry]['normal_mean_read_depth'] = float(float(normal_depth_counter) / float(base_counter))
                intervals_with_depths[entry]['tumor_mean_read_depth'] = float(float(tumor_depth_counter) / float(base_counter))
            # If prev.end != item.start, then the new interval becomes a single base
            else:
                chr = item.chr
                start = item.start
                end = item.end
                entry = "%s:%s:%s" % (chr, start, end)
                base_counter = 1
                normal_depth_counter += item.normal_depth
                intervals_with_depths[entry]['normal_mean_read_depth'] = float(float(normal_depth_counter) / float(base_counter))
                intervals_with_depths[entry]['tumor_mean_read_depth'] = float(float(tumor_depth_counter) / float(base_counter))


    return intervals_with_depths


def calculate_coverage_statistics():
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    
    for file in files:
        
        depths_to_assess = [1000,500,200,100,50,20,10,5,1,0]

        normal_depth_dict = {}
        tumor_depth_dict = {}
        line_counter = 0
        
        if re.search("normal-tumor",file) and re.search("mean_depth",file) and re.search(opts.base_output, file) and os.stat(file).st_size != 0:
            normal_depth_dict = normal_depth_dict.fromkeys(depths_to_assess,0)
            tumor_depth_dict = tumor_depth_dict.fromkeys(depths_to_assess,0)
            with open(file) as f:
                for line in f.readlines():
                    line_counter += 1
                    split_line = line.split("\t")
                    normal_depth, tumor_depth = split_line[-2], split_line[-1]
                    for depths in depths_to_assess:
                        if float(normal_depth) > float(depths):
                            normal_depth_dict[depths] += 1
                        if float(tumor_depth) > float(depths):
                            tumor_depth_dict[depths] += 1
                            
            print "------------------------"
            print "FILENAME = %s" % file
            print "------------------------"
            print "Coverage\tTumor\tNormal"
            for k1 in sorted(tumor_depth_dict.keys(),reverse=True):
                for k2 in sorted(normal_depth_dict.keys(),reverse=True):
                    if k1 == k2:
                        print "%s\t%s\t%s" % (k2,float(tumor_depth_dict[k2])/float(line_counter),float(normal_depth_dict[k2])/float(line_counter))              

        if re.search("normal-tumor",file) and re.search("coverage_per_base",file) and re.search(opts.base_output, file) and os.stat(file).st_size != 0:
            normal_depth_dict = normal_depth_dict.fromkeys(depths_to_assess,0)
            tumor_depth_dict = tumor_depth_dict.fromkeys(depths_to_assess,0)
            with open(file) as f:
                for line in f.readlines():
                    line_counter += 1
                    split_line = line.split("\t")
                    normal_depth, tumor_depth = split_line[-2], split_line[-1]
                    for depths in depths_to_assess:
                        if float(normal_depth) > float(depths):
                            normal_depth_dict[depths] += 1
                        if float(tumor_depth) > float(depths):
                            tumor_depth_dict[depths] += 1
                            
            print "------------------------"
            print "FILENAME = %s" % file
            print "------------------------"
            print "Coverage\tTumor\tNormal"
            for k1 in sorted(tumor_depth_dict.keys(),reverse=True):
                for k2 in sorted(normal_depth_dict.keys(),reverse=True):
                    if k1 == k2:
                        print "%s\t%s\t%s" % (k2,float(tumor_depth_dict[k2])/float(line_counter),float(normal_depth_dict[k2])/float(line_counter))                 

OCP_HOTSPOT_GENES = ['ABL1','AKT1','ALK','AR','ARAF','BRAF','BTK','CBL','CDK4','CHEK2','CSF1R',
                     'CTNNB1','DDR2','DNMT3A','EGFR','ERBB2','ERBB3','ERBB4','ESR1','EZH2',
                     'FGFR1','FGFR2','FGFR3','FLT3','FOXL2','GATA2','GNA11','GNAQ','GNAS',
                     'HNF1A','HRAS','IDH1','IDH2','IFITM1','IFITM3','JAK1','JAK2','JAK3','KDR',
                     'KIT','KNSTRN','KRAS','MAGOH','MAP2K1','MAP2K2','MAPK1','MAX','MED12','MET',
                     'MLH1','MPL','MTOR','MYD88','NFE2L2','NPM1','NRAS','PAX5','PDGFRA','PIK3CA',
                     'PPP2R1A','PTPN11','RAC1','RAF1','RET','RHEB','RHOA','SF3B1','SMO','SPOP',
                     'SRC','STAT3','U2AF1','XPO1']

OCP_TUMOR_SUPPRESSOR_GENES = ['APC','ATM','BAP1','BRCA1','BRCA2','CDH1','CDKN2A','FBXW7',
                              'GATA3','MSH2','NF1','NF2','NOTCH1','PIK3R1','PTCH1','PTEN',
                              'RB1','SMAD4','SMARCB1','STK11','TET2','TP53','TSC1','TSC2',
                              'VHL','WT1']

OCP_AMPLIFIED_GENES = ['ACVRL1','AKT1','APEX1','AR','ATP11B','BCL2L1','BCL9','BIRC2','BIRC3',
                       'CCND1','CCNE1','CD274','CD44','CDK4','CDK6','CSNK2A1','DCUN1D1',
                       'EGFR','ERBB2','FGFR1','FGFR2','FGFR3','FGFR4','FLT3','GAS6',
                       'IGF1R','IL6','KIT','KRAS','MCL1','MDM2','MDM4','MET','MYC',
                       'MYCL','MYCN','MYO18A','NKX2-1','NKX2-8','PDCD1LG2','PDGFRA','PIK3CA',
                       'PNP','PPARG','RPS6KB1','SOX2','TERT','TIAF1','ZNF217']

normal_bam, tumor_bam, target_bed, base_output, target_type, dir_output = opts.normal_bam, opts.tumor_bam, opts.target_bed, opts.base_output, opts.target_type, opts.dir_output

bams_to_process = {
                   'normal' : normal_bam,
                   'tumor' : tumor_bam
                   }
if normal_bam is None or normal_bam == "None":
    bams_to_process.pop('normal')


for sample_type in bams_to_process.keys():
    bedtools_mean_depth_command(target_bed, bams_to_process[sample_type], base_output, sample_type)
    bedtools_per_base_depth_command(target_bed, bams_to_process[sample_type], base_output, sample_type)

if normal_bam is None or normal_bam == "None":
    sys.exit("ERROR: Option is unsupported")
    #coverage_files = create_empty_normal_column_in_file(base_output)
else:
    coverage_files = awk_add_last_column(base_output)



gene_subtypes = {
                 "hotspots" : OCP_HOTSPOT_GENES,
                 "tumor_suppressors" : OCP_TUMOR_SUPPRESSOR_GENES,
                 "amplified" : OCP_AMPLIFIED_GENES
                 }


split_based_on_gene_category(coverage_files)
#plot_low_coverage_regions(coverage_files)
calculate_coverage_statistics()

ensembl_gene_name_to_canonical_transcript = parse_canonical_ensembl_transcripts()
ensembl_genes = initialize_biomart()

intervals_with_depths = calculate_dropouts() 

open_low_target_dropout_file(base_output)

for interval in sorted(intervals_with_depths.keys()):
    query_biomart_and_collect_overlapping_gene_information(interval)


if dir_output == base_output:
    subprocess.call("mkdir %s; mv %s* %s/" % (base_output, base_output, base_output), shell=True)
else:
    subprocess.call("mkdir %s; mv %s* %s/" % (dir_output, base_output, dir_output), shell=True)

        