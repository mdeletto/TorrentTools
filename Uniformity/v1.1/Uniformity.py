#/usr/bin/python

#from collections import OrderedDict
import sys, os, re, subprocess

def define_regions(regions,gender):
    if regions=="CHPv2":
        print "Unsupported option -- please try again."
        sys.exit(1)
    elif regions=="CCP":
        bed_filepath = "/results/plugins/Uniformity/BEDs/AmpliSeq_CCP.bed"
    elif regions=="OCP":
        bed_filepath = "/results/plugins/Uniformity/BEDs/AmpliSeq_OCP.bed"
    elif regions=="WhEx" and gender=="Male":
        bed_filepath = "/results/plugins/Uniformity/BEDs/AmpliSeq_Exome.male.bed"
    elif regions=="WhEx" and gender=="Female":
        bed_filepath = "/results/plugins/Uniformity/BEDs/AmpliSeq_Exome.female.bed"
    elif regions is None:
        print "No argument passed for defined regions."
        sys.exit(1)
    else:
        print "Unsupported option -- please try again."
        sys.exit(1)
    return bed_filepath

def write_html_header(output):
    header = """
                        <!-- CSS Code -->
                    <style type="text/css" scoped>
                    table.GeneratedTable {
                        width:90%;
                        border-top:1px solid #e5eff8;
                        border-right:1px solid #e5eff8; 
                        border-collapse:collapse;
                        margin:1em auto;
                    }
        
                    a:link {
                        color:#d42945;
                        text-decoration:none;
                        border-bottom:1px dotted #ffbac8;
                    }    
                    a:visited {
                        color:#d42945;
                        border-bottom:none;
                        text-decoration:none;
                    }        
                    a:hover,
                    a:focus {
                        color:#f03b58;
                        border-bottom:1px solid #f03b58;
                        text-decoration:none;
                    }
                    table a,
                    table a:link,
                    table a:visited {
                        border:none;
                    }
            
            
                    tr.odd td    {
                        background:#f7fbff
                    }
            
                        
            
            
                    td {
                        color:#678197;
                        border-bottom:1px solid #e5eff8;
                        border-left:1px solid #e5eff8;
                        padding:.3em 1em;
                        text-align:center;
                    }                
                    th {
                        font-weight:normal;
                        color: #678197;
                        text-align:left;
                        border-bottom: 1px solid #e5eff8;
                        border-left:1px solid #e5eff8;
                        padding:.3em 1em;
                    }                            
                    thead th {
                        background:#f4f9fe;
                        text-align:center;
                        font:bold 1.2em/2em "Century Gothic","Trebuchet MS",Arial,Helvetica,sans-serif;
                        color:#66a3d3
                    }
                    </style>
                    
                    <!-- HTML Code -->
                    <table class="GeneratedTable">
                    <thead>
                    <tr>
                    <th>Barcode</th>
                    <th>Summary</th>
                    </tr>
                    </thead>
                    <tbody>"""
    output.write(header)

def write_html_per_bam(output,results_dir,basename,text):
    html = """
            <tr>
            <td><a href="%s/%s.coverage.tsv">%s</a></td>
            <td>%s</td>
            </tr>""" % (results_dir[17:],basename,basename,text)
    output.write(html)

def write_html_footer(output):
    html = """</tbody>
            </table>"""
    output.write(html)

def determine_base_filename(bam):
    match = re.search("(IonXpress_\d*).*bam",bam)
    basename = match.group(1)
    return basename

def bedtools_coverage(bam_dir,bam,bed,results_dir,basename):
    subprocess.call("bedtools coverage -d -abam %s/%s -b %s > %s/%s.coverage.tsv" % (bam_dir,bam,bed,results_dir,basename),shell=True)

def select_BAMs_for_analysis(BAM_dir):

    IonTorrent_bams = []
    for bam in os.listdir(sys.argv[1]):
            if re.search("Ion.*bam",bam) and bam.endswith(".bam"):
                    IonTorrent_bams.append(bam)
    return IonTorrent_bams


def open_output(output_dir):
    output = open("%s/index_block.html" % output_dir,"w")
    return output


def uniformity_plugin(results_dir,basename):
    
    uniformity_counter = 0
    read_depth_sum = 0
    chrom_loc_depth_dict = {}
    chrom_loc_depth_per_pool = []
    pools = []
    pool_read_depth_dict={}

    f = open("%s/%s.coverage.tsv" % (results_dir,basename),"r")
    
    for line in f.readlines():
        line = line.rstrip().split("\t")
    
        class fields:
            chrom_loc = "%s\t%s" % (str(line[0]),int(line[1])+int(line[5])-1)
            depth = int(line[6])
            pool = line[4]
            
        pools.append(line[4]) # add all pool numbers to a list (which will be reduced to a set later)
        chrom_loc_depth_dict[fields.chrom_loc] = fields.depth  # unique locations matched to depth --- used to calculate total uniformity
        chrom_loc_depth_per_pool.append([fields.chrom_loc,fields.depth,fields.pool])  # a list of all locations with depth and pool number --- used to calculate uniformity per pool
    
    #---------UNIFORMITY CALCULATION----------
    
    for key,value in chrom_loc_depth_dict.iteritems():
        read_depth_sum += value
       
    mean_read_depth = read_depth_sum/len(chrom_loc_depth_dict) # mean read depth = total depth / total bp
       
    for chrom_loc in chrom_loc_depth_dict:
        if chrom_loc_depth_dict[chrom_loc] >= (0.2 * mean_read_depth): # calculate number of bases that meet uniformity threshold (0.2x mean read depth)
            uniformity_counter += 1                            # count bases      
        else:
            pass
           
    uniformity = (float(uniformity_counter)/float(len(chrom_loc_depth_dict)))*100   # calculate uniformity = bases at 0.2x mean read depth / total evaluated bases

    text = "<b>Mean read depth: %s<br>" % mean_read_depth
    text += "Uniformity: %s%%" % (int(uniformity))
    text += "</b><br><hr>"
    text += """      <table class="GeneratedTable">
                    <thead>
                    <tr>
                    <th>Pool</th>
                    <th>Details</th>
                    </tr>
                    </thead>
                    <tbody>"""
    

    
    #----------POOL-SPECIFIC UNIFORMITY CALCULATION-----------
    
    for pool in set(pools):
        read_depth_sum_per_pool = 0
        base_number_sum = 0
        pool_read_depth_dict[pool] = [read_depth_sum_per_pool,base_number_sum]
    
    for chrom_loc in chrom_loc_depth_per_pool:
        pool = chrom_loc[2]
        depth = chrom_loc[1]
        pool_read_depth_dict[pool][0] += depth
        pool_read_depth_dict[pool][1] += 1
    
    
    
    for key,value in pool_read_depth_dict.iteritems():
        uniformity_counter = 0
        for chrom_loc in chrom_loc_depth_per_pool:
            pool_number = key
            read_depth_sum = value[0]
            base_number_sum = value[1]
            depth_at_pos = chrom_loc[1]
            pool = chrom_loc[2]
            if pool_number == pool and depth_at_pos >= (0.2*(read_depth_sum/base_number_sum)):
                uniformity_counter += 1
        uniformity_per_pool = (float(uniformity_counter)/float(base_number_sum))*100
        text += """<tr><th rowspan="2"><center>%s</center></th>""" % pool_number
        text += """<td>Uniformity: %s%%</td></tr>""" % (int(uniformity_per_pool))
        text += """<tr><td>Mean read depth: %s</td></tr>""" % (int(read_depth_sum/base_number_sum))
    text += """</tbody>
                </table>"""
    return text

def autodetermine_target_regions():
    base_url = 'http://10.80.157.20'
    #base_url = 'http://localhost'
    resp = requests.get('%s/rundb/api/v1/experiment?format=json&order_by=date&expName__iexact=S15-26300' % base_url, auth=('downstream','downstream')).json()
    plan_url = resp['objects'][0]['plan']
    #print json.dumps(resp, sort_keys=True, indent=4, separators=(',', ': '))
    plan_resp = requests.get('%s/%s?format=json' % (base_url,plan_url), auth=('downstream','downstream')).json()
    bedfile = plan_resp['bedfile']
    if re.search("OCP",bedfile):
        bed = "OCP"
    elif re.search("CCP",bedfile):
        bed = "CCP"
    else:
        bed = None
        
    return bed
    #print json.dumps(resp, sort_keys=True, indent=4, separators=(',', ': '))

## Initialize variables

bam_dir = sys.argv[1]
results_dir = sys.argv[2]
analysis_name = sys.argv[3]
print analysis_name

if len(sys.argv) == 4:
    regions_selection = autodetermine_target_regions()
    bed = define_regions(regions_selection,"Female")
else:
    try:
        regions_selection = sys.argv[4]
        gender = sys.argv[5]
        bed = define_regions(regions_selection,gender)
    except:
        print "ERROR: Unable to parse arguments successfully"
        sys.exit(1)


## MAIN

output = open_output(results_dir)
write_html_header(output)
bams = select_BAMs_for_analysis(bam_dir)
for bam in bams:
    basename = determine_base_filename(bam)
    bedtools_coverage(bam_dir, bam, bed, results_dir, basename)
    text = uniformity_plugin(results_dir, basename)
    write_html_per_bam(output, results_dir, basename, text)
