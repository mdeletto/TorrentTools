#

import subprocess,sys,json,optparse,urllib,re

global opts

desc="""Downloads a VCF file for a given IR analysis.  Compatible with IR v4.0 through v4.4.  Later IR versions are untested."""
usage="python IR_download --analysis_name <analysis_name>"
parser = optparse.OptionParser(description=desc,usage=usage)

group = optparse.OptionGroup(parser, "General options")

group.add_option("--authorization","-a", help="Authorization API key for IonReporter",dest='auth_key',action='store')
group.add_option("--ionreporter_url","-i", help="IP address or hostname of IR server (ex: 10.80.157.179)",dest='url',action='store')
group.add_option("--keep_tsv", help="Keep .tsv file that is download in the IR.zip.  Default is to remove it.",dest='keep_tsv',action='store_true',default=False)
group.add_option("--keep_entire_zip", help="Keep the entire IR.zip file.  Default is to remove it.",dest='keep_zip',action='store_true',default=False)
group.add_option("-v", help="verbose <DEFAULT=[%default]>", dest='verbose',action='store_true',default=False)
parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Single mode",
                    "Downloads VCF for a single analysis")

group.add_option('--analysis_name', help="IR analysis name", dest='input_analysis_name', action='store')
parser.add_option_group(group)

group = optparse.OptionGroup(parser, "Batch mode",
                    "Downloads VCFs from multiple analyses.")
 
group.add_option("-b", help="Path to batch file.  Expects tabular file with 2 columns [<analysis_name>tab<analysis_id_if_known>].  If you are adding column titles, please add # to the start of those lines.", dest='batch_input',action='store')
parser.add_option_group(group)

(opts, args) = parser.parse_args()

def IR_locate_variant_zip(url_encoded_basename,ionreporter_id,IR_url,authorization_key):    
    try:
        if ionreporter_id is None:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&name=%s" 2> /dev/null""" % (authorization_key,IR_url,url_encoded_basename)],shell=True,stdout=subprocess.PIPE)
        else:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "https://%s/webservices_42/rest/api/analysis?format=json&name=%s&id=%s" 2> /dev/null""" % (authorization_key,IR_url,url_encoded_basename,ionreporter_id)],shell=True,stdout=subprocess.PIPE)
        output, err = proc.communicate()
    except:
        print "ERROR: Unable to communicate with server.  Check Authorization key, server address, and your network connectivity."
        sys.exit(1)
    try:
        try:
            data = json.loads(output)
            if opts.verbose is True:
                print json.dumps(data,indent=4)
        except:
            print "Could not load json string"
        unfiltered_variants = data[0]['data_links']['unfiltered_variants']
        filtered_variants = data[0]['data_links']['filtered_variants']
        return unfiltered_variants
    except Exception,e:
        #print str(e)
        print "ERROR: Unable to process IonReporter links.  Is this a valid analysis name?  Aborting..."
        return 0

def IR_download_variant_zip(basename,variant_link,authorization_key):
    if bool(variant_link) is False:
        print "ERROR: Unable to process case %s" % basename
        sys.exit(1)
    else:
        try:
            proc = subprocess.Popen(["""curl -k -H "Authorization:%s" "%s" 2> /dev/null -o IR.zip; unzip IR.zip && unzip %s.zip; cp ./Variants/*/*.vcf %s.vcf && cp ./Variants/*/*.tsv %s.tsv; rm -rf %s.zip QC Variants Workflow_Settings""" % (authorization_key,variant_link,basename,basename,basename,basename)],shell=True,stdout=subprocess.PIPE)
            output, err = proc.communicate()
        except:
            print "Unable to download and/or unzip IonReporter files.  Aborting..."
            sys.exit(1)

def command_line_parsing(opts):
    
    mandatory_options = ['auth_key','url']
    for m in mandatory_options:
        # Making sure all mandatory options appeared
        if not opts.__dict__[m]:
            print "Mandatory option is missing!\n"
            parser.print_help()
            sys.exit()       

    if opts.input_analysis_name is not None:
        print "------------------------------------"
        print "OK: Proceeding in single sample mode"
        print "------------------------------------"
        return "SINGLE"
    elif opts.batch_input is not None:
        print "------------------------------------"
        print "OK: Proceeding in batch sample mode"
        print "------------------------------------"
        return "BATCH"
    elif opts.batch_input is not None and opts.input_analysis_name is not None:
        print "ERROR: Cannot run both single sample and batch sample modes simulataneously"    
        parser.print_help()
        sys.exit() 
    elif opts.batch_input is None and opts.input_analysis_name is None:
        print "ERROR: No input analysis name provided"
        parser.print_help()
        sys.exit()
    else:
        print "ERROR: Command line parsing failed"
        parser.print_help()
        sys.exit()        

      
mode = command_line_parsing(opts)
analysis_name_id = {}

if mode == "SINGLE":
    analysis_name_id[opts.input_analysis_name] = None
elif mode == "BATCH":
    with open(opts.batch_input,'r') as f:
        for line in f.readlines():
            line = line.strip()
            if re.match("^#",line):
                pass
            else:
                split_line = line.split("\t")
                analysis_name = split_line[0]
                try:
                    analysis_id = split_line[1]
                except:
                    analysis_id = None
                    print "WARNING: No analysis ID provided w/ %s" % analysis_name
                    print "WARNING: Proceeding with only using analysis name.  This shouldn't be a problem."
                    
                if analysis_id is None:
                    analysis_name_id[analysis_name] = None
                else:
                    analysis_name_id[analysis_name] = analysis_id

for k,v in analysis_name_id.iteritems():
    analysis_name,analysis_id = k,v 
    url_encoded_basename = urllib.quote_plus(analysis_name) # catch whitespace or other problematic characters and encode for url
    modified_basename = analysis_name.replace(" ","_") # IR converts any whitespace characters in <analysis_name> to underscores
    variant_link = IR_locate_variant_zip(url_encoded_basename,None,opts.url,opts.auth_key)
    IR_download_variant_zip(modified_basename, variant_link,opts.auth_key)
    print "SUCCESS: %s (%s) was downloaded sucessfully" % (analysis_name,url_encoded_basename)
