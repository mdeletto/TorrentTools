#/usr/bin/python

import re,sys

class BED_fields():
    def __init__(self,fields):
        self.entry = ""
        for field in fields:
            self.entry += field
        self.chrom = fields[0]
        self.start_pos = fields[1]
        self.end_pos = fields[2]
        self.amplicon_id = fields[3]
        pool = re.search("Pool=(.+?);",self.entry)
        self.pool = pool.group(1)

def write_temp_bed(bed_fields):

    def check_multiple_pools(bed_fields):
        bed_fields.pool = bed_fields.pool.split(",")
    
    check_multiple_pools(bed_fields)
    
    for pool in bed_fields.pool:
        print "%s\t%s\t%s\t%s\t%s" % (bed_fields.chrom,bed_fields.start_pos,bed_fields.end_pos,bed_fields.amplicon_id,pool)

if len(sys.argv) != 2:
    print "ERROR: Script requires input AmpliSeq BED file with primer pool information!!!"
    sys.exit(1)
else:
    if not sys.argv[1].endswith(".bed"):
        print "ERROR: Script requires input AmpliSeq BED file with primer pool information!!!"
        sys.exit(1)
    else:
        input_file = sys.argv[1]
 
#with open("/home/michael/Desktop/IonReporter/AmpliSeq_Exome/AmpliSeqExome.20141113.designed.bed","r") as f:
with open(input_file,"r") as f:
    for line in f.readlines():
        if not re.match("^chr",line):
            print line.strip()
        else:
            line = line.strip()
            fields = line.split("\t")
            #print fields
            bed_fields = BED_fields(fields)
            write_temp_bed(bed_fields)
            #print bed_fields.pool