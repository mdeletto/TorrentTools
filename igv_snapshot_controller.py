#!/usr/bin/python2.7

import socket
import sys
import time
import argparse
import os
global args

def output_directory_handling(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)



def set_snapshot_dir(socket, dir):
    """Set IGV snapshot directory for this session."""
    socket.send("snapshotDirectory %s" % dir)
    socket.send("\n")
    time.sleep(2)

def go_to_locus_and_snapshot(socket, loci):
    """Go to a locus and take a snapshot."""
    if loci is None or loci=="None":
        # if locus is not specified, default to whole-genome view
        print "WARNING: No loci were supplied.  Defaulting to whole-genome view."
        socket.send("goto all")
        socket.send("\n")
        time.sleep(2)
        socket.send("snapshot")
        socket.send("\n")
        time.sleep(2)      
    else:
        for locus in loci:
            socket.send("goto %s" % str(locus))
            socket.send("\n")
            time.sleep(2)
            socket.send("snapshot %s" % (str(locus)+".png"))
            socket.send("\n")
            time.sleep(2)

def autodetect_input_file_format(socket, files, force_loading):
    compatible_extensions = {'bam' : 'BAM format',
                             'bed' : 'BED format',
                             'bedgraph' : 'BedGraph format',
                             'bb' : 'bigBed format',
                             'bigWig' : 'bigWig format',
                             'bw' : 'bigWig format',
                             'broadPeak' : 'broadPeak format',
                             'cbs' : 'CBS format',
                             'seg' : 'CBS format',
                             'cn' : 'CN format',
                             'gct' : 'GCT format',
                             'gff' : 'GFF format',
                             'gff3' : 'GFF3 format',
                             'gtf' : 'GTF format',
                             'gistic' : 'GISTIC format',
                             'linear' : 'GWAS format',
                             'logistic' : 'GWAS format',
                             'assoc' : 'GWAS format',
                             'qassoc' : 'GWAS format',
                             'gwas' : 'GWAS format',
                             'igv' : 'IGV format',
                             'loh' : 'LOH format',
                             'maf' : 'MAF (Multiple Alignment Format or Mutation Annotation Format) format',
                             'mut' : 'MUT format',
                             'narrowPeak' : 'narrowPeak format',
                             'psl' : 'PSL format',
                             'res' : 'RES format',
                             'sam' : 'SAM format',
                             'snp' : 'SNP format',
                             'tab' : 'Tab-delimited format',
                             'tdf' : 'TDF format',
                             'vcf' : 'VCF format',
                             'wig' : 'WIG format'
                            }

    for input_file in files:
        file_extension = input_file.split('.')[-1]
        if file_extension in compatible_extensions.keys():
            print "INPUT FILE ::: OK: AUTODETECTED FORMAT: %s as %s" % (input_file, compatible_extensions[file_extension])
            load_file_in_igv(socket, input_file)
        else:
            if force_loading is True:
                print "INPUT FILE ::: WARNING: Format not autodetected.  Attempting to load file %s.  This may cause errors in IGV if the file cannot be loaded.  You have been warned!" % input_file
                load_file_in_igv(socket, input_file)
            else:
                print "INPUT FILE ::: WARNING: Format not autodetected.  %s will not be loaded.  If you would like to force this file to load, invoke --force_loading." % input_file

def load_file_in_igv(socket, input_file):
    """Given a file and socket information, load a file into IGV."""
    socket.send("load %s" % input_file)
    socket.send("\n")
    time.sleep(2)

def manage_tracks(socket, track_mode):
    """Expand, collapse, or squish tracks."""
    socket.send("%s" % track_mode)
    socket.send("\n")
    time.sleep(2)

def parse_input_files(socket, files, force_loading):
    """Load files in IGV."""
    if files is None or files=="None":
        print "WARNING: No input files were provided.  No tracks will be loaded"
    else:
        autodetect_input_file_format(socket, files, force_loading)

def main():

    global socket
    
    desc="""When given a list of files and coordinates, this script will attempt to load all files in a local IGV instance and take snapshots at select loci. Requires IGV already be open and able to listen on port 60151."""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    general_arguments = parser.add_argument_group('GENERAL OPTIONS')
    general_arguments.add_argument("--files","-f", help="Files to load into IGV as tracks.  Please see IGV's website for compatible filetypes",dest='input_files',action='append',nargs='+', default=None)
    general_arguments.add_argument("--loci","-l", help="Files to load into IGV as tracks.  Please see IGV's website for compatible filetypes",dest='loci',action='store',nargs='+',default=None)
    general_arguments.add_argument("--output_dir","-d", help="Full path to output directory.  The script will create this directory if it does not exist.",dest='output_dir',action='store',default=os.getcwd()+"/igv")
    general_arguments.add_argument("--port","-p", help="Port that IGV communicates over",dest='port',action='store',nargs=1,default=60151)
    general_arguments.add_argument("--ip_address","-a", help="IP address of IGV instance",dest='ip_address',action='store',nargs=1, default="127.0.0.1")
    general_arguments.add_argument("--force_loading", help="Force loading of all input files (i.e. load even if file extension check fails)",dest='force_loading',action='store_true',default=False)
    general_arguments.add_argument("--track_mode", help="Choose how IGV will display tracks.",dest='track_mode',action='store',choices=['collapse','expand','squish'], default='expand')
    
    args = parser.parse_args()
    
    if args.input_files is not None:
        args.input_files = [item for sublist in args.input_files for item in sublist]

    TCP_IP = args.ip_address
    TCP_PORT = args.port
    
    ### CONNECT TO IGV PORT ###
    print "ADDRESS:", TCP_IP
    print "PORT:", TCP_PORT
    # Connect to socket
    socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    socket.connect( (TCP_IP, TCP_PORT) )
    
    ### INPUT FILE PROCESSING ###
    print "INPUT FILES:", args.input_files
    # load files
    parse_input_files(socket, args.input_files, args.force_loading)
    manage_tracks(socket, args.track_mode)
    
    ### IGV SNAPSHOT DIRECTORY HANDLING ###
    print "IGV SNAPSHOT OUTPUT DIR:", args.output_dir              
    # Change cwd to output directory.  Directory is created if it doesn't exist
    output_directory_handling(args.output_dir)
    # Set snapshot directory.
    set_snapshot_dir(socket, args.output_dir)
    
    ### LOAD LOCI ###
    print "LOCI TO LOAD:", args.loci
    # Loop over loci and take a snapshot
    go_to_locus_and_snapshot(socket, args.loci)

if __name__ == '__main__':
    main()




    
