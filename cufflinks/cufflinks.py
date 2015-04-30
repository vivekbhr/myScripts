__version__ = "Cufflinks v0.1.1"

__description__ = """
    {version}
    Abdullah H. Sahyoun - April 27, 2015
    email: sahyoun@ie-freiburg.mpg.de
    ---------------------------------
    Cufflinks pipeline for processing RNA sequencing data from high throughput sequencing.

    This software is distributed WITHOUT ANY WARRANTY!

    Example:
        python cufflinks.py -i path/to/main/bam/input/directory -o path/to/output/dir -t 10 -g dm3 -sampleinfo sampleinfo.txt


    The organisation of the bam files in the input directory should be as follow:

    input
        ---sample1
            ---- sample1.bam
            ---- sample1.bam.bai
        ---sample2
            ---- sample2.bam
            ---- sample2.bam.bai
        ---sample3
            ---- sample3.bam
            ---- sample3.bam.bai


    The sampleInfo file should be something like this:

        nb_replicates	2
        input	GFP_kd_total_unique
        sample1	Msl3_kd_total_unique
        sample2	Msl3_Rrp6_kd_total_unique
        ...
        ...

    """.format(version=__version__)



from os import listdir
from os.path import isfile, join
import argparse
import os
import sys
import re
from collections import OrderedDict
import time
import datetime
import shutil
import subprocess
import textwrap

'''


TODO:
-----

Cuffcompare function needs to be implemented to do the comparison and see if we found any novel transcripts


As mentioned in the cufflinks paper:
-------------------------------------
find . -name transcripts.gtf > gtf_out_list.txt ## get all pathes of transcripts.gt

cuffcompare -i gtf_out_list.txt -r ref_genes.gtf

'''


script_path = os.path.dirname(os.path.realpath(__file__)) + "/"


## Path defaults to all needed scripts and programs
cufflinks_path = '/package/cufflinks-2.2.1.Linux_x86_64/'



### Default variables

samtools_mem = 2
samtools_threads = 2

def parseArguments(args=None):

    requiredArgs = getRequiredArgs()
    parser = argparse.ArgumentParser(parents=[requiredArgs],prog='Cufflinks.py',
    formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(__description__),conflict_handler='resolve')

    args = parser.parse_args(args)


    ## Get reference data paths from config file
    ref_cfg_file_path = os.path.join(script_path, "config_files/{}.cfg".format(args.genome))

    if not os.path.isfile(ref_cfg_file_path):
        print "Error! Configuration file NOT found for {}: {}".format(args.genome, ref_cfg_file_path)
        exit(1)
    configs = parse_config(ref_cfg_file_path)
    try:
        args.fasta_index = configs["fasta_index"]
        args.genome_index = configs["genome_index"]
        args.transcriptome_index = configs["transcriptome_index"]
        args.gtf = configs["gtf"]
        args.bed = configs["bed"]
    except:
        print "Error! Unable to read paths from config file:", ref_cfg_file_path
        exit(1)


    return(args)


def getRequiredArgs():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('Required arguments')

    # define the arguments

    required.add_argument('-i',dest="indir", required=True, help='Path to bam files')
    required.add_argument("-o", "--output-dir", dest="outdir", required=True, help="Output directory")
    required.add_argument("-g", "--genome", dest="genome", required=True, help="Reference genome build")
    required.add_argument('-sampleinfo',dest="sampleinfo",required=True, help='Path to the config file. This file contains the mapping for each sample and replicates i.e., which is the input and which is the treatment for each replicate.')
    parser.add_argument('-t',dest="threads",required=True,help='Number of threads (default: 1)', default=1)
    parser.add_argument("--overwrite", dest="overwrite", action = "store_true", default=False, help="Overwrite results in existing folders!")
    parser.add_argument("--norm", dest="norm", action = "store_true", help="cuffdiff library normalisation method", default="geometric")
    parser.add_argument('-libtype',dest='libtype',help='Library type used in sequencing', default="fr-unstranded",choices=['ff-firststrand','ff-secondstrand','ff-unstranded','fr-firststrand','fr-secondstrand','fr-unstranded','transfrags'])

    return parser


def parse_config(file_path):
    """
    Parse a configuration text file, e.g. *.cfg
    """
    options = {}
    for line in file(file_path):
        line = line.rstrip("\r\n")
        if not line.startswith("#"):
            if "=" in line:
                key, value = line.split('=', 1)
                value = value.strip()
                value = value.strip('"')
                value = value.strip("'")
                options[key] = value
    return options





def run_Cufflinks(args):

    """
    Run Cufflinks! with user specified options.
    :rtype : object
    """

    infiles = []
    analysis_name="Cufflinks"
    args.analysis_counter += 1
    cufflinks_outdir = os.path.join(args.outdir,analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)


    #print "Outdir:", outdir
    if args.overwrite and os.path.isdir(cufflinks_outdir):
        shutil.rmtree(cufflinks_outdir)

    if os.path.isdir(cufflinks_outdir):
        print "Output folder already present: {}".format(cufflinks_outdir)
    else:
        os.mkdir(cufflinks_outdir)
        os.chdir(cufflinks_outdir)
        cwd = os.getcwd()

    aligner="other"

    ### get all bam files within the main directory
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(args.indir)):
        for filename in [f for f in filenames if f.endswith(".bam")]:
            infiles.append(os.path.join(dirpath, filename))

    if any("accepted_hits.bam" in s for s in infiles):
        aligner = "tophat"


    print "Aligner identified: ",aligner
    for infile in infiles:
        if aligner=="tophat":
            if "accepted_hits" in infile:

                basename = infile.split("/")[-2]        ### get the sample name to use as output folder name for each sample
                sample_outpath = os.path.join(cufflinks_outdir,basename)

                os.system("mkdir {}".format(sample_outpath))
                print
                print "{}cufflinks -p {} -g {} -b {} -u --library-type {} -o {} {}".format(cufflinks_path,args.threads,args.gtf,args.fasta_index,args.libtype,sample_outpath,infile)
                print


                os.system("{}cufflinks -p {} -g {} -b {} -u --library-type {} -o {} {}".format(cufflinks_path,args.threads,args.gtf,args.fasta_index,args.libtype,sample_outpath,infile))
                os.system("echo '"+sample_outpath+"/transcripts.gtf' >> "+args.outdir+"assemblies.txt")


            ### other mode
            #os.system("/package/cufflinks-2.2.1.Linux_x86_64/cufflinks -p "+threads+" -G "+GTF_file+" -b "+ ref_fasta+" -u --library-type "+library_type+" -o "+full_outpath+" "+ full_bam_path)		### reference based assembly (RABT)
            #print "/package/cufflinks-2.2.1.Linux_x86_64/cufflinks -p " +threads+" -G "+GTF_file+" -b "+ ref_fasta+" -u --library-type "+library_type+" -o "+full_outpath+" "+ full_bam_path
            #os.system("/package/cufflinks-2.2.1.Linux_x86_64/cufflinks -p "+threads+ " --library-type "+library_type+" -g "+GTF_file+" -b "+ ref_fasta+" -u -o "+full_outpath+" "+ full_bam_path)		### Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.

        else:

            basename = infile.split("/")[-1].split(".")[0]
            sample_outpath = os.path.join(cufflinks_outdir,basename)

            os.system("mkdir {}".format(sample_outpath))

            print
            print "{}cufflinks -p {} -g {} -b {} -u --library-type {} -o {} {}".format(cufflinks_path,args.threads,args.gtf,args.fasta_index,args.libtype,sample_outpath,infile)
            print

            os.system("{}cufflinks -p {} -g {} -b {} -u --library-type {} -o {} {}".format(cufflinks_path,args.threads,args.gtf,args.fasta_index,args.libtype,sample_outpath,infile))
            os.system("echo '"+sample_outpath+"/transcripts.gtf' >> "+args.outdir+"assemblies.txt")

            #### other mode to run cufflinks
            #os.system("/package/cufflinks-2.2.1.Linux_x86_64/cufflinks -p "+threads+ " --library-type "+library_type+" -g "+GTF_file+" -b "+ ref_fasta+" -u -o "+full_outpath+" "+ full_bam_path)		### Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.







def run_Cuffmerge(args):

    '''
       @Input:
            bam_output_folder: the folder that contains the sorted/indexed bam files
            GTF_file: reference transcriptome annotation file in GTF format
            ref_fasta : reference genome fasta file
            output_folder: folder location where the merged trascripts should be merged
    '''

    analysis_name="Cuffmerge"
    args.analysis_counter += 1
    cuffmerge_outdir = os.path.join(args.outdir,analysis_name)
    args.cuffmerge_outdir = cuffmerge_outdir
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)


    #print "Outdir:", outdir
    if args.overwrite and os.path.isdir(cuffmerge_outdir):
        shutil.rmtree(cuffmerge_outdir)

    if os.path.isdir(cuffmerge_outdir ):
        print "Output folder already present: {}".format(cuffmerge_outdir )
    else:
        os.mkdir(cuffmerge_outdir )
        os.chdir(cuffmerge_outdir )
        cwd = os.getcwd()


    assemblies_path = os.path.join(args.outdir,"assemblies.txt")

    print
    print "{}cuffmerge -p {} -g {} -s {} -o {} {}".format(cufflinks_path,args.threads,args.gtf,args.fasta_index,cuffmerge_outdir,assemblies_path)
    #os.system("{}cuffmerge -p {} -g {} -s {} -o {} {}".format(cufflinks_path,args.threads,args.gtf,args.fasta_index,cuffmerge_outdir,assemblies_path))
    print









def run_Cuffdiff(args):
    '''
        @Input:
          An example of a config file should be something like that:

            nb_replicates	2
            input	GFP_kd_total_unique
            sample1	Msl3_kd_total_unique
            sample2	Msl3_Rrp6_kd_total_unique


        The replicates of each sample (i.e., sam and final bam files) should be named in the following way, e.g., Rep1_GFP_kd_total_unique, Rep2_GFP_kd_total_unique, ..., RepN_GFP_kd_total_unique


        The function parses the config file and extract the needed sample names and save it in a dict.
        Then it loops through all the bamfiles, and run Cufflinks, then Cuffmerge.

        After merging the transcripts, the scripts run cuffdiff for the input against all other samples.



        example:
            python cufflinks.py all -config_file cufflinks_config.txt -bam_files /Users/abdullah/Documents/work/Giuseppe/20141204_polyA_RNA_HiSeq/20141219_fraction_RNAseq_analysis//unique_bams/total/ -gtf_file ../../genomes/GTF/dm3_transcriptome_annotation.gtf -ref_fasta ../../genomes/fasta/dm3.fa -threads 4 -o ~/Downloads/

    '''


    infiles = []
    analysis_name="Cuffdiff"
    args.analysis_counter += 1
    cuffdiff_outdir = os.path.join(args.outdir,analysis_name)
    print "\n{} {}) {}".format(datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]'), args.analysis_counter, analysis_name)


    #print "Outdir:", outdir
    if args.overwrite and os.path.isdir(cuffdiff_outdir):
        shutil.rmtree(cuffdiff_outdir)

    if os.path.isdir(cuffdiff_outdir):
        print "Output folder already present: {}".format(cuffdiff_outdir)
    else:
        os.mkdir(cuffdiff_outdir)
        os.chdir(cuffdiff_outdir)
        cwd = os.getcwd()


    configs = {}		### contains the input of the config file
    replicates_dict = {} 	### contains all the data of the replicates, for each sample in the config file, contains a list of all path to replicates files

    ### parse the config file and get informations about the samples
    sample_config = open(args.sampleinfo, "r")
    for l in sample_config:
        if not l.startswith("#"):
            l = l.strip().split()
            if not configs.has_key(l[0]):
                configs[l[0]] = l[1]


    for path, subdirs, files in os.walk(args.indir):
        for name in files:
            full_bam_path = os.path.join(path, name)
            sample_name_replicate_number_removed = re.sub('Rep\d_', '', name).split(".")[0]			### remove the RepN_ and remove the .bam

            #print full_bam_path
            if name.endswith(".bam"):
                basename = name.split(".")[0]

                for sample in configs:

                    if configs[sample] == sample_name_replicate_number_removed:
                        if not replicates_dict.has_key(sample):
                            replicates_dict[sample] = []
                        replicates_dict[sample].append((full_bam_path))			### contains all data about samples and their replicates in a non-sorted way


    ### the sorting step is important to maintain that replicates are put in the same order in different sample in the cuffdiff command line and not Rep1 with Rep3 or so.
    for k in replicates_dict:
        replicates_dict[k] = sorted(replicates_dict[k])



    print
    for sample in replicates_dict:
        if sample=='input':

            continue

        #os.system("mkdir {}{}".format(cuffdiff_outdir,args.norm))
        #os.system("mkdir {}{}{}".format(cuffdiff_outdir,args.norm,configs['input']+"_vs_"+configs[sample]))

        print "mkdir {}.{}".format(cuffdiff_outdir,args.norm)
        print "mkdir {}.{}/{}".format(cuffdiff_outdir,args.norm,configs['input']+"_vs_"+configs[sample])
        sample_output_path = "{}.{}/{}".format(cuffdiff_outdir,args.norm,configs['input']+"_vs_"+configs[sample])

        #os.system("mkdir "+output_folder+"cuffdiff."+normalization_method+"/")
        #os.system("mkdir "+output_folder+"cuffdiff."+normalization_method+"/"+configs['input']+"_vs_"+configs[sample])


        print
        print configs['input']+"_vs_"+configs[sample]
        print '---------------------------------------------------'
        print

        cuffdiff_cmd = ""
        cuffdiff_cmd = "/package/cufflinks-2.2.1.Linux_x86_64/cuffdiff -u -p {} --library-norm-method {} -b {} -o {} {}/merged.gtf ".format(args.threads,args.norm,args.fasta_index,sample_output_path,args.cuffmerge_outdir)


        for i in range(len(replicates_dict['input'])):          ### the input replicates bam to the cuffmerge command separated by a comma

            cuffdiff_cmd+=replicates_dict['input'][i]+","


        cuffdiff_cmd = cuffdiff_cmd[:-1]+" "        ### remove the additional comma and replace it by a space (to separate the bam files of input and treatment)


        for i in range(len(replicates_dict[sample])):

            cuffdiff_cmd+=replicates_dict[sample][i]+","


        cuffdiff_cmd = cuffdiff_cmd[:-1]
        print
        print cuffdiff_cmd.strip(",")
        print
        #os.system(cuffdiff_cmd.strip(","))


    print















if __name__ == "__main__":
    args = parseArguments()
    args_dict=vars(args)

    args.analysis_counter = 0   # Counter for analyzes conducted

    #run_Cufflinks(args)
    run_Cuffmerge(args)
    run_Cuffdiff(args)

    print "Done!"





