#!/usr/bin/env python
__description__="""
SGEeasy v0.4.1 - A launcher executing a single line of BASH commands on the grid engine
Fabian Kilpert <- March 11, 2015
Vivek Bhardwaj <- June 2019

SGEeasy is a wrapper for easy execution of user shell commands on the
grid engine. Technically, a temporary job script (*.SGE) is generated
that is submitted to qsub. 
SGE has been designed to be used on MPI-IE deep servers!
Examples:
  SGE 'date'                                      submit shell command i.e. 'date'
  SGE 'sh a.sh'                                   submit a script
  SGE --qsub '-pe smp 2' 'sh a.sh'                2 threads per job
  SGE --qsub '-pe smp 2 -t 1:3' 'sh a.sh'         2 threads per job, 3 array jobs
                                                    (available: $SGE_TASK_ID)
  SGE --qsub '-j y' 'sh a.sh'                     merge STDOUT and STDERR
  SGE --qsub '-q deep.q' 'sh a.sh'                submit to a specific queue 
  SGE --qsub '-M my@email.com -m ea' 'sh a.sh'    email notification when job ends or aborts
"""
import argparse
import datetime
import os
import subprocess
import sys
import textwrap
## User defaults
QSUB_OPTS = ""
SGE_LOG_DIR = ""
EMAIL = ""
## other defaults passed to the SGE
SGE_SCRIPT_HEADER = """
#!/bin/bash
TMPDIR=/hpc/hub_oudenaarden/vbhardwaj/tmpdir

scratch=$(mktemp -p $TMPDIR -d -t tmp.XXXXXXXXXX)
export TMPDIR=$scratch
export TMP=$scratch
export TEMP=$scratch
function cleanup {
    rm -rf $scratch
}
trap cleanup SIGUSR2
cd $scratch
""".rstrip()
SGE_SCRIPT_FOOTER = """
cleanup 
""".rstrip()
## Internal default variables (Do not change!!!)
this_script = os.path.basename(sys.argv[0])
cwd = os.getcwd()
now = datetime.datetime.now()
## Functions
def headline(h, max_len=64):
    """
    Output headline on stdout.
    """
    return " ".join( [(max_len-len(h))/2*"-", h, (max_len-len(h))/2*"-"] )
def parse_args():
    """
    Parse arguments from the command line.
    """
    
    parser = argparse.ArgumentParser(
    prog="SGE",
    formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(__description__))
    ## optional
    parser.add_argument("--qsub", dest="qsub_opts", metavar="STR", help="qsub option string (default: '{}')".format(QSUB_OPTS), type=str, default=QSUB_OPTS)
    parser.add_argument("-l", "--logdir", dest="log_dir", metavar="LOGDIR", help="log directory for grid engine output (default: '{}')".format(SGE_LOG_DIR), type=str, default='')
    parser.add_argument("-n", "--name", dest="name", metavar="JOBNAME", help="job name (default: your user name i.e.'{}')".format(os.getlogin()), type=str, default=os.getlogin())
    parser.add_argument("-k", "--keepscript", dest="keep_script", action="store_true", default=False, help="keep auto-generated job script (*.SGE)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False, help="verbose output")
    ## positional
    parser.add_argument('cmd', metavar="'BASH COMMANDS'", action="store", help="BASH commands to be executed by the grid engine")
    
    args = parser.parse_args()
    
    ## Variable sanity checking
    
    ## Log dir
    if args.log_dir:
        args.log_dir = os.path.join(os.getcwd(), args.log_dir)
    else:
        if SGE_LOG_DIR:
            args.log_dir = SGE_LOG_DIR
        else:
            args.log_dir = os.getcwd()
    if not os.path.isdir(args.log_dir):
        os.makedirs(args.log_dir)
    
    return args
def main():
    args = parse_args()
    
    outfile = "{}.SGE".format(args.name)
    
    ## Generating the shell script to be passed with qsub
    script = []
     
    ## Header
    script.extend( [line for line in SGE_SCRIPT_HEADER.splitlines() if line != ''] )
    
    # if "-N" not in args.qsub_opts.split():
        # script.insert(1,"#$ -N {}".format(args.name))
    
    if "-t" in args.qsub_opts.split():
        if not "-o" in args.qsub_opts.split():
            script.insert(1,"#$ -o {}/{}.{}.{}.{}.{}.o".format(args.log_dir, now.strftime('%Y%m%d_%H%M%S'), "$JOB_ID", "$TASK_ID", this_script, args.name))
        if not "-e" in args.qsub_opts.split():
            script.insert(2,"#$ -e {}/{}.{}.{}.{}.{}.e".format(args.log_dir, now.strftime('%Y%m%d_%H%M%S'), "$JOB_ID", "$TASK_ID", this_script, args.name))
    else:
        if not "-o" in args.qsub_opts.split():
            script.insert(1,"#$ -o {}/{}.{}.{}.{}.o".format(args.log_dir, now.strftime('%Y%m%d_%H%M%S'), "$JOB_ID", this_script, args.name))
        if not "-e" in args.qsub_opts.split():
            script.insert(2,"#$ -e {}/{}.{}.{}.{}.e".format(args.log_dir, now.strftime('%Y%m%d_%H%M%S'), "$JOB_ID", this_script, args.name))
    
    ## Command line
    script.append(args.cmd)
    ## Footer
    script.extend( [line for line in SGE_SCRIPT_FOOTER.splitlines() if line != ''] )
    
    ## Output to stdout
    if args.verbose:
        print headline(outfile)
        for line in script:
            print line
    
    ## to file
    with open(outfile, "w") as f:
        for line in script:
            f.write(line+"\n")
    
    ## Run qsub
    if args.verbose:
        print headline("qsub")
    qsub = ["qsub", "-notify"]
    qsub += args.qsub_opts.split()
    qsub.append(outfile)
    
    if args.verbose:
        print " ".join(qsub), "\n"
    
    ## Send to qsub
    subprocess.call(qsub)
    ## remove job script
    if not args.keep_script:
        os.remove(outfile)
        if args.verbose:
            print "\nAuto clean-up:", outfile
        
if __name__ == "__main__":
    main()

