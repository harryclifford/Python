###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#
#   Copyright (C) 2013 Harry Clifford
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#################################################################################
"""
=================================
RNA-Seq Quality Control Pipeline
=================================

Requires the following input files in the working directory:
    1. fastq files - CURRENTLY THIS SCRIPT ONLY TAKES GZIPPED INPUT! (must end with .fastq.gz)
    2. the configuration file - "pipeline.ini"
    (3. a fasta file of contamination sequences to be removed, named "contaminant_remove.fasta"
            if this is not provided, pipeline will make one from fastqc results)
    (4. an alternative to fastqc's standard contaminants file for annotation of overrepresented sequences
            if added, this must be named "contaminant_list.txt"; if not provided, fastqc's standard file is used)


Overview
========

This pipeline will perform quality control and correction for fastq files.

This pipeline performs the following:
    - fastqc analysis before and after
    - trims read ends based on quality
    - removes contaminants
    - filters reads based on average quality

Because each stage produces fastq files for each sample, there will be a build up of files which may take up a large amount of disk space.

Ruffus output (telling you progress of the pipeline and any failings) will be directed to stdout. Outputs of the third-party software will be redirected to summary files.


Usage
=====

The command line used when running this pipeline should consist of: the python command, followed by the name of this file, followed by the stage to run to (if not specified or "all", pipeline runs through completely; if "test", pipeline produces graph of modules to be run and lists tasks to be run in std.out). All other settings can be changed in the pipeline.ini file.

e.g.1. python harryc_QC_pipeline.py all
e.g.2. python harryc_QC_pipeline.py test

You may also specify the number of parallels you wish to use in the command line, which will override the pipeline.ini file. This is simply provided as an integer following the specified stage to run to.

e.g.  python harryc_QC_pipeline.py all 5


Third-party software
====================

This pipeline requires the following third-party software to be installed:
    - fastqc
    - cutadapt
    - fastx_toolkit


Code
====

"""

###################################################################
## For use throughout pipeline
###################################################################

# modules
from ruffus import *
import sys
import os
import Pipeline
import glob
import errno
import bisect
from rpy2.robjects import r as R

# options from ini file
params = Pipeline.getParameters(["pipeline.ini"])

## create directory if does not exist function (bypasses race condition, but still raises any other errors)
def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise



############
## fastqc1
############

@transform(glob.glob('*.fastq.gz'), regex(r'(.+)\.fastq.gz'), r'1_fastqc1/\1_fastqc/fastqc_data.txt')
def fastqc1(infile, outfiles):
    
    makedirs('1_fastqc1')
    
    # builds fastqc commandline
    command = params['fastqc1_fastqc1_bashcommand']
    outdir = "-o 1_fastqc1"
    if 'contaminant_list.txt' in os.listdir(os.getcwd()):
        contams = "-c contaminant_list.txt"
    else: contams = ""
    
    # runs commandline
    os.system(" ".join([command,infile,outdir,contams,"-q"]))
    
    # deletes zip files produced
    zipname = "1_fastqc1/" + infile.split(".fastq.gz")[0] + "_fastqc.zip"
    os.remove(zipname)

## Example of command line produced:
# fastqc test.fastq.gz -o 1_fastqc1 -c contaminants_list.txt



#################
## quality_trim
#################

@follows(fastqc1)
@transform(glob.glob('*.fastq.gz'), regex(r'(.+)\.fastq.gz'), [r'2_quality_trim/outputs/\1.fastq.gz',r'2_quality_trim/summaries/\1.summary'])
def quality_trim(infile, outfiles):
    
    makedirs('2_quality_trim/outputs')
    makedirs('2_quality_trim/summaries')
    
    # builds quality trim commandline
    infiles = params['general_gunzip_bashcommand'] + " -c " + infile + " |"
    command = params['quality_trim_quality_trim_bashcommand']
    options = "-Q 33 -v -z -t " + str(params['quality_trim_quality_threshold']) + "-l " + str(params['quality_trim_minimum_length'])
    outfile1 = "-o " + outfiles[0]
    outfile2 = "> " + outfiles[1]
    
    # runs commandline
    os.system(" ".join([infiles,command,options,outfile1,outfile2]))


## Example of command line produced:
# fastq_quality_trimmer -Q 33 -v -z -t 20 -l 30 -i test.fastq.gz -o 2_quality_trim/outputs/test.fastq.gz > 2_quality_trim/summaries/test.summary




#########################
## contaminants_removal
#########################


## main function - used with differing decoraters below
def cutadapt(infiles,outfiles):
    
    makedirs('3_contaminants_removal/outputs')
    makedirs('3_contaminants_removal/summaries')
    
    # builds contaminants_removal commandline
    command = params['contaminants_removal_cutadapt_bashcommand']
    tmp_contams = open('contaminant_remove.fasta',"r").read().split('\n')
    if '' in tmp_contams: tmp_contams.remove('')
    contams = " ".join( [ " ".join([params['contaminants_removal_remove_from'],x]) for x in tmp_contams ] )
    options = params['contaminants_removal_wildcards'] + " -O " + str(params['contaminants_removal_minimum_overlap']) + " -m " + str(params['contaminants_removal_minimum_length'])
    outfile1 = "-o " + outfiles[0]
    outfile2 = "> " + outfiles[1]
    
    # runs commandline
    os.system(" ".join([command,contams,infiles[0],options,outfile1,outfile2]))

    

## Example of command line produced:
# cutadapt -b AACGAGG -b AAGTCGA -b CGTTAGC --no-match-adapter-wildcards -O 20 -m 30 2_quality_trim/outputs/test.fastq.gz -o 3_contaminants_removal/outputs/test.fastq.gz > 3_contaminants_removal/summaries/test.summary



## secondary function - used with differing decoraters below
@follows(fastqc1)
@merge(fastqc1, 'contaminant_remove.fasta')
def build_contaminant_list(infiles, outfile):
    
    makedirs('RESULTS')
    
    detailed_results = []
    forCutAdapt_results = []
    # loop through folders
    for infile in infiles:
        
        # take data from folder
        stats = open(infile,"r").read().split('\n')
        
        # extract overrepresented sequences
        if '>>Overrepresented sequences\tfail' in stats or '>>Overrepresented sequences\twarn' in stats:
            # if any, extracts all
            if '>>Overrepresented sequences\tfail' in stats: pos1 = stats.index('>>Overrepresented sequences\tfail')
            if '>>Overrepresented sequences\twarn' in stats: pos1 = stats.index('>>Overrepresented sequences\twarn')
            ends = [i for i, x in enumerate(stats) if x == '>>END_MODULE']
            pos2 = ends[bisect.bisect(ends, pos1)]
            
            # adds to results
            tmp_result = [stats[i] for i in range(pos1+2,pos2)]
            tmp_result2 = ['\t'.join([i.split('\t')[0],i.split('\t')[-1].split(' (')[0]]) for i in tmp_result]
            tmp_result3 = [i.split('\t')[0] for i in tmp_result]
            
            detailed_results = detailed_results + tmp_result2
            forCutAdapt_results = forCutAdapt_results + tmp_result3
    
    unique_detailed_results = sorted(list(set(detailed_results)))
    unique_forCutAdapt_results = sorted(list(set(forCutAdapt_results)))
    
    with open("RESULTS/contaminants_removed.txt",'w') as final: final.write('\n'.join(unique_detailed_results))
    with open(outfile,'w') as final: final.write('\n'.join(unique_forCutAdapt_results))



# only does secondary function if file produced does not already exist
if 'contaminant_remove.fasta' not in os.listdir(os.getcwd()):
    @follows(quality_trim)
    @follows(build_contaminant_list)
    @transform(quality_trim, regex(r'2_quality_trim/outputs/(.+)\.fastq.gz'), [r'3_contaminants_removal/outputs/\1.fastq.gz',r'3_contaminants_removal/summaries/\1.summary'])
    def contaminants_removal(infiles,outfiles):
        cutadapt(infiles,outfiles)
else:
    @follows(quality_trim)
    @transform(quality_trim, regex(r'2_quality_trim/outputs/(.+)\.fastq.gz'), [r'3_contaminants_removal/outputs/\1.fastq.gz',r'3_contaminants_removal/summaries/\1.summary'])
    def contaminants_removal(infiles,outfiles):
        cutadapt(infiles,outfiles)



###################
## quality_filter
###################


@follows(contaminants_removal)
@transform(contaminants_removal, regex(r'3_contaminants_removal/outputs/(.+)\.fastq.gz'), [r'4_quality_filter/outputs/\1.fastq.gz',r'4_quality_filter/summaries/\1.summary'])
def quality_filter(infiles, outfiles):
    
    makedirs('4_quality_filter/outputs')
    makedirs('4_quality_filter/summaries')
    
    # builds quality trim commandline
    infile_command = params['general_gunzip_bashcommand'] + " -c " + infiles[0] + " |"
    command = params['quality_filter_quality_filter_bashcommand']
    options = "-Q 33 -v -z -q " + str(params['quality_filter_quality_threshold']) + "-p " + str(params['quality_filter_minimum_percent'])
    outfile1 = "-o " + outfiles[0]
    outfile2 = "> " + outfiles[1]
    
    # runs commandline
    os.system(" ".join([infile_command,command,options,outfile1,outfile2]))


## Example of command line produced:
# fastq_quality_filter -v -q 10 -p 100 -i test.fastq.gz -o 4_quality_filter/outputs/test.fastq.gz > 4_quality_filter/summaries/test.summary




############
## fastqc2
############


@follows(quality_filter)
@transform(quality_filter, regex(r'4_quality_filter/outputs/(.+)\.fastq.gz'), r'5_fastqc2/\1_fastqc/fastqc_data.txt')
def fastqc2(infiles, outfiles):
    
    makedirs('5_fastqc2')
    
    # builds fastqc commandline
    command = params['fastqc2_fastqc2_bashcommand']
    outdir = "-o 5_fastqc2"
    if 'contaminant_list.txt' in os.listdir(os.getcwd()):
        contams = "-c contaminant_list.txt"
    else: contams = ""
    
    # runs commandline
    os.system(" ".join([command,infiles[0],outdir,contams,"-q"]))
    
    # deletes zip files produced
    zipname = "5_fastqc2/" + infiles[0].split(".fastq.gz")[0].rsplit("/",1)[-1] + "_fastqc.zip"
    os.remove(zipname)

## Example of command line produced:
# fastqc 5_fastqc2/test.fastq.gz -o 5_fastqc2 -c contaminants_list.txt




############
## results
############


@follows(quality_filter)
@transform(quality_filter, regex(r'4_quality_filter/outputs/(.+)\.fastq.gz'), r'RESULTS/final_outputs/\1.fastq.gz')
def results1(infiles, outfile):
    
    makedirs('RESULTS/final_outputs')
    # makes (or deletes then makes if already exists) symlink to final output
    if os.path.exists(outfile): os.remove(outfile)
    os.symlink( os.getcwd()+"/"+infiles[0] , outfile )



@follows(fastqc1)
@follows(quality_trim)
@follows(contaminants_removal)
@follows(quality_filter)
@follows(fastqc2)
@follows(results1)
@merge(glob.glob('*.fastq.gz'), 'RESULTS/summary.pdf')
def results2(infile, outfile):
    
    R('''
    
    # loads relevant packages
    suppressMessages(library(gplots))
    suppressMessages(library(png))
    suppressMessages(library(Hmisc))
    
    # acquires filenames
    files <- sort(Sys.glob('5_fastqc2/*_fastqc/fastqc_data.txt'))
    
    # begins pdf (with 5 plots to each page
    pdf('RESULTS/summary.pdf',height=4,width=20)
    par(mfrow=c(1,5))
    
    ## begins loop through files
    for(file in files){
        
        # makes individual filename
        filename <- strsplit( strsplit(file,"/")[[1]][2] ,"_fastqc" )[[1]]
        
        # changes margin settings for upcoming plots
        par(mar=c(0.1,0.1,2,0.1),oma=c(0,0,1,0))
        
        ### plot 1 - summary ###
        
        # obtains fastqc stats
        after_fastqc_stats <- read.delim(file,header=F,stringsAsFactors=F)
        before_fastqc_stats <- read.delim(paste("1_fastqc1/",filename,"_fastqc/fastqc_data.txt",sep=""),header=F,stringsAsFactors=F)
        
        # builds text for plot
        text_tmp <- before_fastqc_stats[grep(">>Basic Statistics",before_fastqc_stats[,1]):min(grep(">>END_MODULE",before_fastqc_stats[,1])),]
        text_readsbefore <- text_tmp[grep("Filename",text_tmp[,1]):(length(text_tmp[,1])-1),]
        
        # plots
        textplot(text_readsbefore,valign="top",cex=1.2, show.rownames=F, show.colnames=F)
        title("\nSummary\n(before processing)")
        box("figure")
        
        
        ### plot 2 - quality_trim ###
        
        # obtains quality trim stats
        trim_stats_tmp <- read.delim(paste("2_quality_trim/summaries/",filename,".summary",sep=""),header=F,stringsAsFactors=F)
        trim_stats <- capitalize(trim_stats_tmp[length(trim_stats_tmp[,1]),])
        
        # obtains plot for before and after
        trim_before <- readPNG(paste("1_fastqc1/",filename,"_fastqc/Images/per_base_quality.png",sep=""))
        trim_after <- readPNG(paste("5_fastqc2/",filename,"_fastqc/Images/per_base_quality.png",sep=""))
        
        # plots
        textplot(paste(trim_stats,"\nBefore and after pipeline:\n",sep=""),valign="top", cex=1.2)
        rasterImage(trim_before, 0,   0.1, 0.5, 0.75)
        rasterImage(trim_after,  0.5, 0.1, 1,   0.75)
        title("\nQuality Trim")
        box("figure")
        
        
        ### plot 3 - contaminants_removal ###
        
        # obtains cutadapt stats
        contams_stats <- read.delim(paste("3_contaminants_removal/summaries/",filename,".summary",sep=""),header=F,stringsAsFactors=F)
        
        # obtains overrepresented sequence data from first fastqc, if exists
        contams_before_test <- before_fastqc_stats[grep(">>Overrepresented sequences",before_fastqc_stats[,1]),2]
        if(contams_before_test == "fail" | contams_before_test == "warn"){
            tmp_subset <- before_fastqc_stats[grep(">>Overrepresented sequences",before_fastqc_stats[,1]):length(before_fastqc_stats[,1]),]
            tmp_subset2 <- tmp_subset[4:(grep(">>END_MODULE",tmp_subset[,1])[1]-1),1]
            contamsbefore <- length(tmp_subset2)/2
        }else{ contamsbefore <- 0 }
        
        # obtains overrepresented sequence data from second fastqc, if exists
        contams_after_test <- after_fastqc_stats[grep(">>Overrepresented sequences",after_fastqc_stats[,1]),2]
        if(contams_after_test == "fail" | contams_after_test == "warn"){
            tmp_subset <- after_fastqc_stats[grep(">>Overrepresented sequences",after_fastqc_stats[,1]):length(after_fastqc_stats[,1]),]
            tmp_subset2 <- tmp_subset[4:(grep(">>END_MODULE",tmp_subset[,1])[1]-1),1]
            contamsafter <- length(tmp_subset2)/2
        }else{ contamsafter <- 0 }
        
        # obtains specific data from cutadapt stats
        trimmed_reads <- unlist(contams_stats)[grep("Trimmed reads:",unlist(contams_stats))]
        too_short_reads <- paste( strsplit(
                                unlist(contams_stats)[grep("Too short reads:",unlist(contams_stats))],
                                " of processed reads")[[1]][1], ")", sep="")
        
        contams_before_plot <- readPNG(paste("1_fastqc1/",filename,"_fastqc/Images/sequence_length_distribution.png",sep=""))
        contams_after_plot <- readPNG(paste("5_fastqc2/",filename,"_fastqc/Images/sequence_length_distribution.png",sep=""))
        
        # plots
        textplot(paste("Overrepresented sequences:\n",contamsbefore," before pipeline, ",contamsafter," after pipeline.\n\n",trimmed_reads,"\n",too_short_reads,"\n\nBefore and after pipeline:",sep=""),valign="top",cex=1)
        rasterImage(contams_before_plot, 0,   0.1, 0.5, 0.75)
        rasterImage(contams_after_plot,  0.5, 0.1, 1,   0.75)
        title("\nContaminants Removal")
        box("figure")
        
        
        ### plot 4 - quality_filter ###
        
        # obtains quality filter stats
        quality_stats_tmp <- read.delim(paste("4_quality_filter/summaries/",filename,".summary",sep=""),header=F,stringsAsFactors=F)
        quality_stats <- capitalize(quality_stats_tmp[length(quality_stats_tmp[,1]),])
        
        # obtains plot for before and after
        quality_before_plot <- readPNG(paste("1_fastqc1/",filename,"_fastqc/Images/per_sequence_quality.png",sep=""))
        quality_after_plot <- readPNG(paste("5_fastqc2/",filename,"_fastqc/Images/per_sequence_quality.png",sep=""))
        
        # plots
        textplot(paste(quality_stats,"\nBefore and after pipeline:\n",sep=""),valign="top", cex=1.2)
        rasterImage(quality_before_plot, 0,   0.1, 0.5, 0.75)
        rasterImage(quality_after_plot,  0.5, 0.1, 1,   0.75)
        title("\nQuality Filter")
        box("figure")
        
        
        ### plot 5 - final_summary ###
        
        # build text for plot (from second fastqc stats)
        text_tmp <- after_fastqc_stats[grep(">>Basic Statistics",after_fastqc_stats[,1]):min(grep(">>END_MODULE",after_fastqc_stats[,1])),]
        text_readsafter <- text_tmp[grep("Filename",text_tmp[,1]):(length(text_tmp[,1])-1),]
        
        # calculates total reads removed
        reads_before <- as.numeric(before_fastqc_stats[min(grep("Total Sequences",before_fastqc_stats[,1])),2])
        reads_after <- as.numeric(after_fastqc_stats[min(grep("Total Sequences",after_fastqc_stats[,1])),2])
        reads_removed <- round( ((reads_before - reads_after) / reads_before)*100 ,digits=2)
        
        # plots
        textplot(text_readsafter,valign="top",cex=1.2, show.rownames=F, show.colnames=F)
        par(new=TRUE)
        textplot(paste("\nTOTAL READS REMOVED:\n","       ",reads_removed,"%",sep=""),cex=2)
        title("\nSummary\n(after processing)")
        par(new=FALSE)
        box("figure")
        
    }
    
    # closes pdf
    dev.off()
    
    
    ''')





###########################
###    Running lines    ###
###########################

## obtains instructions from command line, raises errors if not provided

if len(sys.argv) < 2: sys.exit("\nError:\n\tNo task name provided\n\tTo run all tasks in pipeline use 'all'\n\t\tE.G. python harryc_QC_pipeline.py all\n\tTo produce a graph of pipeline and task names use 'test'\n\t\tE.G. python harryc_DE_pipeline.py test\n")

if sys.argv[1] == "all": run = 'results2'
else: run = sys.argv[1]

# number of parallels can also be specified in the command line if want to override pipeline.ini
if len(sys.argv) > 2: parallels = int(sys.argv[2])
else: parallels = int(params['general_parallels'])


## runs pipeline

if run == "test":
    
    pipeline_printout_graph("pipeline_printout.svg", "svg", ['results2'], no_key_legend=False, minimal_key_legend=True, user_colour_scheme={"colour_scheme_index":6}, draw_vertically=True, pipeline_name='Quality Control Pipeline by Harry Clifford')
    
    pipeline_printout(sys.stdout, ['results2'], verbose = 4)
    
else:
    
    pipeline_printout_graph("pipeline_printout.svg", "svg", [run], no_key_legend = False, minimal_key_legend = True, user_colour_scheme = {"colour_scheme_index":6}, draw_vertically=True, pipeline_name='Quality Control Pipeline by Harry Clifford')
    
    pipeline_run([run], multiprocess=parallels, verbose = 4)
    
    pipeline_printout_graph("pipeline_printout.svg", "svg", [run], no_key_legend = False, minimal_key_legend = True, user_colour_scheme = {"colour_scheme_index":6}, draw_vertically=True, pipeline_name='Quality Control Pipeline by Harry Clifford')



###################################################################
###################################################################
###################################################################





