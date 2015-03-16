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
================
RNA-Seq Pipeline
================

Performs alignment of paired-end reads, duplicate removal, and differential expression/
Requires the following input files:
    1. Read files - pairs of fastq files (.fastq or .fastq.gz)
    2. A reference genome - single fasta file (.fa)
    3. A reference genome annotations - single gtf file (.gtf)
    4. The configuration file - pipeline.ini
    (5. If required - Sequences to be removed - single fasta file (.fa))


Overview
========

This pipeline will perform differential expression analysis between multiple tissues/conditions. which may have one or more biological replicate.

This pipeline performs the following:
    - aligns to sequences to be removed and retains unmapped using Bowtie2
    - perform an alignment for each sample using GSnap
    - remove duplicates using PicardTools
    - compute readcounts using HTSeq
    - perform differential expression analysis using DESeq and EdgeR
        - differential expression will occur through grouping of replicates, then comparisons between tissues/conditions

Ruffus output (telling you progress of the pipeline and any failings) will be directed to stdout. Outputs of the third-party software will be redirected to .log files in your working directory.

Usage
=====

Command line
------------

The command line used when running this pipeline should consist of: the python command, followed by the name of this file, followed by the stage you wish the pipeline to run until (if you would like to run the whole pipeline, use the argument "all"), followed by the number of parallel processes allowed to run at any one time (this is advised to be kept low if performing duplicate removal, as some of the temporary files produced can be 2-3 times as big as the standard bam files).

E.g.    python pipeline.py all 10


Read files
----------

Reads are imported through identification of the .fastq extension.
The default file format for these should be as follows:

    <tissue/condition>_<replicate>_<pair>.fastq


<tissue/condition> should correspond to the grouping of replicates
<replicate> should be 'R' followed by a numerical value
<pair> will be either '1' or '2' corresponding to the paired-end fastq files

For example 2 tissues with 2 replicates each will have 8 files:

    parietal_R1_1.fastq
    parietal_R1_2.fastq
    parietal_R2_1.fastq
    parietal_R2_2.fastq
    occipital_R1_1.fastq
    occipital_R1_2.fastq
    occipital_R2_1.fastq
    occipital_R2_2.fastq

This pipeline does not currently support gunzipped input files

Reference genome
----------------

Imported through identification of a file with a .fa extension


Reference genome annotations
----------------------------

Should be an Ensembl GTF file
Imported through identification of a file with a .gtf extension


Configuration file
------------------

The pipeline.ini file provided with this pipeline, which can be customized by the user before running the pipeline. Imported through identification of file names pipeline.ini.


Sequences to be removed
-----------------------

A file of sequences to be filtered out before Gsnap alignment (e.g. rRNA). Imported through identification of a file with a .remove.fa extension


Third-party software
====================

This pipeline requires the following third-party software to be installed:
    - Bowtie2
    - Gmap/Gsnap
    - PicardTools
    - Samtools
    - HTSeq
    - R packages - gplots, edgeR, DESeq


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
import errno
import re
import Pipeline
import glob
import HTSeq
import IOTools
import copy
from time import gmtime, strftime
from rpy2.robjects import r as R


# options from ini file
params = Pipeline.getParameters(["pipeline.ini"])


## startlog function
# opens log file for specific task, under name "LOGS/"+dirname+logname+".log"
def startlog(dirname,logname):
    # begins pipeline log if doesn't exist, and writes in logname and time
    makedirs("LOGS/"+dirname)
    if not os.path.exists("LOGS/"+dirname+"/"+logname+".log"):
        open("LOGS/"+dirname+"/"+logname+".log",'w').close()
    with open("LOGS/"+dirname+"/"+logname+".log",'a') as log:
        log.write('\n\n##### Starting - ' + logname +'\n### ' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + '\n\n')
        # (using 'with' makes 'close' function redundant)


## endlog function
# closes log file for specific task, under name "LOGS/"+logname+".log"
# (startlog function must have been run before)
def endlog(dirname,logname):
    with open("LOGS/"+dirname+"/"+logname+".log",'a') as log:
        log.write('\n\n##### Finished - ' + logname +'\n### ' + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + '\n\n')


## determine boolean from string function
def str2bool(string): return string in ["True","true","T","t","Yes","yes",1]

## create directory if does not exist function (bypasses race condition, but still raises any other errors)
def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise



###################################################################
## Sequence Removal
###################################################################


## bowtie indexing

@split(glob.glob("*.remove.fa")[0], regex(r"(.+)\.remove.fa"), [r'sequence_removal/bowtie_index/\1.index.*'])
def bowtie(infile,outfile):
    
    # makes directory if doesn't already exist
    makedirs("sequence_removal/bowtie_index")
    
    # starts pipeline log
    startlog('bowtie','bowtie')
    
    # builds bowtie indexing commandline
    command = params['sequence_removal_bowtie2_bashcommand'] + '-build'
    output = 'sequence_removal/bowtie_index/'+infile.rstrip('.remove.fa').split("/")[-1]+'.index'
    stdouterr = '>> LOGS/bowtie/bowtie.log 2>&1'
    
    # runs commandline
    os.system(" ".join([command,infile,output,stdouterr]))
    
    # ends pipeline log
    endlog('bowtie','bowtie')


## Example of command line produced:
# bowtie2-build rRNA.remove.fa sequence_removal/bowtie_index/rRNA.index >> LOGS/bowtie/bowtie.log 2>&1


###################################################################


## bowtie alignment and subsequent removal

@follows(bowtie)
@transform(glob.glob('*.fastq'), regex(r'(.+)\_1.fastq'), add_inputs(r'\1_2.fastq', r'sequence_removal/bowtie_index/*.index*'), [r'sequence_removal/outputs/\1_1.fastq', r'sequence_removal/outputs/\1_2.fastq'])
def sequence_removal(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("sequence_removal/outputs")
    
    # starts pipeline log
    logname = "_".join(infiles[0].split(".")[0].split("_")[:-1])
    startlog('sequence_removal','sequence_removal-'+logname)
    
    # builds bowtie commandline
    command = params['sequence_removal_bowtie2_bashcommand']
    outfile = "--un-conc " + re.sub(r'(.*).1', r'\1', outfiles[0])
    index = "-x " + infiles[2].split('.index', 1)[0] + '.index'
    infile1 = "-1 "+infiles[0]
    infile2 = "-2 "+infiles[1]
    throwaway = "-S /dev/null"
    stdouterr = ">> LOGS/sequence_removal/sequence_removal-" + logname + ".log 2>&1"
    
    # runs commandline
    os.system(" ".join([command,outfile,index,infile1,infile2,throwaway,stdouterr]))
    
    # renames outputs
    os.system("mv " + re.sub(r'(.*).1', r'\1.1', outfiles[0]) + " " + outfiles[0])
    os.system("mv " + re.sub(r'(.*).2', r'\1.2', outfiles[1]) + " " + outfiles[1])
    
    # ends pipeline log
    endlog('sequence_removal','sequence_removal-'+logname)


## Examples of command lines produced:
# bowtie2 --un-conc sequence_removal/outputs/nodose_R1.fastq -x trigeminal_R3_2.fastq.index -1 petrosal_R1_2.fastq -2 profundal_R3_1.fastq -S /dev/null >> LOGS/sequence_removal/sequence_removal-petrosal_R1.log 2>&1
# mv sequence_removal/outputs/nodose_R1.1.fastq sequence_removal/outputs/nodose_R1_1.fastq
# mv sequence_removal/outputs/nodose_R1.2.fastq sequence_removal/outputs/nodose_R1_2.fastq



###################################################################
## Alignment
###################################################################


## gmap build

def gmap_build(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("alignment/gmap")
    
    # starts pipeline log
    startlog('gmap','gmap')
    
    # build gmap commandline
    command = params['alignment_gmap_build_bashcommand']
    refgenname = "-d "+infiles.rstrip(".fa")
    stdouterr = '>> LOGS/gmap/gmap.log 2>&1'
    
    # runs commandline
    os.system(" ".join([command,refgenname,"-D alignment/gmap",infiles,stdouterr]))
    
    # ends pipeline log
    endlog('gmap','gmap')



# varying decorators for ruffus depending on sequence_removal true or false
if str2bool(params["general_sequence_removal"]):
    @follows(sequence_removal)
    @split([f for f in glob.glob('*.fa') if '.remove.fa' not in f], regex(r'(.+)\.fa'), r'alignment/gmap/\1/\1.*')
    def gmap(infiles, outfiles):
        gmap_build(infiles, outfiles)
else:
    @split([f for f in glob.glob('*.fa') if '.remove.fa' not in f], regex(r'(.+)\.fa'), r'alignment/gmap/\1/\1.*')
    def gmap(infiles, outfiles):
        gmap_build(infiles, outfiles)



## Example of command line produced:
# gmap_build -d Galgal.68.dna.toplevel -D alignment/gmap Galgal.68.dna.toplevel.fa >> LOGS/gmap/gmap.log 2>&1



###################################################################


## preliminary gsnap run

# function
def short_gsnap(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("alignment/preliminary_gsnap")
    
    # starts pipeline log
    logname = "_".join(infiles[0].split("/")[-1].split(".")[0].split("_")[:-1])
    startlog('preliminary_gsnap','preliminary_gsnap-'+logname)
    
    # build gsnap commandline
    refgenname = "-d " + [f for f in glob.glob('*.fa') if '.remove.fa' not in f][0].rstrip(".fa")
    pairexpect = "--pairexpect=" + str(params["alignment_long_exon_length"]/2)
    pairdev = "--pairdev=" + str(params["alignment_long_exon_length"]/2)
    samtools = "| " + params['general_samtools_bashcommand']  + " view -bS - >"
    stderr = "2>> LOGS/preliminary_gsnap/preliminary_gsnap-" + logname + ".log"
    
    # runs commandline
    os.system(" ".join([params['alignment_gsnap_bashcommand'], refgenname, "-D alignment/gmap", infiles[0], infiles[1], "-A sam --nthreads="+str(params['alignment_gsnap_nthreads']), pairexpect, pairdev, stderr, samtools, outfiles, stderr]))
    
    # ends pipeline log
    endlog('preliminary_gsnap','preliminary_gsnap-'+logname)


# varying decorators for ruffus depending on sequence_removal true or false
if str2bool(params["general_sequence_removal"]):
    @follows(sequence_removal)
    @follows(gmap)
    @transform(sequence_removal, regex(r'sequence_removal/outputs/(.+)_1.fastq'), add_inputs(r'sequence_removal/outputs/\1_2.fastq'), r'alignment/preliminary_gsnap/\1.bam')
    def preliminary_gsnap(infiles, outfiles):
        
        short_gsnap(infiles, outfiles)
        
else:
    @follows(gmap)
    @transform(glob.glob('*.fastq'), regex(r'(.+)_1.fastq'), add_inputs(r'\1_2.fastq'), r'alignment/preliminary_gsnap/\1.bam')
    def preliminary_gsnap(infiles, outfiles):
        
        short_gsnap(infiles, outfiles)



## Example of command line produced:
# gsnap -d Galgal.68.dna.toplevel -D alignment/gmap nodose_R1_1.fastq nodose_R1_2.fastq -A sam --nthreads=2 --pairexpect=500 --pairdev=500 2>> LOGS/preliminary_gsnap/preliminary_gsnap-nodose_R1.log | samtools view -bS - > alignment/preliminary_gsnap/nodose_R1.bam 2>> LOGS/preliminary_gsnap/preliminary_gsnap-nodose_R1.log



###################################################################


## convert gtf to gff3 ### MAY NEED IMPROVING - UNTESTED

@transform(glob.glob('*.gtf')[0], regex(r'(.+)\.gtf'), r'\1.gff3')
def gtf_to_gff3(infiles, outfiles):
    
    if len(glob.glob(infiles.replace("gtf","gff3"))) == 0:
        
        # starts pipeline log
        startlog('gtf_to_gff3','gtf_to_gff3')
        
        # redirects stdout and stderr before running main script below
        sys.stdout = open("LOGS/gtf_to_gff3/gtf_to_gff3.log", "a")
        sys.stderr = open("LOGS/gtf_to_gff3/gtf_to_gff3.log", "a")
        
        ### GTF2GFF function ### MAY NEED IMPROVING - UNTESTED
        # cycles through each line of gtf
        outstream = IOTools.openFile(outfiles, "w")
        gene = ""
        gtf = HTSeq.GFF_Reader(infiles)
        for line in gtf:
            splitline = line.get_gff_line().split("\t")
            ## if new gene
            if gene != line.attr["gene_id"]:
                # prints final details from old gene
                if len(gene) > 0:
                    outstream.write("\t".join(geneline)+"\n")
                    if len(str(mRNA).split("], [")) == 1:
                        components[0] = components[0][:-1] 
                        outstream.write("\t".join(mRNA)+"\n")
                        outstream.write("\t".join(components)+"\n")
                    else:
                        components[-1][0] = components[-1][0][:-1] 
                        for n in range(len(str(mRNA).split("], ["))):
                            outstream.write("\t".join(mRNA[n])+"\n")
                            outstream.write("\t".join(components[n])+"\n")
                # updates gene and transcript and starts new geneline (with gene name if available)
                gene = line.attr["gene_id"]
                transcript = line.attr["transcript_id"] 
                if "gene_name" in splitline[-1]: details = ["ID="+line.attr["gene_id"]+";Name="+line.attr["gene_name"]+";"]
                else: details = ["ID="+line.attr["gene_id"]+";"]
                geneline = splitline[0:2]+["gene"]+splitline[3:8]+details
                # creates mRNA line
                mRNA = splitline[0:2]+["mRNA"]+splitline[3:8]+["ID="+transcript+";Parent="+gene+";"]
                # creates first component line
                components = ["\t".join(splitline[:-1]+["ID="+splitline[1]+":"+transcript+":"+line.attr["exon_number"]+";Parent="+gene+";"])+"\n"]
            ## if same gene
            else:
                ## if new mRNA
                if transcript != line.attr["transcript_id"]:
                    # updates transcript
                    transcript = line.attr["transcript_id"]
                    # updates geneline based on new component line
                    for i in [0,1,5,6,7]:
                        if splitline[i] != geneline[i]: geneline[i] = "."
                    if splitline[3] < geneline[3]: geneline[3] = splitline[3]
                    if splitline[4] > geneline[4]: geneline[4] = splitline[4]
                    # additional mRNA and component line creations
                    if len(str(mRNA).split("], [")) == 1:
                        mRNA = [mRNA]+[splitline[0:2]+["mRNA"]+splitline[3:8]+["ID="+transcript+";Parent="+gene+";"]]
                        components[0] = components[0][:-1]                
                        components = [components]+[["\t".join(splitline[:-1]+["ID="+splitline[1]+":"+transcript+":"+line.attr["exon_number"]+";Parent="+gene+";"])+"\n"]]
                    else:
                        mRNA = [i for i in mRNA] + [splitline[0:2]+["mRNA"]+splitline[3:8]+["ID="+transcript+";Parent="+gene+";"]]
                        components[-1][0] = components[-1][0][:-1]
                        components = [i for i in components] + [["\t".join(splitline[:-1]+["ID="+splitline[1]+":"+transcript+":"+line.attr["exon_number"]+";Parent="+gene+";"])+"\n"]]
                ## if same mRNA
                else:
                    # updates geneline and mRNA based on new component line:
                    # introduces "." into slots with multiple possible answers
                    for i in [0,1,5,6,7]:
                        if splitline[i] != geneline[i]:
                            geneline[i] = "."
                            if i == 0: geneline[3] = "." ; geneline[4] = "."
                        if len(str(mRNA).split("], [")) == 1:
                            if splitline[i] != mRNA[i]:
                                mRNA[i] = "."
                                if i == 0: mRNA[3] = "." ; mRNA[4] = "."
                        if len(str(mRNA).split("], [")) > 1:
                            if splitline[i] != mRNA[-1][i]:
                                mRNA[-1][i] = "."
                                if i == 0: mRNA[-1][3] = "." ; mRNA[-1][4] = "."
                    # changes position slots if new component positions lie outside of these
                    if splitline[3] < geneline[3] and geneline[0] != ".": geneline[3] = splitline[3]
                    if len(str(mRNA).split("], [")) == 1:
                        if splitline[3] < mRNA[3] and mRNA[0] != ".": mRNA[3] = splitline[3]
                    if len(str(mRNA).split("], [")) > 1:
                        if splitline[3] < mRNA[-1][3] and mRNA[-1][0] != ".": mRNA[-1][3] = splitline[3]
                    if splitline[4] > geneline[4] and geneline[0] != ".": geneline[4] = splitline[4]
                    if len(str(mRNA).split("], [")) == 1:
                        if splitline[4] > mRNA[4] and mRNA[0] != ".": mRNA[4] = splitline[4]
                    if len(str(mRNA).split("], [")) > 1:
                        if splitline[4] > mRNA[-1][4] and mRNA[-1][0] != ".": mRNA[-1][4] = splitline[4]  
                    # adds to current component line
                    if len(str(components).split("], [")) == 1: components = [components[0]+"\t".join(splitline[:-1]+["ID="+splitline[1]+":"+line.attr["transcript_id"]+":"+line.attr["exon_number"]+";Parent="+gene+";"])+"\n"]
                    else: components[-1] = ["".join(components[-1]+["\t".join(splitline[:-1]+["ID="+splitline[1]+":"+line.attr["transcript_id"]+":"+line.attr["exon_number"]+";Parent="+gene+";"])+"\n"])]
        # prints final gene
        outstream.write("\t".join(geneline)+"\n")
        if len(str(mRNA).split("], [")) == 1: outstream.write("\t".join(mRNA)+"\n") ; outstream.write("\t".join(components)+"\n")
        else:
            for n in range(len(str(mRNA).split("], ["))): outstream.write("\t".join(mRNA[n])+"\n") ; outstream.write("\t".join(components[n])+"\n") 
        
        
        # restores stdout/stderr
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        
        # ends pipeline log
        endlog('gtf_to_gff3','gtf_to_gff3')




###################################################################

## extract long exons

@transform(gtf_to_gff3, regex(r'(.+)\.gff3'), "alignment/determine_insert_size/"+glob.glob("*.gtf")[0].rstrip(".gtf")+"3.min_"+str(params['alignment_long_exon_length'])+".const_exons.gff")
def extract_long_exons(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("alignment/determine_insert_size/outputs")
    
    # starts pipeline log
    startlog('extract_long_exons','extract_long_exons')
    
    # build miso commandline
    command = params['alignment_exon_utils_bashcommand']
    min_exon = "--min-exon-size " + str(params['alignment_long_exon_length'])
    stdouterr = '>> LOGS/extract_long_exons/extract_long_exons.log 2>&1'
    
    # runs commandline
    os.system(" ".join([command, "--get-const-exons", infiles, min_exon, "--output-dir alignment/determine_insert_size", stdouterr]))
    
    # ends pipeline log
    endlog('extract_long_exons','extract_long_exons')


## Example of command line produced:
# python /net/isi-software/src/misopy-0.4.5/misopy/exon_utils.py --get-const-exons Galgal.68.gff3 --min-exon-size 1000 --output-dir alignment/determine_insert_size >> LOGS/extract_long_exons/extract_long_exons.log 2>&1



###################################################################


## determine insert length

@follows(extract_long_exons)
@follows(preliminary_gsnap)
@transform(preliminary_gsnap, regex(r'alignment/preliminary_gsnap/(.+).bam'), add_inputs(r'alignment/determine_insert_size/*.gff'), r'alignment/determine_insert_size/outputs/\1.bam.insert_len')
def determine_insert_size(infiles, outfiles):
    
    # starts pipeline log
    logname = infiles[0].split("/")[-1].split(".")[0]
    startlog('determine_insert_size','determine_insert_size-'+logname)
    
    # runs miso commandline
    command = params['alignment_pe_utils_bashcommand']
    stdouterr = " >> LOGS/determine_insert_size/determine_insert_size-" + logname + ".log 2>&1"
    os.system(" ".join([command, "--compute-insert-len", infiles[0], infiles[1], "--output-dir alignment/determine_insert_size/outputs", stdouterr]))
    
    # removes resultant bam file that is not required
    rmfilename = glob.glob('alignment/determine_insert_size/*.const_exons.gff')[0].replace('determine_insert_size/','determine_insert_size/outputs/bam2gff_')+'/'+logname+'.bam'
    os.system('rm '+rmfilename+stdouterr)
    
    # ends pipeline log
    endlog('determine_insert_size','determine_insert_size-'+logname)



## Examples of command lines produced:
# python /net/isi-software/src/misopy-0.4.5/misopy/pe_utils.py --compute-insert-len alignment/preliminary_gsnap/nodose_R1.bam alignment/determine_insert_size/Galgal.683.min_1000.const_exons.gff --output-dir alignment/determine_insert_size/outputs >> LOGS/determine_insert_size/determine_insert_size-nodose_R1.log 2>&1
# 'rm -r alignment/determine_insert_size/outputs/bam2gff_Galgal.683.min_1000.const_exons.gff/nodose_R1.bam >> LOGS/determine_insert_size/determine_insert_size-nodose_R1.log 2>&1'



###################################################################


## gsnap alignment


# input objects decided based on whether sequence_removal was performed

if str2bool(params["general_sequence_removal"]):
    alignment_fastqinputs = ["sequence_removal/outputs/" + i for i in glob.glob("*.fastq")]
    alignment_regex = regex(r'sequence_removal/outputs/(.+)_1.fastq')
    alignment_add_inputs = r'sequence_removal/outputs/\1_2.fastq'
else:
    alignment_fastqinputs = glob.glob('*.fastq')
    alignment_regex = regex(r'(.+)_1.fastq')
    alignment_add_inputs = r'\1_2.fastq'


# decorators and functions defined depending on whether determine_insert_size was performed

if str2bool(params["alignment_determine_insert_size"]):
    @follows(determine_insert_size)
    @transform(alignment_fastqinputs, alignment_regex, add_inputs(alignment_add_inputs, r'alignment/determine_insert_size/outputs/\1.bam.insert_len'), r'alignment/gsnap/\1.bam')
    def gsnap(infiles, outfiles):
        
        # makes directory if doesn't already exist
        makedirs("alignment/gsnap")
        
        # starts pipeline log
        logname = "_".join(infiles[0].split("/")[-1].split(".")[0].split("_")[:-1])
        startlog('gsnap','gsnap-'+logname)
        
        # obtains insert length from miso results
        insertline = open(infiles[2],"r").readline()
        exec insertline.split(",")[0].replace("#","")
        exec insertline.split(",")[1].replace("#","")
        mean = str(round(mean)).split(".")[0]
        sdev = str(round(sdev)).split(".")[0]
        
        # builds gsnap commandline
        refgenname = "-d " + [f for f in glob.glob('*.fa') if '.remove.fa' not in f][0].rstrip(".fa")
        pairexpect = "--pairexpect=" + str(mean)
        pairdev = "--pairdev=" + str(sdev)
        stderr = "2>> LOGS/gsnap/gsnap-" + logname + ".log"
        samtools = "| " + params['general_samtools_bashcommand'] + " view -bS - >"
        
        # runs commandline
        os.system(" ".join([params['alignment_gsnap_bashcommand'], refgenname, "-D alignment/gmap", infiles[0], infiles[1], "-n 1 -A sam --nthreads="+str(params['alignment_gsnap_nthreads']), pairexpect, pairdev, "--novelsplicing=1", stderr, samtools, outfiles, stderr]))
        
        # ends pipeline log
        endlog('gsnap','gsnap-'+logname)
        
        
else:
    @follows(gmap)
    @transform(alignment_fastqinputs, alignment_regex, add_inputs(alignment_add_inputs), r'alignment/gsnap/\1.bam')
    def gsnap(infiles, outfiles):
        
        # makes directory if doesn't already exist
        makedirs("alignment/gsnap")
        
        # starts pipeline log
        logname = "_".join(infiles[0].split("/")[-1].split(".")[0].split("_")[:-1])
        startlog('gsnap','gsnap-'+logname)
        
        # obtains insert length from pipeline.ini
        pairexpect = "--pairexpect=" + str(params['alignment_average_insert'])
        pairdev = "--pairdev=" + str(params['alignment_insert_dev'])
        
        # builds gsnap commandline
        refgenname = "-d " + [f for f in glob.glob('*.fa') if '.remove.fa' not in f][0].rstrip(".fa")
        stderr = "2>> LOGS/gsnap/gsnap-" + logname + ".log"
        samtools = "| " + params['general_samtools_bashcommand'] + " view -bS - >"
        
        # runs commandline
        os.system(" ".join([params['alignment_gsnap_bashcommand'], refgenname, "-D alignment/gmap", infiles[0], infiles[1], "-n 1 -A sam --nthreads="+str(params['alignment_gsnap_nthreads']), pairexpect, pairdev, "--novelsplicing=1", stderr, samtools, outfiles, stderr]))
        
        # ends pipeline log
        endlog('gsnap','gsnap-'+logname)



## Examples of command lines produced:
# gsnap -d Galgal.68.dna.toplevel -D alignment/gmap nodose_R1_1.fastq nodose_R1_2.fastq -A sam --nthreads=2 --pairexpect=116 --pairdev=72 --novelsplicing=1 2>> LOGS/gsnap/gsnap-nodose_R1.log | samtools view -bS - > alignment/gsnap/nodose_R1.bam 2>> LOGS/gsnap/gsnap-nodose_R1.log




###################################################################
## Post-Alignment Processes
###################################################################


## alignment stats

@follows(gsnap)
@transform(gsnap, regex(r'alignment/gsnap/(.+).bam'), r'alignment/alignment_stats/alignment_stats_\1.txt')
def alignment_stats(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("alignment/alignment_stats")
    
    # starts pipeline log
    logname = infiles.split("/")[-1].split(".")[0]
    startlog('alignment_stats','alignment_stats-'+logname)
    
    # runs commandline
    command = params['general_samtools_bashcommand']
    stderr = "2>> LOGS/alignment_stats/alignment_stats-" + logname + ".log"
    os.system(" ".join([command,"flagstat",infiles,">",outfiles,stderr]))
    
    # ends pipeline log
    endlog('alignment_stats','alignment_stats-'+logname)



## Examples of command lines produced:
# samtools flagstat alignment/gsnap/nodose_R1.bam > alignment/alignment_stats/alignment_stats_nodose_R1.txt 2>> LOGS/alignment_stats/alignment_stats-nodose_R1.log




###################################################################

## duplicate removal


## functions

# ubiquitous removal and sort bams
def ubiquitous_renaming(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("duplicate_removal/sorted_ubiquitous_renamed_bams")
    
    # starts pipeline log
    logname = infiles.split("/")[-1].split(".")[0]
    startlog('ubiquitous_renaming','ubiquitous_renaming-'+logname)
    
    # builds temporary python script
    # (built as separate script, so can be used in bash command line, and piped directly into samtools view)
    tmp_script = """# imports relevant modules
import sys
# function for determining if one read matches another
def match(read1,read2):
    if read1.split("\t")[2] == read2.split("\t")[6] and read2.split("\t")[2] == read1.split("\t")[6]: return True
    else: return False
## reads in data
inputfile = sys.argv[1]
env = open(inputfile,"r")
# prints any headers
datastream = env.readline()
while "@SQ" in datastream:
    print datastream.strip()
    datastream = env.readline()
readname = datastream.split("\t",1)[0]
readvector = []
while True:
    ### if same readname, appends read to old readvector
    if readname == datastream.split("\t",1)[0]:
        readvector.append(datastream.split("\t",1)[1])
    ### if new readname, processes old readvector, then creates new
    else:
        ## Processing - if only 2 reads (or even 1 read), does not rename, prints
        if len(readvector) <= 2:
            print "".join([readname+"\t"+i for i in readvector]).rstrip()
        ## Processing - if greater than 2 reads, matches reads and pops out of vector until less than 2 left
        else:
            finalreadvector  = []
            while True:
                oldlength = len(readvector)
                # tries to find a match with first read, and pops matching reads out
                for n in range(1,len(readvector)):
                    if match(readvector[0],readvector[n]):
                        finalreadvector.extend([readname+"\t"+readvector.pop(0)])
                        finalreadvector.extend([readname+"\t"+readvector.pop(n-1)])
                        break
                # if length has not changed (no pop occurred), just deletes first read
                if len(readvector) == oldlength: del readvector[0]
                # breaks loop only when less than 2 reads remain
                if len(readvector) < 2: break
            # renames those in the finalreadvector, then prints
            if len(finalreadvector) > 0:
                for i in range((len(finalreadvector)-2)/2):
                    indices = range(0,len(finalreadvector),2)[1:]
                    finalreadvector[indices[i]] = str(i+2)+finalreadvector[indices[i]]
                    finalreadvector[indices[i]+1] = str(i+2)+finalreadvector[indices[i]+1]
            print "".join(finalreadvector).rstrip()
        ## creates new readname and readvector
        readname = datastream.split("\t",1)[0]
        readvector = [datastream.split("\t",1)[1]]
    ### next line becomes datastream
    datastream = env.readline()
    if datastream == "": break"""
    
    with open("duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_script_"+logname+".py",'w') as tmp: tmp.write(tmp_script)
    
    ## Produce main command line
    command = params['general_samtools_bashcommand']
    stderr = "2>> LOGS/ubiquitous_renaming/ubiquitous_renaming-" + logname + ".log"
    stdouterr = ">> LOGS/ubiquitous_renaming/ubiquitous_renaming-" + logname + ".log 2>&1"
    
    # first sort command
    sort_cmd1 = " ".join([command, 'sort -no', infiles, 'duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_' +logname, stderr])
    # first view command (bam to sam)
    view_cmd1 = " ".join([command, 'view -h - > duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_' +logname +'.sam', stderr])
    # python command
    py_cmd = " ".join([params['duplicate_removal_python_bashcommand'], 'duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_script_' +logname +'.py', 'duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_' +logname +'.sam', stderr])
    # second view command (sam to bam)
    view_cmd2 = " ".join([command, 'view -bS - > duplicate_removal/sorted_ubiquitous_renamed_bams/' +logname +'.bam', stderr])
    # second sort command
    sort_cmd2 = " ".join([command, 'sort duplicate_removal/sorted_ubiquitous_renamed_bams/' +logname +'.bam duplicate_removal/sorted_ubiquitous_renamed_bams/' +logname, stdouterr])
    
    
    ## Runs complete commandlines
    
    # sorts and creates temporary sam
    os.system(" | ".join([sort_cmd1, view_cmd1]))
    
    # runs python script for renaming ubiquitous reads, followed by a sort
    os.system(" | ".join([py_cmd, view_cmd2]))
    os.system(sort_cmd2)
    
    # deletes temporary script and temporary sam
    os.system("rm duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_script_"+logname+".py "+stdouterr)
    os.system("rm duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_"+logname+".sam "+stdouterr)
    
    # ends pipeline log
    endlog('ubiquitous_renaming','ubiquitous_renaming-'+logname)



# Examples of command lines produced:
# samtools sort -no duplicate_removal/sorted_ubiquitous_renamed_bams/nodose_R1.bam duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_nodose_R1 2>> LOGS/ubiquitous_renaming/ubiquitous_renaming-nodose_R1.log | samtools view -h - > duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_nodose_R1.sam 2>> LOGS/ubiquitous_renaming/ubiquitous_renaming-nodose_R1.log
# python duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_script_nodose_R1.py duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_nodose_R1.sam 2>> LOGS/ubiquitous_renaming/ubiquitous_renaming-nodose_R1.log | samtools view -bS - > duplicate_removal/sorted_ubiquitous_renamed_bams/nodose_R1.bam 2>> LOGS/ubiquitous_renaming/ubiquitous_renaming-nodose_R1.log
# samtools sort duplicate_removal/sorted_ubiquitous_renamed_bams/nodose_R1.bam duplicate_removal/sorted_ubiquitous_renamed_bams/nodose_R1 >> LOGS/ubiquitous_renaming/ubiquitous_renaming-nodose_R1.log 2>&1
# rm duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_script_nodose_R1.py 2>> LOGS/ubiquitous_renaming-nodose_R1.log
# rm duplicate_removal/sorted_ubiquitous_renamed_bams/tmp_nodose_R1.sam 2>> LOGS/ubiquitous_renaming-nodose_R1.log



# picardtools
def picardtools(infiles, outfiles):
    
    # makes directories if don't already exist
    makedirs("duplicate_removal/outputs")
    makedirs("duplicate_removal/reports")
    
    # starts pipeline log
    logname = infiles.split("/")[-1].split(".")[0]
    startlog('duplicate_removal','duplicate_removal-'+logname)
    
    # runs commandline
    command = params['duplicate_removal_picard_bashcommand']
    inputs = "INPUT="+infiles
    outputs = "OUTPUT="+outfiles
    metrics = "METRICS_FILE=duplicate_removal/reports/"+logname+".report.txt"
    max_files = "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP="+str(params['duplicate_removal_max_filehandles'])
    stdouterr = ">> LOGS/duplicate_removal/duplicate_removal-" + logname + ".log 2>&1"
    os.system(" ".join([command,inputs,outputs,metrics,"REMOVE_DUPLICATES=true ASSUME_SORTED=true",max_files,stdouterr]))
    
    # ends pipeline log
    endlog('duplicate_removal','duplicate_removal-'+logname)



# Example of command line produced:
# java -jar /net/isi-software/src/picard-tools-1.49/MarkDuplicates.jar INPUT=duplicate_removal/sorted_ubiquitous_renamed_bams/nodose_R1.bam OUTPUT=duplicate_removal/outputs/nodose_R1.bam METRICS_FILE=duplicate_removal/reports/nodose_R1.report.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 >> LOGS/duplicate_removal/duplicate_removal-nodose_R1.log 2>&1




## main function

@follows(alignment_stats)
@transform(gsnap, regex(r'alignment/gsnap/(.+).bam'), r'duplicate_removal/outputs/\1.bam')
def duplicate_removal(infiles,outfiles):
    
    tmp_filename = 'duplicate_removal/sorted_ubiquitous_renamed_bams/' + infiles.split("/")[-1]
    
    ubiquitous_renaming(infiles,tmp_filename)
    
    picardtools(tmp_filename,outfiles)




###################################################################


## duplicate removal stats

@follows(duplicate_removal)
@transform(duplicate_removal, regex(r'duplicate_removal/outputs/(.+).bam'), r'duplicate_removal/updated_alignment_stats/updated_alignment_stats_\1.txt')
def updated_alignment_stats(infiles, outfiles):
    
    # makes directory if doesn't already exist
    makedirs("duplicate_removal/updated_alignment_stats")
    
    # starts pipeline log
    logname = infiles.split("/")[-1].split(".")[0]
    startlog('updated_alignment_stats','updated_alignment_stats-'+logname)
    
    # runs commandline
    command = params['general_samtools_bashcommand']
    stderr = "2>> LOGS/updated_alignment_stats/updated_alignment_stats-" + logname + ".log"
    os.system(" ".join([command,"flagstat",infiles,">",outfiles,stderr]))
    
    # ends pipeline log
    endlog('updated_alignment_stats','updated_alignment_stats-'+logname)



# Example of command line produced:
# samtools flagstat duplicate_removal/outputs/nodose_R1.bam > duplicate_removal/updated_alignment_stats/nodose_R1.stats.txt 2>> LOGS/updated_alignment_stats/updated_alignment_stats-nodose_R1.log




###################################################################

## HTSeq


## functions

def htseq(infiles,outfiles):
    
    # makes directory if doesn't already exist
    makedirs("RESULTS/htseq_counts")
    
    # starts pipeline log
    logname = infiles.split("/")[-1].split(".")[0]
    startlog('htseq_counts','htseq_counts-'+logname)
    
    # runs commandline
    
    stderr = "2>> LOGS/htseq_counts/htseq_counts-" + logname + ".log"
    
    view_cmd = " ".join([params['general_samtools_bashcommand'],'view -h',infiles,stderr])
    sort_cmd = " ".join([params['HTSeq_sort_bashcommand'],'-s -k 1,1',stderr])
    htseq_cmd = " ".join([params['HTSeq_htseq_bashcommand'],'-q -s no -m intersection-strict - ',glob.glob('*.gtf')[0],'>',outfiles,stderr])
    
    os.system(" | ".join([view_cmd,sort_cmd,htseq_cmd]))
    
    # ends pipeline log
    endlog('htseq_counts','htseq_counts-'+logname)


# Example of command line produced:
# samtools view -h nodose_R1.bam 2>> LOGS/htseq_count/htseq_count-nodose_R1.log | sort -s -k 1,1 2>> LOGS/htseq_count/htseq_count-nodose_R1.log | htseq-count -q -s no -m intersection-strict -  Gallus_gallus.WASHUC2.70.gtf > nodose_R1_counts.txt 2>> LOGS/htseq_count/htseq_count-nodose_R1.log




def htseq2(infiles,outfiles):
    
    # starts pipeline log
    startlog('htseq_merge','htseq_merge')
    
    # finds files and changes working directories
    R('''
    # finds files
    origwd <- getwd()
    setwd('RESULTS/htseq_counts')
    Rinfiles <- list.files(pattern='*.txt')
    ''')
    
    # reads in files and merges
    R('''
    # reads in first file
    results <- read.delim(Rinfiles[1],header=F,stringsAsFactors=F,col.names=c("gene_id",gsub('.txt','',Rinfiles[1])))
    results <- results[1:(nrow(results)-5),]
    # reads in and merges rest of files
    for(file in Rinfiles[2:length(Rinfiles)]){
        tmp_results <- read.delim(file,header=F,stringsAsFactors=F,col.names=c("gene_id",gsub('.txt','',file[1])))
        tmp_results <- tmp_results[1:(nrow(tmp_results)-5),]
        results <- merge(results,tmp_results)
    }
    ''')
    
    # outputs resulting table and restores working directories
    R('''
    colnames(results) <- gsub("[.]","-",colnames(results))
    write.table(results, file='ALL_COUNTS.txt', quote=F, sep="\t", row.names=F)
    setwd(origwd)
    ''')
    
    # ends pipeline log
    endlog('htseq_merge','htseq_merge')





## varying decorators for ruffus depending on duplicate_removal true or false

if str2bool(params["general_duplicate_removal"]):
    @follows(updated_alignment_stats)
    @transform(duplicate_removal, regex(r'duplicate_removal/outputs/(.+).bam'), r'RESULTS/htseq_counts/\1.txt')
    def htseq_count(infiles,outfiles):
        
        htseq(infiles,outfiles)
        
else:
    @follows(alignment_stats)
    @transform(gsnap, regex(r'alignment/gsnap/(.+).bam'), r'RESULTS/htseq_counts/\1.txt')
    def htseq_count(infiles,outfiles):
        
        htseq(infiles,outfiles)


@follows(htseq_count)
@merge(htseq_count, 'RESULTS/htseq_counts/ALL_COUNTS.txt')
def htseq_merge(infiles,outfiles):
    
    htseq2(infiles,outfiles)




###################################################################
## Differential Expression Analysis
###################################################################


## DESeq


# DESeq's output is an EdgeR output, to ensure both are always run together
@follows(htseq_merge)
@split(htseq_merge,['RESULTS/deseq/*.png','RESULTS/deseq/*.png', 'RESULTS/deseq/deseq_R_Workspace.RData'])
def deseq(infile,outfile):
    
    # makes directory if doesn't already exist
    makedirs("RESULTS/deseq")
    
    # starts pipeline log
    startlog('deseq','deseq')
    
    # loads expansive level and FDR from pipeline.ini
    expansive_level = params['DESeq_EdgeR_expansive']
    FDR_threshold = params['DESeq_EdgeR_fdr']
    R('''expansive <- as.integer(%(expansive_level)s)'''% locals())
    R('''FDR <- %(FDR_threshold)s'''% locals())
    
    
    ## DESeq R script
    
    # redirects stdout and stderr before running any R
    sys.stdout = open("LOGS/deseq/deseq.log", "a")
    sys.stderr = open("LOGS/deseq/deseq.log", "a")
    
    # loads libraries and files, and changes working directory
    R('''
    # load DESeq
    library("DESeq")
    library("gplots")
    # read in data
    countTable <- read.delim('RESULTS/htseq_counts/ALL_COUNTS.txt',header=T,row.names=1,stringsAsFactors=F)
    origwd <- getwd()
    setwd('RESULTS/deseq')
    ''')
    
    # initial DESeq analysis
    R('''
    ### Initial analysis (only requires one run, does not need to be looped)
    # renames condition names if likely to cause clash later
    if("excluded" %in% sapply(strsplit(colnames(countTable),"_",fixed=T),"[",1)){
        colnames(countTable) <- gsub("excluded","excluded1",colnames(countTable))
    }
    ## create metadata
    Design <- data.frame(
        row.names = colnames(countTable),
        condition = sapply(strsplit(colnames(countTable),"_",fixed=T),"[",1)
        )
    ## process data in DESeq
    condition <- Design$condition
    cds <- newCountDataSet(countTable,condition)
    cds <- estimateSizeFactors(cds)
    ''')
    
    # writes intial analysis to files
    R('''
    ## write size factors to file
    sizefactors <- sizeFactors(cds)
    write.table(sizefactors, file="size_factors.txt", quote=F, sep="\t",col.names=F)
    
    ## plot size factors
    png(file="size_factors.png",height=650,width=650)
    margin_extension <- max(strwidth(names(sizeFactors(cds)),"inches")) + 0.3
    par(mai=c(margin_extension,0.82,0.82,0.42)) # mai to fit bottom labels
    
    barplot(sizeFactors(cds),main="Size Factors",las=2)
    
    par(mai=c(1.02,0.82,0.82,0.42)) # reset back to default mai
    invisible(dev.off())
    
    ## write normalized counts to file
    normalized_cds <- counts(cds,normalized=T)
    write.table(normalized_cds, file="normalized_counts.txt", quote=F, sep="\t")
    ''')
    
    # second initial DESeq analysis
    R('''
    ## variance stabilization
    
    # change of methods if only one replicate
    if(length(condition) == length(unique(condition))){
        cds <- estimateDispersions(cds,method="blind",sharingMode="fit-only")
        }else{
        cds <- estimateDispersions(cds)
    }
    
    vsd <- getVarianceStabilizedData(cds)
    pcaResult<-prcomp(t(vsd))
    dists <- dist(t(vsd))
    ''')
    
    # writes second analysis to files
    R('''
    ## write variance stabilized counts to file
    write.table(vsd, file="variance_stabilized_counts.txt", quote=F, sep="\t")
    
    ## plot heatmap
    
    png(file="heatmap_distances_allgenes.png",height=650,width=650)
    margin_extension <- max(strwidth(names(sizeFactors(cds)),"inches")) - 0.4
    if(margin_extension < 0){ margin_extension <- 0 }
    par(omi=c(margin_extension,0.1,0.1,margin_extension)) # omi to fit labels
    
    heatmap.2(as.matrix(dists), trace="none")
    
    par(omi=c(0,0,0,0)) # reset back to default omi
    invisible(dev.off())
    
    ## plot PCAs
    
    # function
    PCAplots <- function(prcomp_output, condition, firstPC=1, secondPC=2){
        
        prcomp_subset <- prcomp_output$x[,c(firstPC,secondPC)]
        
        margin_extension <- max(strwidth(as.vector(condition),"inches"))
        png(file=paste("Raw_PCA_PC",firstPC,"PC",secondPC,"_allgenes.png",sep=""),height=650,width=650)
        
        plot(prcomp_subset,
            main="Variance Stabilized Principal Component Analysis",
            xlim=c(min(prcomp_subset[,1])*1.1,max(prcomp_subset[,1])*1.1),
            ylim=c(min(prcomp_subset[,2])*1.1,max(prcomp_subset[,2])*1.1),
            pch=20,
            xlab=paste("PC",firstPC,"   (",round(summary(prcomp_output)$importance[2,firstPC]*100),"% of Variance)",sep=""),
            ylab=paste("PC",secondPC,"   (",round(summary(prcomp_output)$importance[2,secondPC]*100),"% of Variance)",sep="")
            )
        text(prcomp_subset, rownames(prcomp_subset), pos=3, cex=0.9)
        
        invisible(dev.off())
        
        png(file=paste("Coloured_PCA_PC",firstPC,"PC",secondPC,"_allgenes.png",sep=""),height=650,width=650)
        
        cols <- as.vector(condition)
        for(i in 1:length(unique(condition))){cols[which(cols==unique(condition)[i])] <- rainbow(length(unique(condition)))[i]}
        
        par(omi=c(0.1, 0.1, 0.1, margin_extension + 0.9),xpd=NA) # omi to fit legend
        plot(prcomp_subset,
            main="Variance Stabilized Principal Component Analysis",
            pch=19, col=cols,
            xlab=paste("PC",firstPC,"   (",round(summary(prcomp_output)$importance[2,firstPC]*100),"% of Variance)",sep=""),
            ylab=paste("PC",secondPC,"   (",round(summary(prcomp_output)$importance[2,secondPC]*100),"% of Variance)",sep="")
            )
        
        legend(par("usr")[2],
            par("usr")[4],
            as.vector(unique(condition)),
            cex=0.9, pch=19,
            col=rainbow(length(unique(condition)))
            )
        
        par(omi=c(0,0,0,0),xpd=FALSE) # reset back to default omi
        invisible(dev.off())
    }
    
    # runs function
    PCAplots(pcaResult, condition, firstPC=1, secondPC=2)
    PCAplots(pcaResult, condition, firstPC=2, secondPC=3)
    PCAplots(pcaResult, condition, firstPC=3, secondPC=4)
    ''')
    
    # Generates group listings for pairwise comparisons (dependant on expansive level)
    R('''
    ### Generate groups to compare (will be attached to metadata)
    
    samples <- sapply(strsplit(colnames(countTable),"_",fixed=T),"[",1)
    types <- unique(samples)
    groups <- list()
    
    ## expansive = 1. Simply compares every single tissue/condition pairwise
    if(expansive >= 1 & length(types > 1)){
        # takes all unique pairwise permutations
        perms <- t(combn(types,2))
        # cycles through pairs, changing names to numbers and appending to list
        for(i in 1:length(perms[,1])){
            tmp_samples <- samples
            tmp_samples[which(samples !=perms[i,1] & samples !=perms[i,2])] <- "excluded"
            groups[[i]] <- tmp_samples
        }
    }
    
    ## expansive = 2. Does above, plus compares each tissue/condition versus all others together
    if(expansive >= 2 & length(types > 2)){
        # cycles through each individual, changing names to numbers and appending to list
        for(i in 1:length(types)){
            tmp_samples <- samples
            tmp_samples[which(tmp_samples!=types[i])] <- "others"
            groups[[length(groups)+1]] <- tmp_samples
        }
    }
    
    ## expansive = 3. Does all of above, plus compares all different possible groupings with one another
    if(expansive >= 3 & length(types > 3)){
        # first loop generates group permutations for each group size
        for(i in 2:floor(length(types)/2)){
            tmp_perms <- t(combn(types,i))
            # second loop cycles through group permutations, changing names to numbers and appending to list
            for(j in 1:length(tmp_perms[,1])){
                tmp_samples <- samples
                tmp_samples[which(samples %in% tmp_perms[j,])] <- paste(tmp_perms[j,],collapse="_")
                tmp_samples[which(! samples %in% tmp_perms[j,])] <- paste(types[which(! types %in% tmp_perms[j,])],collapse="_")
                groups[[length(groups)+1]] <- tmp_samples
            }
        }
    }
    ''')
    
    # defines a function for capitalising strings (for use in graphs)
    R('''    
    # function for use in upcoming loop
    cap <- function(x){
        s <- strsplit(x, "_")[[1]]
        paste(toupper(substring(s,1,1)),substring(s,2),sep="",collapse=" + ")
    }
    ''')
    
    # performs pairwise DESeq analysis on each group (decided above based on expansive level)
    R('''
    ### loop through groups as decided above
    for(groupvector in groups){
        
        # determine groups to compare
        tmp_conditions <- as.vector(unique(groupvector))
        if("excluded" %in% tmp_conditions){ tmp_conditions <- tmp_conditions[-which(tmp_conditions=="excluded")] }
        
        # creates directory for results and changes working directory
        tmp_origwd <- getwd()
        tmp_comparison <- paste(tmp_conditions[1],"VS",tmp_conditions[2],sep="_")
        dir.create(tmp_comparison, showWarnings = FALSE)
        setwd(tmp_comparison)
        
        ## create metadata
        Design <- data.frame(
            row.names = colnames(countTable),
            condition = groupvector
            )
        
        ## process data in DESeq
        condition <- Design$condition
        cds <- newCountDataSet(countTable,condition)
        cds <- estimateSizeFactors(cds)
        
        # change of methods if only one replicate
        if("excluded" %in% condition){
            tmp_reps <- condition[-which(condition=="excluded")]
            }else{
            tmp_reps <- condition
        }
        
        if(length(tmp_reps) != length(unique(tmp_reps))){
            cds <- estimateDispersions(cds)
            }else{
            cds <- estimateDispersions(cds,method="blind",sharingMode="fit-only")
        }
        
        vsd <- getVarianceStabilizedData(cds)
        
        ## differential gene calculations
        res <- nbinomTest(cds,tmp_conditions[1],tmp_conditions[2])
        colnames(res)[3] <- paste("baseMean",tmp_conditions[1],sep="_")
        colnames(res)[4] <- paste("baseMean",tmp_conditions[2],sep="_")
        # writes to files
        write.table(res,
            file=paste(tmp_comparison,"ALLRESULTS.txt",sep="."),
            quote=F, sep="\t", row.names=F
            )
        
        dir.create("specific_results", showWarnings = FALSE)
        
        write.table(res[which(res$padj < FDR),],
            file=paste("specific_results/",tmp_comparison,".SIGRESULTS.txt",sep=""),
            quote=F, sep="\t", row.names=F
            )
        write.table(res[which(res$foldChange < 1),],
            file=paste("specific_results/",tmp_conditions[1],"(UPREG)_VS_",tmp_conditions[2],"(DOWNREG).ALLRESULTS.txt",sep=""),
            quote=F, sep="\t", row.names=F
            )
        write.table(res[which(res$foldChange > 1),],
            file=paste("specific_results/",tmp_conditions[1],"(DOWNREG)_VS_",tmp_conditions[2],"(UPREG).ALLRESULTS.txt",sep=""),
            quote=F, sep="\t", row.names=F
            )
        write.table(res[which(res$foldChange < 1 & res$padj < FDR),], file=paste("specific_results/",tmp_conditions[1],"(UPREG)_VS_",tmp_conditions[2],"(DOWNREG).SIGRESULTS.txt",sep=""), quote=F, sep="\t", row.names=F)
        write.table(res[which(res$foldChange > 1 & res$padj < FDR),], file=paste("specific_results/",tmp_conditions[1],"(DOWNREG)_VS_",tmp_conditions[2],"(UPREG).SIGRESULTS.txt",sep=""), quote=F, sep="\t", row.names=F)
        
        ## plots of differential genes
        png(file=paste(tmp_comparison,"expressionplot.png",sep="."),height=650,width=650)
        suppressWarnings(plot(
            res$baseMean,
            res$log2FoldChange,
            main=paste(cap(tmp_conditions[1]),"vs.",cap(tmp_conditions[2]),sep=" "),
            ylab=expression("log"[2]*" (Fold Change)"),
            xlab=expression("log"[2]*" (Mean Count Value)"),
            log="x",pch=20,cex=.8,
            col=ifelse(res$padj < FDR,"red","black")
            ))
        legend(
            x="topright",
            "Significantly Differentially Expressed",
            pch=16,
            col=c("red")
            )
        invisible(dev.off())
        
        png(file=paste("specific_results/",tmp_comparison,".expressionplot.png",sep=""),height=650,width=650)
        suppressWarnings(plot(
            res$baseMean,
            res$log2FoldChange,
            main=paste(cap(tmp_conditions[1]),"vs.",cap(tmp_conditions[2]),sep=" "),
            ylab=expression("log"[2]*" (Fold Change)"),
            xlab=expression("log"[2]*" (Mean Count Value)"),
            log="x",pch=20,cex=.8,
            col=ifelse(res$padj < FDR, ifelse(res$foldChange > 1, "red", "green"), "black")
            ))
        legend(
            x="topright",
            c(paste(cap(tmp_conditions[2]),"Up /",cap(tmp_conditions[1]),"Down"),paste(cap(tmp_conditions[1]),"Up /",cap(tmp_conditions[2]),"Down")),
            pch=16,
            col=c("red","green")
            )
        invisible(dev.off())
        
        ## plot of p-values
        png(file=paste(tmp_comparison,"pvalplot.png",sep="."),height=650,width=650)
        hist(
            res$pval,
            breaks=100,
            col="skyblue",border="slateblue",
            main=paste(cap(tmp_conditions[1]),"vs.",cap(tmp_conditions[2]),sep=" "),
            xlab="P-values"
            )
        invisible(dev.off())
        
        ## heatmap of top 30 genes
        ordered_res <- res[order(res$padj),]
        top30_res <- rbind(ordered_res[which(ordered_res$foldChange > 1),][1:15,],ordered_res[which(ordered_res$foldChange < 1),][1:15,])
        top30 <- rowcolours <- top30_res[,1]
        rowcolours[top30_res$padj < FDR] <- "red" ; rowcolours[top30_res$padj >= FDR] <- "black"
        heatcols <- row.names(Design)[condition %in% tmp_conditions]
        
        png(file=paste(tmp_comparison,"top30heatmap.png",sep="."),height=650,width=650)
        margin_extension1 <- max(strwidth(names(sizeFactors(cds)),"inches"))
        if(margin_extension1 < 0){ margin_extension1 <- 0 }
        margin_extension2 <- max(strwidth(top30,"inches")) - 1
        if(margin_extension2 < 0){ margin_extension2 <- 0 }
        par(omi=c(margin_extension1,0.2,0.2,margin_extension2)) # omi to fit labels
        
        heatmap.2(as.matrix(vsd[top30,heatcols]),
            trace="none", dendrogram="column", Rowv=F,
            rowsep=15, sepcolor="black", labRow="",
            add.expr=mtext(side=4, text=rev(top30), at=iy, las=2, line=0.5, col=rev(rowcolours), cex=0.7),
            main=bquote(
                atop(.(cap(tmp_conditions[1])) ~ "vs." ~ .(cap(tmp_conditions[2])),
                    atop("Top 15 Upregulated Genes and Top 15 Downregulated Genes.",
                        atop("Variance Stabilized Data.",
                            "Significant Genes Labelled in Red.")))))
        
        par(omi=c(0,0,0,0)) # reset back to default omi
        invisible(dev.off())
        
        ## changes back to original working directory before next iteration
        setwd(tmp_origwd)
    }
    ''')
    
    # saves workspace and changes back working directory
    R('''
    ## saves workspace for use with EdgeR
    save.image(file="deseq_R_Workspace.RData")
    ## changes back to original python working directory
    setwd(origwd)
    ''')
    
    # restores stdout/stderr
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    
    # ends pipeline log
    endlog('deseq','deseq')


###################################################################

## EdgeR


@follows(deseq)
@split(htseq_merge,['RESULTS/edger/*.txt','RESULTS/edger/*.png','RESULTS/edger/edger_deseq_R_Workspace.RData'])
def edger(infile,outfile):
    
    # makes directory if doesn't already exist
    makedirs("RESULTS/edger")
    
    # starts pipeline log
    startlog('edger','edger')
    
    
    ## EdgeR R script
    
    # redirects stdout and stderr before running any R
    sys.stdout = open("LOGS/edger/edger.log", "a")
    sys.stderr = open("LOGS/edger/edger.log", "a")
    
    # loads libraries and files, and changes working directory
    R('''
    # load libraries and previous DESeq data
    library("edgeR")
    library("DESeq")
    library("gplots")
    load("RESULTS/deseq/deseq_R_Workspace.RData")
    
    # set working directory
    setwd('RESULTS/edger')
    ''')
    
    # runs initial EdgeR analysis
    R('''
    ### Initial analysis
    # input into EdgeR
    group = sapply(strsplit(colnames(countTable),"_",fixed=T),"[",1)
    dge <- DGEList(counts=countTable, group=group)
    
    ## normalization (to effective library size)
    dge2 <- calcNormFactors(dge)
    ''')
    
    # writes intial EdgeR analysis to files
    R('''
    # write size factors to file
    sizefactors <- dge2$samples$norm.factors ##### PRINT LIBRARY SIZES AND NORM FACTORS FROM THIS
    names(sizefactors) <- rownames(dge2$samples)
    write.table(sizefactors, file="normalization_factors.txt", quote=F, sep="\t",col.names=F)
    
    # plot size factors
    png(file="normalization_factors.png",height=650,width=650)
    margin_extension <- max(strwidth(names(sizefactors),"inches")) + 0.3
    par(mai=c(margin_extension,0.82,0.82,0.42)) # mai to fit bottom labels
    
    barplot(sizefactors,main="Normalization Factors",las=2)
    
    par(mai=c(1.02,0.82,0.82,0.42)) # reset back to default mai
    invisible(dev.off())
    
    # write normalized counts to file
    rpm <- cpm(dge2, normalized.lib.sizes=F)
    normalized_rpm <- cpm(dge2, normalized.lib.sizes=F)
    write.table(rpm, file="reads_per_million.txt", quote=F, sep="\t")
    write.table(normalized_rpm, file="normalized_reads_per_million.txt", quote=F, sep="\t")
    ''')
    
    # Multidimensional Scaling (and plot)
    R('''
    ## MDS
    png("/dev/null")
    MDS <- plotMDS(dge2)
    invisible(dev.off())
    
    ## MDS plots
    png(file="Raw_MDS_Dim1Dim2_top500genes.png",height=650,width=650)
    margin_extension <- max(strwidth(as.vector(dge2$samples$group),"inches"))
        
    plot(MDS$x,MDS$y,
        main="Top 500 genes, Multidimensional Scaling Plot",
        xlim=c(min(MDS$x)*1.1,max(MDS$x)*1.1),
        ylim=c(min(MDS$y)*1.1,max(MDS$y)*1.1),
        pch=20,
        xlab="Dimension 1",
        ylab="Dimension 2"
        )
    text(MDS$x, MDS$y, names(MDS$x), pos=3, cex=0.9)
    
    invisible(dev.off())

    png(file="Coloured_MDS_Dim1Dim2_top500genes.png",height=650,width=650)
    cols <- as.vector(dge2$samples$group)
    for(i in 1:length(unique(dge2$samples$group))){
        cols[which(cols==unique(dge2$samples$group)[i])] <- rainbow(length(unique(dge2$samples$group)))[i]
        }
    
    par(omi=c(0.1, 0.1, 0.1, margin_extension + 0.9),xpd=NA) # omi to fit legend
    plot(MDS$x,MDS$y,
        main="Top 500 genes, Multidimensional Scaling Plot",
        pch=19, col=cols,
        xlab="Dimension 1",
        ylab="Dimension 2"
        )
    
    legend(par("usr")[2],
        par("usr")[4],
        as.vector(unique(dge2$samples$group)),
        cex=0.9, pch=19,
        col=rainbow(length(unique(dge2$samples$group)))
        )
    
    par(omi=c(0,0,0,0),xpd=FALSE) # reset back to default omi
    invisible(dev.off())
    ''')
    
    # Performs pairwise Edge R analysis on each pairing (based on expansive level)
    # does both classic approach and GLM approach
    R('''
    ### Pairwise tests
    # creates directories
    dir.create("classic_pairwise", showWarnings = FALSE)
    dir.create("GLM_pairwise", showWarnings = FALSE)
    
    # loop through groupings dependant on expansive level
    for(groupvector in groups){
        
        # determine groups to compare
        tmp_conditions <- as.vector(unique(groupvector))
        if("excluded" %in% tmp_conditions){
            tmp_conditions <- tmp_conditions[-which(tmp_conditions=="excluded")]
            }
        
        design <- model.matrix(~0 + as.factor(groupvector), data=DGEList(counts=countTable,group=groupvector)$samples)
        colnames(design) <- sapply(strsplit(colnames(design),"as.factor(groupvector)",fixed = TRUE),"[[",2)
        
        # creates directories for results
        tmp_comparison <- paste(tmp_conditions[1],"VS",tmp_conditions[2],sep="_")
        dir.create(paste("classic_pairwise",tmp_comparison,sep="/"), showWarnings = FALSE)
        dir.create(paste("GLM_pairwise",tmp_comparison,sep="/"), showWarnings = FALSE)
        dir.create(paste("classic_pairwise",tmp_comparison,"specific_results",sep="/"), showWarnings = FALSE)
        dir.create(paste("GLM_pairwise",tmp_comparison,"specific_results",sep="/"), showWarnings = FALSE)
        
        # input into EdgeR
        dge <- DGEList(counts=countTable, group=groupvector)
        dge2 <- calcNormFactors(dge)
        
        
        ## Estimate dispersions and results
        
        # options changed to skip functions that produce warnings
        old_opts <- options()$warn
        options(warn=2)
        
        # classic dispersions
        disp <- estimateCommonDisp(dge2)
        try(disp <- estimateTagwiseDisp(disp), silent=T)
        
        et <- exactTest(disp,pair=tmp_conditions)
        class_res <- as.data.frame( topTags(et, n=length(rownames(et$table)) ) )
        
        # GLM dispersions
        
        GLMdisp <- estimateGLMCommonDisp(dge2, design)
        try(GLMdisp <- estimateGLMTrendedDisp(GLMdisp, design), silent=T)
        try(GLMdisp <- estimateGLMTagwiseDisp(GLMdisp, design), silent=T)
        
        fit <- glmFit(GLMdisp, design)
        
        contrast <- colnames(design)
        contrast <- gsub(tmp_conditions[1],-1,contrast)
        contrast <- gsub(tmp_conditions[2],1,contrast)
        contrast <- gsub("excluded",0,contrast)
        
        lrt <- glmLRT(dge, fit, contrast=as.numeric(contrast))
        
        GLM_res <- as.data.frame( topTags(lrt, n=length(rownames(lrt$table)) ) )
        
        #####
        ## fix for old version of edger (if newer version, can remove the 2 lines below)
        # (replaces logFC and logCPM of GLM, which is bugged due to inability to cope with zero readcounts)
        GLM_res$logFC <- class_res[match(rownames(GLM_res),rownames(class_res)),]$logFC
        GLM_res$logCPM <- class_res[match(rownames(GLM_res),rownames(class_res)),]$logCPM
        #####
        
        # options changed back
        options(warn=old_opts)
        
        ### function for plots and printing results
        multiplots <- function(wd, res, tmp_comparison, tmp_conditions, dge2, groupvector){
            
            # converts row names to first column
            tmp_res <- cbind(rownames(res), res)
            colnames(tmp_res) <- c("gene_id",colnames(res))
            
            ## writes tables to files
            write.table(tmp_res,
                file=paste(wd,"/",tmp_comparison,"/",tmp_comparison,".ALLRESULTS.txt",sep=""),
                quote=F, sep="\t", row.names=F
                )
            write.table(tmp_res[which(tmp_res$FDR < FDR),],
                file=paste(wd,"/", tmp_comparison, "/", tmp_comparison, ".SIGRESULTS.txt", sep=""),
                quote=F, sep="\t", row.names=F
                )
            write.table(tmp_res[which(tmp_res$logFC < 0),],
                file=paste(wd,"/",tmp_comparison,"/specific_results/",tmp_conditions[1],"(UPREG)_VS_",tmp_conditions[2],"(DOWNREG).ALLRESULTS.txt",sep=""),
                quote=F, sep="\t", row.names=F
                )
            write.table(tmp_res[which(tmp_res$logFC > 0),],
                file=paste(wd,"/",tmp_comparison,"/specific_results/",tmp_conditions[1],"(DOWNREG)_VS_",tmp_conditions[2],"(UPREG).ALLRESULTS.txt",sep=""),
                quote=F, sep="\t", row.names=F
                )
            write.table(tmp_res[which(tmp_res$logFC < 0 & tmp_res$FDR < FDR),],
                file=paste(wd,"/",tmp_comparison,"/specific_results/",tmp_conditions[1],"(UPREG)_VS_",tmp_conditions[2],"(DOWNREG).SIGRESULTS.txt",sep=""),
                quote=F, sep="\t", row.names=F
                )
            write.table(tmp_res[which(tmp_res$logFC > 0 & tmp_res$FDR < FDR),],
                file=paste(wd,"/",tmp_comparison,"/specific_results/",tmp_conditions[1],"(DOWNREG)_VS_",tmp_conditions[2],"(UPREG).SIGRESULTS.txt",sep=""),
                quote=F, sep="\t", row.names=F
                )
            
            ## plots of differential genes
            png(file=paste(wd,"/",tmp_comparison,"/",tmp_comparison,".expressionplot.png",sep=""),height=650,width=650)
            suppressWarnings(plot(
                res$logCPM,
                res$logFC,
                main=paste(cap(tmp_conditions[1]),"vs.",cap(tmp_conditions[2]),sep=" "),
                ylab=expression("log"[2]*" (Fold Change)"),
                xlab=expression("Mean log"[2]*" (Counts-per-Million)"),
                pch=20,cex=.8,
                col=ifelse(res$FDR < FDR,"red","black")
                ))
            legend(
                x="topright",
                "Significantly Differentially Expressed",
                pch=16,
                col=c("red")
                )
            invisible(dev.off())
            
            png(file=paste(wd,"/",tmp_comparison,"/specific_results/",tmp_comparison,".expressionplot.png",sep=""),height=650,width=650)
            suppressWarnings(plot(
                res$logCPM,
                res$logFC,
                main=paste(cap(tmp_conditions[1]),"vs.",cap(tmp_conditions[2]),sep=" "),
                ylab=expression("log"[2]*" (Fold Change)"),
                xlab=expression("Mean log"[2]*" (Counts-per-Million)"),
                pch=20,cex=.8,
                col=ifelse(res$FDR < FDR, ifelse(res$logFC > 0, "red", "green"), "black")
                ))
            legend(
                x="topright",
                c(paste(cap(tmp_conditions[2]),"Up /",cap(tmp_conditions[1]),"Down"),paste(cap(tmp_conditions[1]),"Up /",cap(tmp_conditions[2]),"Down")),
                pch=16,
                col=c("red","green")
                )
            invisible(dev.off())
            
            ## plot of p-values
            png(file=paste(wd,"/",tmp_comparison,"/",tmp_comparison,".pvalplot.png",sep=""),height=650,width=650)
            hist(
                res$PValue,
                breaks=100,
                col="skyblue",border="slateblue",
                main=paste(cap(tmp_conditions[1]),"vs.",cap(tmp_conditions[2]),sep=" "),
                xlab="P-values"
                )
            invisible(dev.off())
            
            ## heatmap of top 30 genes
            top30_res <- rbind(res[which(res$logFC > 0),][1:15,], res[which(res$logFC < 0),][1:15,])
            top30 <- rowcolours <- rownames(top30_res)
            rowcolours[top30_res$FDR < FDR] <- "red" ; rowcolours[top30_res$FDR >= FDR] <- "black"
            heatcols <- rownames(dge2$samples)[groupvector %in% tmp_conditions]
            
            png(file=paste(wd,"/",tmp_comparison,"/",tmp_comparison,".top30heatmap.png",sep=""),height=650,width=650)
            margin_extension1 <- max(strwidth(rownames(dge2$samples),"inches"))
            if(margin_extension1 < 0){ margin_extension1 <- 0 }
            margin_extension2 <- max(strwidth(top30,"inches")) - 1
            if(margin_extension2 < 0){ margin_extension2 <- 0 }
            par(omi=c(margin_extension1,0.2,0.2,margin_extension2)) # omi to fit labels
            
            heatmap.2(as.matrix(cpm(dge2,normalized.lib.sizes=T)[top30,heatcols]),
                scale="row", trace="none", dendrogram="column", Rowv=F,
                rowsep=15, sepcolor="black", labRow="",
                add.expr=mtext(side=4, text=rev(top30), at=iy, las=2, line=0.5, col=rev(rowcolours), cex=0.7),
                main=bquote(
                    atop(.(cap(tmp_conditions[1])) ~ "vs." ~ .(cap(tmp_conditions[2])),
                        atop("Top 15 Upregulated Genes and Top 15 Downregulated Genes.",
                            atop("Data Scaled by Row.",
                                "Significant Genes Labelled in Red.")))))
            
            par(omi=c(0,0,0,0)) # reset back to default omi
            invisible(dev.off())
        }
        
        multiplots("classic_pairwise", class_res, tmp_comparison, tmp_conditions, dge2, groupvector)
        multiplots("GLM_pairwise", GLM_res, tmp_comparison, tmp_conditions, dge2, groupvector)
        
    }
    ''')
    
    ## performs ANOVA-like GLM test
    # reads in data required and does base calculations
    R('''    
    ### ANOVA-like test
    # creates directories
    dir.create("ANOVA_like", showWarnings = FALSE)
    
    # reads in countdata
    dge <- DGEList(counts=countTable, group=group)
    dge2 <- calcNormFactors(dge)
    
    # creates design matrix
    design <- model.matrix(~as.factor(group), data=dge$samples)
    ''')
    
    # calculates dispersions
    R('''
    ## GLM dispersions
    # options changed to skip functions that produce warnings
    old_opts <- options()$warn
    options(warn=2)
    
    GLMdisp <- estimateGLMCommonDisp(dge2, design)
    try(GLMdisp <- estimateGLMTrendedDisp(GLMdisp, design), silent=T)
    try(GLMdisp <- estimateGLMTagwiseDisp(GLMdisp, design), silent=T)
    
    fit <- glmFit(GLMdisp, design)
    
    lrt <- glmLRT(dge, fit, coef=2:length(colnames(design)))
    
    GLM_res <- as.data.frame( topTags(lrt, n=length(rownames(lrt$table)) ) )
    
    # options changed back
    options(warn=old_opts)
    ''')
    
    # writes results to files
    R('''
    ## Plots and tables
    # tables of results
    write.table(GLM_res, file="ANOVA_like/ANOVA_GLM.ALLRESULTS.txt", quote=F, sep="\t", row.names=F)
    write.table(GLM_res[which(GLM_res$FDR < FDR),], file="ANOVA_like/ANOVA_GLM.SIGRESULTS.txt", quote=F, sep="\t", row.names=F)
    
    # plot of p-values
    png(file="ANOVA_like/ANOVA_GLM.pvalplot.png",height=650,width=650)
    hist(
        GLM_res$PValue,
        breaks=100,
        col="skyblue",border="slateblue",
        main="ANOVA-like GLM Testing",
        xlab="P-values"
        )
    invisible(dev.off())
    
    # heatmap of top 30 genes
    top30_res <- GLM_res[1:30,]
    top30 <- rowcolours <- rownames(top30_res)
    rowcolours[top30_res$FDR < FDR] <- "red" ; rowcolours[top30_res$FDR >= FDR] <- "black"
    heatcols <- rownames(dge2$samples)
    
    png(file="ANOVA_like/ANOVA_GLM.top30heatmap.png",height=650,width=650)
    margin_extension1 <- max(strwidth(rownames(dge2$samples),"inches"))
    if(margin_extension1 < 0){ margin_extension1 <- 0 }
    margin_extension2 <- max(strwidth(top30,"inches")) - 1
    if(margin_extension2 < 0){ margin_extension2 <- 0 }
    par(omi=c(margin_extension1,0.2,0.2,margin_extension2)) # omi to fit labels
    
    heatmap.2(as.matrix(cpm(dge2,normalized.lib.sizes=T)[top30,heatcols]),
        scale="row", trace="none", dendrogram="column", Rowv=F,
        rowsep=15, sepcolor="black", labRow="",
        add.expr=mtext(side=4, text=rev(top30), at=iy, las=2, line=0.5, col=rev(rowcolours), cex=0.7),
        main=bquote(
            atop("ANOVA-like GLM Testing.",
                atop("Top 30 Differentially Expressed Genes.",
                    atop("Data Scaled by Row.",
                        "Significant Genes Labelled in Red.")))))
        
        
        par(omi=c(0,0,0,0)) # reset back to default omi
        invisible(dev.off())
    ''')
    
    # saves workspace and changes back working directory
    R('''
    ## saves workspace for use with EdgeR
    save.image(file="edger_deseq_R_Workspace.RData")
    ## changes back to original python working directory
    setwd(origwd)
    ''')
    
    # restores stdout/stderr
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    
    # ends pipeline log
    endlog('edger','edger')
    



###################################################################
###################################################################
###################################################################



###########################
###    Running lines    ###
###########################

## obtains instructions from command line, raises errors if not provided

if len(sys.argv) < 2: sys.exit("\nError:\n\tNo task name provided\n\tTo run all tasks in pipeline use 'all'\n\t\tE.G. python harryc_DE_pipeline.py all 10\n\tTo produce a graph of pipeline and task names use 'graph'\n\t\tE.G. python harryc_DE_pipeline.py graph\n")

if sys.argv[1] == "all": run = edger
else: run = sys.argv[1]

if len(sys.argv) >= 3:
    try: par = int(sys.argv[2])
    except ValueError: sys.exit("\nError:\n\tIf second argument provided, must be integer for number of parallels\n\t\tE.G. python harryc_DE_pipeline.py all 10\n")
else:
    par = 1

## runs pipeline

if run == "graph":
    
    pipeline_printout_graph("pipeline_printout.svg", "svg", [edger], no_key_legend=False, minimal_key_legend=True, user_colour_scheme={"colour_scheme_index":6}, draw_vertically=True, pipeline_name='Differential Expression Pipeline by Harry Clifford')
    
else:
    
    makedirs("LOGS")
    
    pipeline_printout_graph("LOGS/pipeline_printout.svg", "svg", [run], no_key_legend = False, minimal_key_legend = True, user_colour_scheme = {"colour_scheme_index":6}, draw_vertically=True, pipeline_name='Differential Expression Pipeline by Harry Clifford')
    
    pipeline_run([run], multiprocess=int(par), verbose = 4)
    
    pipeline_printout_graph("LOGS/pipeline_printout.svg", "svg", [run], no_key_legend = False, minimal_key_legend = True, user_colour_scheme = {"colour_scheme_index":6}, draw_vertically=True, pipeline_name='Differential Expression Pipeline by Harry Clifford')




###################################################################
###################################################################
###################################################################









##### TO ADD

## RPKM values following edger
# (take the logCPM values and subtract the log2 of the gene length)
# http://www.biostars.org/p/61523/
# https://groups.google.com/forum/?fromgroups=#!topic/rsem-users/iSYkFWaKxz4
# or
# http://cufflinks.cbcb.umd.edu/manual.html#fpkm_track

#####



