############################################################################################
#  Script to plot number of reads in a specified region (normalized to reads per million)  #
############################################################################################


##############
# User Input #
##############


# working directory
import os
os.chdir('/net/isi-scratch/normalized_peak_plot')


# region to be plotted
# (standard format chr:start-stop)
# (can extend this a little either side of the plot if you would like to see whole peak shape)
peak = '1:42739500-42763500'


# bam files to plot peaks for
# (names in a vector)
infiles = ['inputs/N2A-BRG1_KD-rep1.bam','inputs/N2A-BRG1_KD-rep2.bam','inputs/N2A-BRG1_WT-rep1.bam','inputs/N2A-BRG1_WT-rep2.bam','inputs/N2A-INPUT_KD-rep1.bam','inputs/N2A-INPUT_WT-rep1.bam']


# coveragebed command
# (the command you wish to use to run coveragebed - as this can take some time, may want to add cluster line)
coveragebed = "qrsh -cwd -q medium_jobs.q -v BASH_ENV='~/.bashrc' -now n coverageBed"


# maximum number of parallels to run at any one time
parallels = 6



############################################################################################

### Main script

from rpy2.robjects import r as R
from multiprocessing import Pool
import subprocess



## creates bed file for peak

ranges = range(int(peak.split(":")[-1].split("-")[0]),int(peak.split(":")[-1].split("-")[1]))
full_ranges = ['\t'.join([peak.split(":")[0],str(i),str(i+1)]) for i in ranges]
bed = '\n'.join(full_ranges)

with open('TMP_peakfile.bed', 'w') as tmp: tmp.write(bed)



## runs coveragebed on bed file

# function to run coveragebed
def fun(infile):
    coveragebed_cmd = coveragebed + ' -abam ' + infile + ' -b TMP_peakfile.bed > TMP_' + infile.rsplit('/')[-1].rsplit('.bam')[0] + '.bed'
    os.system(coveragebed_cmd)

# runs in parallel
pool = Pool(processes=int(parallels))
pool.map(fun, infiles)



## calculates total number of reads mapped in each bam
mapped_reads = []
for infile in infiles:
    tmp_reads = os.popen('samtools view -c -F 4 ' + infile).read()
    bed_filename = 'TMP_' + infile.rsplit('/')[-1].rsplit('.bam')[0] + '.bed'
    mapped_reads.append([bed_filename,tmp_reads])

mapfile = "".join(["\t".join(i) for i in mapped_reads])
with open('TMP_mapfile.txt', 'w') as tmp: tmp.write(mapfile)



## runs R on tmp files to produce graphs

# calculates reads per million for each file
R('''

r_mapfile <- read.delim('TMP_mapfile.txt',header=FALSE,stringsAsFactors=FALSE)
to_plot <- list()
for(i in 1:length(r_mapfile[,1])){
    infile <- r_mapfile[i,1]
    million_reads <- as.integer(r_mapfile[i,2])/1000000
    rpm <- read.delim(infile,header=FALSE,stringsAsFactors=FALSE)[,4] / million_reads
    to_plot[[infile]] <- rpm
}

''')

# draw initial plot
R('''

png(filename='peak_plot.png',width=1200,height=600)
bases <- read.delim(r_mapfile[1,1],header=FALSE,stringsAsFactors=FALSE)[,2]
maximum <- max(unlist(to_plot))
col_pal <- colorRampPalette(c("red","red","green","green","blue","blue"))
colrs <- col_pal(length(r_mapfile[,1]))
initial_points <- to_plot[[r_mapfile[1,1]]]
plot(cbind(bases,initial_points), type="l", col=colrs[1], ylim=c(0,maximum), xlab="Base Position", ylab="Reads-per-million", lwd=3, lty=1)

''')

# adds line for each file
R('''

for(i in 2:length(r_mapfile[,1])){
    list_component <- r_mapfile[i,1]
    if(i%%2==0){line_type <- 3}else{line_type <- 1}
    lines(cbind(bases,to_plot[[list_component]]), type="l", col=colrs[i], lwd=3, lty=line_type)
}

''')

# adds legend and closes plot
R('''

legend('topright', gsub('TMP_|.bed','',r_mapfile[,1]), col=colrs, lty=1:2)
dev.off()

''')



## deletes tmp files

for infile in infiles:
    tmp_infile = 'TMP_'+infile.rsplit('/')[-1].rsplit('.bam')[0]+'.bed'
    os.remove(tmp_infile)

os.remove('TMP_mapfile.txt')
os.remove('TMP_peakfile.bed')











