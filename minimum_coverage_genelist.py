
#####################################################################################################
# Script to take coverage files (following coveragebed) and give a list of genes that have at least one transcript with an coverage of a user-given percentage over all samples
#####################################################################################################


###################
# user defined info
###################

# working directory must contain the .coverage files, which are gtf files run through coveragebed each with a bam for each sample
#(each sample must have been run on the same gtf file)
# output will be minimum_coverage_genelist.txt

# coverage limit (between 0 and 1)
cov_lim = 0.8


########################################



# imports required packages
import os
import glob

# find infiles
infiles = glob.glob('*.coverage')

### loops through files creating full dictionary of all coverages
fulldic = {}
for infile in infiles:
    
    # read in data
    env = open(infile,"r")
    data = env.read().split("\n")
    del data[-1]
    dic = {}
    
    ## first loop through line by line - build dictionary
    for line in data:
        # split line up and extract transcript
        transcript = line.split('transcript_id "')[-1].split('";')[0]
        # add data to dictionary of transcripts - if already exists, adds values together
        splitline = line.split('\t')
        if transcript not in dic: dic[transcript] = [splitline[-3],splitline[-2]]
        if transcript in dic: dic[transcript] = [int(dic[transcript][0])+int(splitline[-3]), int(dic[transcript][1])+int(splitline[-2])]
    
    ## second loop through line by line - calculate coverage at each transcript
    for key in dic.keys(): dic[key] = float(dic[key][0])/float(dic[key][1])
    
    # adds dic to fulldic
    fulldic[infile] = dic
    # writes dic to outside file
    cov_output = ''
    for key in dic.keys(): cov_output = cov_output + key + '\t' + str(dic[key]) + '\n'
    with open('.genecoverage'.join(infile.rsplit('.coverage',1)),'w') as final: final.write(cov_output)


### uses gtf of the last file through the previous loop to make dictionary of genes:transcripts
genedic = {}
for line in data:
    
    gene = line.split('gene_id "')[-1].split('";')[0]
    transcript = line.split('transcript_id "')[-1].split('";')[0]
    
    # add data to dictionary if doesn't already exist
    if gene not in genedic: genedic[gene] = [transcript]
    if gene in genedic:
        if transcript not in genedic[gene]: genedic[gene] = genedic[gene] + [transcript]


### loops through dictionary of genes deciding whether to include in final list
final_list = []
for gene in genedic:
    
    # assume gene will not be kept
    keep = False
    
    ## loops through each sample
    for sample in fulldic.keys():
        
        # stops looking if already found a coverage above threshold
        if keep == True: break
        
        ## loops through each transcript
        for transcript in genedic[gene]:
            
            # stops looking if already found a coverage above threshold
            if keep == True: break
            
            # changes keep to true if finds a coverage above threshold
            if fulldic[sample][transcript] >= float(cov_lim): keep = True
    
    # if keep has now been changes to false, adds gene to final list
    if keep: final_list.append(gene)


## writes results to file

line1 = [ 'Number of genes before: ' + str( len( genedic.keys() ) ) ]
line2 = [ 'Number of genes after: ' + str( len( final_list ) ) ]
line3 = [ 'Total genes removed due to low coverage: ' + str( len( genedic.keys() ) - len( final_list ) ) ]
line4 = [ "\nRemaining genes:" ]

final_output = "\n".join(line1 + line2 + line3 + line4 + final_list)
with open('minimum_coverage_genelist.txt','w') as final: final.write(final_output)




