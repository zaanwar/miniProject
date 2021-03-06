##MiniProject
import os
from Bio import Entrez
from Bio import SeqIO
import shutil #use to move around test files
from Bio.Blast import NCBIWWW

#create all necessary files, directory will be used to store output
os.system('mkdir miniProject_Zain_Anwar')#make directory
os.chdir('miniProject_Zain_Anwar')#movrm e into directory
os.system('touch miniProject.log')#make log file inside directory

######Step 1######

userInput = input("What set of data would you like to use? Please enter: test OR full\n")#ask user which input data they will use

#if user selects full dataset, retrieve the data from ncbi and coverted into paired-end fastq files
if userInput == "full":
    #Retrieve the transcriptomes from patient donors using the wget command#
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1')#Donor 1(2dpi)
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1')#Donor 1(6dpi)
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1')#Donor 3(2dpi)
    os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1')#Donor 3(6dpi)
    #Convert the trasncriptome files to paired-end fastq files
    os.system('fastq-dump -I --split-files SRR5660030.1')
    os.system('fastq-dump -I --split-files SRR5660033.1')
    os.system('fastq-dump -I --split-files SRR5660044.1')
    os.system('fastq-dump -I --split-files SRR5660045.1')

#if user selects test datast, use the shortened test data provided on github, renameand move the test data so the will work in place of the paired-end fast q files
elif userInput == "test":
    os.chdir('..') 
    shutil.copyfile('test_files/SRR5660030.1_1_test.fastq', 'miniProject_Zain_Anwar/SRR5660030.1_1.fastq')
    shutil.copyfile('test_files/SRR5660030.1_2_test.fastq', 'miniProject_Zain_Anwar/SRR5660030.1_2.fastq')
    shutil.copyfile('test_files/SRR5660033.1_1_test.fastq', 'miniProject_Zain_Anwar/SRR5660033.1_1.fastq')
    shutil.copyfile('test_files/SRR5660033.1_2_test.fastq', 'miniProject_Zain_Anwar/SRR5660033.1_2.fastq')
    shutil.copyfile('test_files/SRR5660044.1_1_test.fastq', 'miniProject_Zain_Anwar/SRR5660044.1_1.fastq')
    shutil.copyfile('test_files/SRR5660044.1_2_test.fastq', 'miniProject_Zain_Anwar/SRR5660044.1_2.fastq')
    shutil.copyfile('test_files/SRR5660045.1_1_test.fastq', 'miniProject_Zain_Anwar/SRR5660045.1_1.fastq')
    shutil.copyfile('test_files/SRR5660045.1_2_test.fastq', 'miniProject_Zain_Anwar/SRR5660045.1_2.fastq')
    os.chdir('miniProject_Zain_Anwar')

else:
    print("please enter: 'test' or 'full' ")
    
######Step 2######

#Build a transcriptome index for HCMV (NCBI accession EF999921). We will use biopython to retirevce and generate the appropiate input and then build the
#the index with kallisto. Useb biopython to retireve the HCMV genome from NCBI
Entrez.email = "zanwar@luc.edu"
handle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype = 'gb', retmode = 'text')
sequences = SeqIO.read(handle, 'genbank')

cds_FASTA = open('HCMVtranscriptome.fasta', 'w')#file to hold CDS's
outfile = open('miniProject.log','w')#log file

count = 0#will be used to track the number of CDS's

#for loop checks feature of every entry, if it is a CD, save the protein name and sequence in FASTA format in the CDS file
for sequence  in sequences.features:
    if sequence.type == 'CDS':
        count = count + 1 
        cds_FASTA.write('>' + sequence.qualifiers['protein_id'][0]  + '\n' +  str(sequence.extract(sequences.seq)) + '\n')
outfile.write('The HCMV genome (EF9999921) has' + str(count) + 'CDS.')
outfile.close()
cds_FASTA.close()
os.system('time kallisto index -i HCMV_index.idx HCMVtranscriptome.fasta')#build  index

######Step 3######

#Quantify the TPM of each CDS through kallisto. Then take the results and use them as an input to find the differentially expressed genes between the two timepoints (2 and 6dpi) for both patients.

#Quantify the TPM of each CDS using kallisto. 
os.system('mkdir TPM_results')#will hold results of the kallisto
os.system('time kallisto quant -i HCMV_index.idx -b 30 -t 2 -o TPM_results/SRR5660030.1 SRR5660030.1_1.fastq SRR5660030.1_2.fastq')
os.system('time kallisto quant -i HCMV_index.idx -b 30 -t 2 -o TPM_results/SRR5660033.1 SRR5660033.1_1.fastq SRR5660033.1_2.fastq')
os.system('time kallisto quant -i HCMV_index.idx -b 30 -t 2 -o TPM_results/SRR5660044.1 SRR5660044.1_1.fastq SRR5660044.1_2.fastq')
os.system('time kallisto quant -i HCMV_index.idx -b 30 -t 2 -o TPM_results/SRR5660045.1 SRR5660045.1_1.fastq SRR5660045.1_2.fastq')

#make table that will be used for sleuth
table = open('table.txt', 'w')
table.write('sample condition path\n')
table.write('SRR5660030.1 2dpi TPM_results/SRR5660030.1\n')
table.write('SRR5660033.1 6dpi TPM_results/SRR5660033.1\n')
table.write('SRR5660044.1 2dpi TPM_results/SRR5660044.1\n')
table.write('SRR5660045.1 6dpi TPM_results/SRR5660045.1\n')
table.close()

#R package sleuth to find differentially expressesd genes between two (2 and 6dpi) timepoints
R_script = open('R_script.R','w')
            
R_script.write('library(sleuth)\n')
R_script.write('table <- read.table("table.txt", header=TRUE, stringsAsFactors=FALSE)\n' +'so <- sleuth_prep(table)\n' + 'so <- sleuth_fit(so, ~condition, "full")\n' +'so <- sleuth_fit(so, ~1, "reduced")\n' + 'so <- sleuth_lrt(so, "reduced", "full")\n')
R_script.write('library(dplyr)\n')
R_script.write('s_table <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)\n' + 'sleuth_significant <- dplyr::filter(s_table, qval <= 0.05) %>% dplyr::arrange(pval)\n' +'write.table(sleuth_significant, file="fdr05_results.txt", quote=FALSE, row.names=FALSE)\n')                    
R_script.close()

os.system('Rscript R_script.R')# run R script:

# read in the output from sleuth, take the necessary data and add to the log file
output = open('fdr05_results.txt', 'r')
ids = [] 
for line in output.readlines():
    ids.append(line.split(' '))

outfile = open('miniProject.log', 'a')

outfile.write(ids[0][0] +'\t'+ ids[0][3] +'\t'+ ids[0][1] +'\t'+ ids[0][2] +'\n')#headers

#line up table appropiately
for i in range(1, len(ids)):
    outfile.write(ids[i][0] + '\t' + ids[i][3] + '\t' + ids[i][1] + '\t' + ids[i][2] + '\n')
outfile.write('\n')
outfile.close()
output.close()

######Step 4######

#assemble transcriptome reads so they can be used in BLAST. Before that, we are going to create an index for HCMV using Bowtie 2. Then we will save the reads that map to ONLY the  HCMV index. 

# use entrez efetch to get the full HCMV genome sequence to make bowtie index
Complete_sequence = open('HCMV_index.fasta', 'w') 

handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='fasta', retmode='text')
sequence = SeqIO.read(handle, 'fasta')

Complete_sequence.write('>' + sequence.description + '\n' + str(sequence.seq))
Complete_sequence.close()

# make the bowtie index using the previously created fasta file (HCMV_complete_seq.fasta)
os.system('bowtie2-build HCMV_index.fasta HCMV_genome')


#run bowtie2 mapping
#--al-conc will output filtered reads
os.system( 'bowtie2 --quiet -x HCMV_genome -1 SRR5660030.1_1.fastq -2 SRR5660030.1_2.fastq -S SRR5660030.1_map.sam --al-conc SRR5660030.1_mapped.fastq')
os.system( 'bowtie2 --quiet -x HCMV_genome -1 SRR5660033.1_1.fastq -2 SRR5660033.1_2.fastq -S SRR5660033.1_map.sam --al-conc SRR5660033.1_mapped.fastq')
os.system( 'bowtie2 --quiet -x HCMV_genome -1 SRR5660044.1_1.fastq -2 SRR5660044.1_2.fastq -S SRR5660044.1_map.sam --al-conc SRR5660044.1_mapped.fastq')
os.system( 'bowtie2 --quiet -x HCMV_genome -1 SRR5660045.1_1.fastq -2 SRR5660045.1_2.fastq -S SRR5660045.1_map.sam --al-conc SRR5660045.1_mapped.fastq')


outfile = open('miniProject.log', 'a')
Day_2_files = ['SRR5660030.1_1.fastq', 'SRR5660033.1_1.fastq', 'SRR5660044.1_1.fastq', 'SRR5660045.1_1.fastq']
Day_6_files = ['SRR5660030.1_mapped.1.fastq', 'SRR5660033.1_mapped.1.fastq', 'SRR5660044.1_mapped.1.fastq', 'SRR5660045.1_mapped.1.fastq']

donor_ids = ['Donor 1 (2dpi)', 'Donor 1 (6dpi)', 'Donor 3 (2dpi)', 'Donor 3 (6dpi)']

#for loop iterates through day 2 and day 6 transcriptomes and adds the number of reads for each before and after to log file
for i in range(len(Day_2_files)):
    #count to keep track of number of reads
    Day_2_count = 0
    Day_6_count = 0
    Day_2_file = open(Day_2_files[i], 'r')
    Day_6_file = open(Day_6_files[i], 'r')
    #iterate through both day 2 and day 6 reads
    for record in SeqIO.parse(Day_2_file, 'fastq'):
        Day_2_count =Day_2_count + 1
    for record in SeqIO.parse(Day_6_file, 'fastq'):
        Day_6_count = Day_6_count + 1
    outfile.write(donor_ids[i] + ' had ' + str(Day_2_count) + ' read pairs before Bowtie2 filtering and ' + str(Day_6_count) +
                   ' read pairs after.\n')


######Step 5######
#Assemble all 4 transcriptomes together using SPAdes
command = ('spades -k 127 -t 2 --only-assembler ' + '--pe1-1 SRR5660030.1_mapped.1.fastq --pe1-2 SRR5660030.1_mapped.2.fastq ' +'--pe2-1 SRR5660033.1_mapped.1.fastq --pe2-2 SRR5660033.1_mapped.2.fastq ' +'--pe3-1 SRR5660044.1_mapped.1.fastq --pe3-2 SRR5660044.1_mapped.2.fastq ' + '--pe4-1 SRR5660045.1_mapped.1.fastq --pe4-2 SRR5660045.1_mapped.2.fastq ' +'-o final_assembly/')
os.system(command)
outfile.write(command)


######Step 6######
#calculate they number of contigs that have a lenght < 1000, write this value to the log_file
contig = open('final_assembly/contigs.fasta', 'r')


count_contigs = 0 #count to number of contigs with length >1000
long_contigs = []#list to hold the longest contigs
for record in SeqIO.parse(contig, "fasta"):
    if len(record.seq) > 1000: #checks length of sequence
        count_contigs = count_contigs + 1
        long_contigs.append(str(record.seq)) #add contig to list
outfile.write('\nThere are ' + str(count_contigs) + ' contigs > 1000 bp in the assembly.\n')

####Step 7######
#Calculate length of assembly
length = 0
for record in SeqIO.parse(contig, "fasta"):
    if len(record.seq) > 1000: #checks length of sequence
              length = length +  len(record.seq) #adds to sum total of bp
outfile.write('\nThere are ' + str(length) + ' bp in the assembly.\n')
              
######Step 8######
#retireve longest contig from assembly and use it as blast+ input. Identify top 10 hits and write to log file

sorted_list = sorted(long_contigs, key = len, reverse = True) #sort the sequence list by length


Input_file = open('blast_input.fasta', 'w')
Input_file.write(str(sorted_list[0])) #store longest sequence
Input_file.close()
            
handle = Entrez.esearch(db = "nucleotide", term = "txid10357[Organism:exp]", retmax = 20000) #search for refseq betaherps seqs
entry = Entrez.read(handle)
ids = entry.get('IdList') #saves ids of the search results

#retrieve fasta sequences for those ids
handle = Entrez.efetch(db = "nucleotide", id = ids, rettype = 'fasta')
entries = list(SeqIO.parse(handle, 'fasta'))

# seqs in fasta format for database
fasta_sequences = open('fasta_sequences.fasta', 'w')
for entry in entries:
    fasta_sequences.write(entry.format('fasta'))
fasta_sequences.close()

#create local database
os.system("makeblastdb -in fasta_sequences.fasta -out HCMVblast -title HCMVblast -dbtype nucl") 
              
# blastn command
os.system('blastn -query blast_input.fasta -db HCMV_blastdb -out blastn_results.txt -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"')

#identify the top 10 hits and write them to the log file
results = open('blastn_results.txt', 'r')
              
#write headers given in step 8 directions
outfile.write('sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')

#write out top 10 reads to log file
for i in range(10): 
    outfile.write(str(results.readline()))
outfile.close()




              
              
