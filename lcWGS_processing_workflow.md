S. orbicularis low-coverage whole genome sequencing processing
================
Bailey Carlson

#### Quick notes/code to remember

``` bash
# Sample M0D039945U only has sequences from Lane 1 and the QC report for that sample is low quality
# M0D040049U_S55_L003 has issues

# Login to cluster
ssh bcarlson4@login.rc.ucmerced.edu

# Determine the size and permissions of each file in a folder
ls -lh

# Determine the number of files present in a folder
num_files=$(ls -l | grep -v ^d | wc -l)
echo "Number of files in the directory: $num_files"

# Copy folders
scp -r

# File permissions *NEED to do this for all shell scripts you create on the cluster!*
chmod +x /path/file_name.sh

# Determine total storage of a folder and folders therein
du -sh .

# Moved tar files onto cluster using this command
scp fastq_23203Daw_N23157_L00*.tar bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/raw_sequences/L00*

# Check format of bam file
module load samtools/1.13
samtools view -H -X M0D039797C_S102_clipped_realigned.bam | head

# Check your bam file
module load samtools/1.13
samtools quickcheck -qvvv your.bam

# Check # of lines in a file
wc -l

# Shows the resources the job used
scontrol show job JOBID

# Create a chromosome specific output file from a bam file
samtools view -bo out.chr21.bam in.bam chr21
samtools view out.chr21.bam | head -n 5

# Install R packages on the cluster
module load ananconda
source activate my-R
conda install r-package

# Check job usage
seff <job_number>

# Run a for loop on chromosome files for line counts
for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24 mtDNA; do
wc -l S_orbicularis_gl_$chr.snps
done
```

#### Before working with the reference nuclear genome (GCF_902148855.1) all single digit chromosome numbers need to be changed to have 0 before the number (i.e. chr01). Additionally, you have to add in the complete mtDNA genome (AP018927.1) to the reference nuclear genome fasta file.

## 1.reference

#### Index your reference genome

- Takes less than 7 minutes
- bwa-mem2 v2.2.1, gatk v4.2.6, samtools v1.13
- Removed NCBI chromosome names and all scaffold sequences from
  GCF_902148855.1_fSphaOr1.1_genomic.fna to create ref.fna
- File ref.fna has chr# and all spaces filled with an underscore
- Indexing requires 28N GB memory where N is the size of the reference
  sequence
- bwa-mem manual <https://bio-bwa.sourceforge.net/bwa.shtml>

``` bash
#!/bin/bash

#SBATCH -J index_reference.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 2G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_index_reference_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_index_reference_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load bwa-mem2 module
module load bwa-mem2/2.2.1

# Load gatk module
module load gatk/4.2.6

# Load samtools module
module load samtools/1.13

# Assuming your files are in the current directory
ref_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/"

      # Indexing the reference sequence
        bwa-mem2 index \
          -p ${ref_dir}ref \
          ${ref_dir}ref.fna

      gatk CreateSequenceDictionary \
          -R ${ref_dir}ref.fna

        samtools faidx \
          ${ref_dir}ref.fna
```

## 2.trim_align_sequences

### 1.load_raw_sequences

#### Trim file names

- Takes less than 10 minutes

``` bash
#!/bin/bash

#SBATCH -J load_raw_sequences.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_load_raw_sequences_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_load_raw_sequences_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/1.raw_sequences/"

tar -xf "$input_dir"/L001/fastq_23203Daw_N23157_L001.tar -C "$input_dir"/L001/
tar -xf "$input_dir"/L002/fastq_23203Daw_N23157_L002.tar -C "$input_dir"/L002/
tar -xf "$input_dir"/L003/fastq_23203Daw_N23157_L003.tar -C "$input_dir"/L003/
tar -xf "$input_dir"/L004/fastq_23203Daw_N23157_L004.tar -C "$input_dir"/L004/
```

#### Use mv 23203Daw_R23157_M0D0\* ../. to move all the sequence files up into the main lane folder L00\*

## 2.trim_align_sequences

### 2.rename_sequences

#### Trim file names

- Takes less than 10 minutes
- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J rename_sequences_L001.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH --ntasks=4 #number of cores per node
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_rename_sequences_L001_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_rename_sequences_L001_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/1.raw_sequences/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/2.rename_sequences/L00*/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for file_R1 in "${input_dir}"*_R1_001.fastq.gz; do
    # Generate the corresponding R2 file name
    file_R2="${file_R1/_R1_/_R2_}"

    # Extract the sample name from the file name
    sample_name=$(echo "$file_R1" | cut -d '_' -f 3-6)
    
    # Create new output files with the desired name structure
    cp "$file_R1" "${output_dir}${sample_name}.fastq.gz"
    cp "$file_R2" "${output_dir}${sample_name}.fastq.gz" 
done
```

## 2.trim_align_sequences

### 3.trimmomatic_trim_filter

#### Trim sequences according to quality and adapters

- Takes 16.5 hours
- trimmomatic v0.39
- trimmomatic manual
  <https://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf>
- The sequences are paired end so we use PE
- Determine phred quality score using nomenclature from raw fastq files
  but trimmomatic can determine automatically zcat sample.fastq.gz \|
  head -n 8 <https://drive5.com/usearch/manual/quality_score.html>
- trimmomatic PE \[- threads <threads>\] \[-phred33 \| -phred64\]
  \[-basein <inputBase> \| \<input 1\> \<input 2\>\] \[-baseout
  <outputBase> \| \<paired output 1\> \<unpaired output 1\> \<paired
  output 2\> \<unpaired output 2\>
  ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
  <minimum length of sequence>
- I used the TruSeq3-PE-2.fa for adapters found on the
  <https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE-2.fa>
- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J trimmomatic_trim_filter_L00*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 0-16:30:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_trimmomatic_trim_filter_L00*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_trimmomatic_trim_filter_L00*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load trimmomatic/0.39

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/2.rename_sequences/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/3.trimmomatic_trim_filter/L00*/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the R1 files in the input directory
for file_R1 in "${input_dir}"*_R1.fastq.gz; do
    # Generate the corresponding R2 file name
    file_R2="${file_R1/_R1/_R2}"

      # Run Trimmomatic to determine paired reads and trim adapters
      trimmomatic PE \
        -threads 28 \
        -phred33 \
        "$file_R1" "$file_R2" \
        "${output_dir}$(basename "${file_R1}" .fastq.gz)_trim_filter_paired.fq.gz" \
        "${output_dir}$(basename "${file_R1}" .fastq.gz)_trim_filter_unpaired.fq.gz" \
        "${output_dir}$(basename "${file_R2}" .fastq.gz)_trim_filter_paired.fq.gz" \
        "${output_dir}$(basename "${file_R2}" .fastq.gz)_trim_filter_unpaired.fq.gz" \
        ILLUMINACLIP:"${input_dir}"TruSeq3-PE-2.fa:2:30:10:1:True MINLEN:40
done
```

## 2.trim_align_sequences

### 4.bwamem2_alignment

#### Align your sequences to your reference genome

- Takes 3 days, less time for certain lanes
- bwa-mem2 v2.2.1
- bwa-mem manual <https://bio-bwa.sourceforge.net/bwa.shtml>
- M is required for using Picard which I use later to mark duplicates
- R names files
- t is threads
- I include the path to the ref files I created above and the R1 and R2
  files that would be mapped onto the reference
- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J bwamem2_alignment_L00*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 48G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_bwamem2_alignment_L00*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_bwamem2_alignment_L00*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load bwa-mem2/2.2.1

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/3.trimmomatic_trim_filter/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/4.bwamem2_alignment/L00*/"
ref_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the R1 files in the input directory
for file_R1 in "${input_dir}"*_R1_trim_filter_paired.fq.gz; do
    # Generate the corresponding R2 file name
    file_R2="${file_R1/_R1/_R2}"
    # Extract the sample name from the file name
    sample=$(basename "$file_R1" | sed 's/_[^_]*_[^_]*_[^_]*$//')

      # Align the reads to the reference genome
      bwa-mem2 mem \
        -M \
        -R "@RG\tID:${sample}_L00*\tSM:${sample}\tPL:illumina" \
        -t 64 \
        ${ref_dir}ref \
        "$file_R1" "$file_R2" > "${output_dir}$(basename "${file_R1}" R1_paired.fq.gz)alignment.sam"
done
```

## 3.edit_alignments

### 1.samtools_view_sort

#### Transform SAM to BAM and sort sequences

- Takes 6 hours
- samtools v1.13
- samtools manual <https://www.htslib.org/doc/samtools.html>
- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J 1.samtools_view_sort_L00*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_1.samtools_view_sort_L00*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_1.samtools_view_sort_L00*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.trim_align_sequences/4.bwamem2_alignment/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/1.samtools_view_sort/"

# Iterate over the files in the input directory
for file in "${input_dir}"*_alignment.sam; do

      # Filter out reads with low mapping quality
      samtools view \
        -h \
        -@64 \
        -b "$file" |
      samtools view -buS - |
      samtools sort \
        -O bam,level=1 \
        -@64 \
        -o "${output_dir}$(basename "${file}" _bwamem2_alignment.sam)view_sort.bam" \
        -T $output_dir
done
```

## 3.edit_alignment

### 2.picard_markdups

#### Mark duplicate sequences

- Takes 1 day
- picard v2.26.2
- Input lane files
- Output library file (required)
- Metrics file (required)

``` bash
#!/bin/bash

#SBATCH -J 2.picard_markdups.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 1-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_2.picard_markdups_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_2.picard_markdups_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load  module
module load picard/2.26.2

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/1.samtools_view_sort/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/2.picard_markdups/"

# Iterate over the files in the input directory
for file_L001 in "${input_dir}"*_L001_view_sort.bam; do
    # Generate file_L002
    file_L002="${file_L001/_L001_/_L002_}"
    # Generate file_L003
    file_L003="${file_L001/_L001_/_L003_}"
    # Generate file_L004
    file_L004="${file_L001/_L001_/_L004_}"

      # Remove duplicates
      picard MarkDuplicates \
          -I "$file_L001" \
        -I "$file_L002" \
        -I "$file_L003" \
        -I "$file_L004" \
        -O "${output_dir}$(basename "$file_L001" L001_view_sort.bam)markdups.bam" \
          -M "${output_dir}$(basename "$file_L001" L001_view_sort.bam)markdups_metrics.txt" \
        -VALIDATION_STRINGENCY SILENT \
        -REMOVE_DUPLICATES true
done
```

## x.check_files

### 1.samtools_flagstat.sh

#### Calculate flag statistics

- Takes 2 hours
- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

#SBATCH -J samtools_flagstat.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_samtools_flagstat_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_samtools_flagstat_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/2.picard_markdups/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/samtools_flagstat/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in ${input_dir}*_marked_dups.bam; do
    samtools flagstat \
        -O tsv $file > ${output_dir}$(basename "$file" _marked_dups.bam)_flagstat.out
done
```

- Download files

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/flagstat /Users/bailey/Documents/research/S_orbicularis/DNA/results/
```

## 3.edit_alignment

### 3.bamutils_clip

#### Clip overlapping sequences

- Takes 13 hours
- bamutil v1.0.15

``` bash
#!/bin/bash

#SBATCH -J 3.bamutils_clip.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-13:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.bamutils_clip_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.bamutils_clip_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load bamutil/1.0.15

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/2.picard_markdups/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/3.bamutils_clip/"

# Iterate over the files in the input directory
for file in "${input_dir}"*_markdups.bam; do

      # Clip overlapping ends of read pairs
      bam clipOverlap \
            --in "$file" \
            --out "${output_dir}$(basename "${file}" markdups.bam)clip.bam" \
            --stats
done
```

## 3.edit_alignment

### 4a.gatk_target_creator

#### Identify where realignment around indels will occur

- Takes 1 day and 12 hours
- gatk v3.7, java v1.8, samtools v1.13
- <https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md>

``` bash
#!/bin/bash

#SBATCH -J 4a.gatk_target_creator.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 1-12:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_4a.gatk_target_creator_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_4a.gatk_target_creator_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/3.bamutils_clip/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
gatk="/home/bcarlson4/local/gatk3.7.0/GenomeAnalysisTK.jar"
java="/home/bcarlson4/local/jdk1.8.0_121/bin/java"

# Make the list file
rm -r /home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list
touch ${input_dir}"bamutils_clip_bam.list"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list"

# Iterate over the files in the input directory
for sample in "${input_dir}"*_clip.bam; do
    # Append the full path to the file list
    echo $sample >> $file_list
done

## Loop over each sample
cd $input_dir
for sample in `cat $file_list`; do
    samtools index $sample
done

      ## Realign around indels
      ## This is done across all samples at once
      ## Create list of potential indels
      $java -Xmx40g -jar $gatk \
        -T RealignerTargetCreator \
        -R $ref \
        -I $file_list \
        -o "${output_dir}indel_realigner.intervals" \
        -drf BadMate
```

## 3.edit_alignment

### 4b.gatk_indel_realigner

#### Realign sequences around indels

- Takes 6 days
- gatk v3.7, java v1.8
- <https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md>

``` bash
#!/bin/bash

#SBATCH -J 4b.gatk_indel_realigner.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 6-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p dept.les # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_4b.gatk_indel_realigner_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_4b.gatk_indel_realigner_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Assuming your files are in the current directory
dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
gatk="/home/bcarlson4/local/gatk3.7.0/GenomeAnalysisTK.jar"
java="/home/bcarlson4/local/jdk1.8.0_121/bin/java"

# Use the list file
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list"

cd $output_dir

      ## Run the indel realigner tool
      $java -Xmx40g -jar $gatk \
        -T IndelRealigner \
        -R $ref \
        -I $file_list \
        -targetIntervals ${dir}indel_realigner.intervals \
        --consensusDeterminationModel USE_READS  \
        --nWayOut _indel_realigner.bam
```

#### Copy output bam files into population folders

``` bash
#!/bin/bash

#SBATCH -J scp_files_to_lists.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_scp_files_to_lists_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_scp_files_to_lists_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"

mkdir ${dir}BAL ${dir}CCM ${dir}CCN ${dir}CLM ${dir}CRF ${dir}FLK ${dir}LMU ${dir}NCO ${dir}NLK ${dir}OCK ${dir}OCM ${dir}OTM ${dir}TKC ${dir}TLN ${dir}ULN

scp M0D039797C* M0D039798D* M0D039920V* M0D039921W* M0D039922X* M0D039923Y* M0D039924Z* M0D039925A* ${dir}BAL
scp M0D040019Q* M0D040020R* M0D040021S* M0D040022T* M0D040023U* M0D040024V* M0D040025W* M0D040026X* M0D040028Z* M0D040031C* ${dir}CCM
scp M0D039950Z* M0D039951A* M0D039952B* M0D039953C* M0D039954D* M0D039955E* M0D039956F* M0D039957G* ${dir}CCN
scp M0D039846Z* M0D039847A* M0D039848B* M0D039852F* M0D039855I* M0D039856J* M0D039970T* M0D040011I* M0D040012J* M0D040014L* M0D040017O* M0D040018P* ${dir}CLM
scp M0D039800F* M0D039859M* M0D039860N* M0D039861O* M0D039862P* M0D039913O* M0D039914P* M0D039915Q* M0D039916R* M0D039917S* M0D039937M* M0D039872Z* ${dir}CRF
scp M0D039864R* M0D039865S* M0D039866T* M0D039867U* M0D039868V* M0D039926B* M0D039927C* M0D039928D* M0D039929E* M0D039931G* M0D039932H* M0D039933I* ${dir}FLK
scp M0D039958H* M0D039959I* M0D039960J* M0D039961K* M0D039962L* M0D039963M* M0D039964N* M0D039965O* M0D039966P* M0D039967Q* ${dir}LMU
scp M0D040047S* M0D040048T* M0D040049U* M0D040050V* M0D040051W* M0D040052X* M0D040053Y* ${dir}NCO
scp M0D039884L* M0D039886N* M0D039888P* M0D039889Q* M0D039891S* M0D039892T* M0D039894V* M0D039895W* M0D039898Z* M0D039900B* M0D039904F* M0D039905G* ${dir}NLK
scp M0D039858L* M0D039907I* M0D039908J* M0D039909K* M0D039910L* M0D039911M* M0D039912N* M0D039934J* M0D039935K* M0D039936L* ${dir}OCK
scp M0D039878F* M0D039879G* M0D039880H* M0D039881I* M0D040034F* M0D040035G* M0D040037I* M0D040038J* M0D040039K* M0D040041M* M0D040042N* M0D040043O* ${dir}OCM
scp M0D039802H* M0D039803I* M0D039804J* M0D039805K* M0D039810P* M0D039811Q* M0D039812R* M0D039813S* M0D039820Z* M0D039873A* M0D039875C* M0D039877E* ${dir}OTM
scp M0D039944T* M0D039946V* M0D039947W* M0D039948X* M0D039949Y* ${dir}TKC
scp M0D039828H* M0D039829I* M0D039830J* M0D039831K* M0D039832L* M0D039833M* M0D039835O* M0D039837Q* M0D039838R* M0D039841U* M0D039842V* M0D039844X* ${dir}TLN
scp M0D039918T* M0D039919U* M0D039938N* M0D039939O* M0D039940P* M0D039941Q* M0D039943S* M0D040027Y* M0D040029A* M0D040030B* M0D040033E* M0D040045Q* ${dir}ULN

mkdir ${dir}Ocean ${dir}Stratified ${dir}Mixed

scp M0D039944T* M0D039946V* M0D039947W* M0D039948X* M0D039949Y* M0D039878F* M0D039879G* M0D039880H* M0D039881I* M0D040034F* M0D040035G* M0D040037I* M0D040038J* M0D040039K* M0D040041M* M0D040042N* M0D040043O* M0D039858L* M0D039907I* M0D039908J* M0D039909K* M0D039910L* M0D039911M* M0D039912N* M0D039934J* M0D039935K* M0D039936L* M0D040052X* M0D040053Y* M0D040047S* M0D040048T* M0D040049U* M0D040050V* M0D040051W* M0D039958H* M0D039959I* M0D039960J* M0D039961K* M0D039962L* M0D039963M* M0D039964N* M0D039965O* M0D039966P* M0D039967Q* M0D039800F* M0D039859M* M0D039860N* M0D039861O* M0D039862P* M0D039913O* M0D039914P* M0D039915Q* M0D039916R* M0D039917S* M0D039937M* M0D039872Z* M0D039950Z* M0D039951A* M0D039952B* M0D039953C* M0D039954D* M0D039955E* M0D039956F* M0D039957G* M0D040019Q* M0D040020R* M0D040021S* M0D040022T* M0D040023U* M0D040024V* M0D040025W* M0D040026X* M0D040028Z* M0D040031C* M0D039797C* M0D039798D* M0D039920V* M0D039921W* M0D039922X* M0D039923Y* M0D039924Z* M0D039925A* ${dir}Ocean
scp M0D039828H* M0D039829I* M0D039830J* M0D039831K* M0D039832L* M0D039833M* M0D039835O* M0D039837Q* M0D039838R* M0D039841U* M0D039842V* M0D039844X* M0D039802H* M0D039803I* M0D039804J* M0D039805K* M0D039810P* M0D039811Q* M0D039812R* M0D039813S* M0D039820Z* M0D039873A* M0D039875C* M0D039877E* M0D039884L* M0D039886N* M0D039888P* M0D039889Q* M0D039891S* M0D039892T* M0D039894V* M0D039895W* M0D039898Z* M0D039900B* M0D039904F* M0D039905G* M0D039846Z* M0D039847A* M0D039848B* M0D039852F* M0D039855I* M0D039856J* M0D039970T* M0D040011I* M0D040012J* M0D040014L* M0D040017O* M0D040018P* ${dir}Stratified
scp M0D039918T* M0D039919U* M0D039938N* M0D039939O* M0D039940P* M0D039941Q* M0D039943S* M0D040027Y* M0D040029A* M0D040030B* M0D040033E* M0D040045Q* M0D039864R* M0D039865S* M0D039866T* M0D039867U* M0D039868V* M0D039926B* M0D039927C* M0D039928D* M0D039929E* M0D039931G* M0D039932H* M0D039933I* ${dir}Mixed

## Make the list files
# Site
for SITE in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN
do
echo $SITE
touch ${dir}$SITE/$SITE'_angsd_site_bam.list'
    # Iterate over the files in the input directory
  for sample in ${dir}$SITE/*_clipped_indel_realigner.bam; do
      # Append the full path to the file list
      echo $sample >> ${dir}$SITE/$SITE'_angsd_site_bam.list'
  done
done

# Type
for TYPE in Ocean Mixed Stratified
do
echo $TYPE
touch ${dir}$TYPE/$TYPE'_angsd_type_bam.list'
    # Iterate over the files in the input directory
  for sample in ${dir}$TYPE/*_indel_realigner.bam; do
      # Append the full path to the file list
      echo $sample >> ${dir}$TYPE/$TYPE'_angsd_type_bam.list'
  done
done
```

## x.check_files

### 2.angsd_dodepth.sh

#### Get depth of loci

- 1 day & 12 hours

``` bash
#!/bin/bash

#SBATCH -J angsd_dodepth.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 1-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_angsd_dodepth_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_angsd_dodepth_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/2.angsd_dodepth/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Make the list file
rm -r ${input_dir}gatk_indel_realigner_bam.list
touch ${input_dir}"gatk_indel_realigner_bam.list"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/gatk_indel_realigner_bam.list"

# Iterate over the files in the input directory
for sample in "${input_dir}"*clipped_indel_realigner.bam; do
    # Append the full path to the file list
    echo $sample >> $file_list
done

        # Determine the per-position depth distribution
        $angsd -b $file_list -ref $ref -out ${output_dir}S_orbicularis_depth \
        -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -GL 1 -doMajorMinor 4 -doMaf 1 -P 10
```

#### Site/Type

``` bash
#!/bin/bash

#SBATCH -J angsd_dodepth_site_type.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 1-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_angsd_dodepth_site_type_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_angsd_dodepth_site_type_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/2.angsd_dodepth/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN; do
        # Determine the per-position depth distribution
        $angsd -b ${input_dir}$POP/$POP'_angsd_site_bam.list' -ref $ref -out ${output_dir}$POP'S_orbicularis_depth' \
        -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -GL 1 -doMajorMinor 4 -doMaf 1 -P 10
done

for POP in Ocean Mixed Stratified; do
        # Determine the per-position depth distribution
        $angsd -b ${input_dir}${POP}/${POP}'_angsd_type_bam.list' -ref $ref -out ${output_dir}$POP'S_orbicularis_depth' \
        -doCounts 1 -doDepth 1 -maxDepth 10000 -dumpCounts 1 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -GL 1 -doMajorMinor 4 -doMaf 1 -P 10
done
```

- Download files

``` bash
# Download Global depth
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/angsd_dodepth/S_orbicularis_depth.depthGlobal .

# Download Sample depth
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/angsd_dodepth/S_orbicularis_depth.depthSample .
```

#### Depth of loci plot

``` r
library(tidyverse)
library(ggplot2)

# Define basedir
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/data/"

# These files will give the depth at 
dg <- as.matrix(read.table(paste0(basedir, "S_orbicularis_depth.depthGlobal"), header = F))

ds <- as.matrix(read.table(paste0(basedir, "S_orbicularis_depth.depthSample"), header = F))

dgt <- t(dg)
dgt <- as.data.frame(dgt)
row.names(dgt) <- NULL
dgt$Depth <- c(0:10000)

dgt$V1 <- as.numeric(dgt$V1)
dgt$Depth <- as.numeric(dgt$Depth)

str(dgt)

# Create the histogram with a log y-axis and overlay the normal distribution
hist_plot <- ggplot(data = dgt, mapping = aes(x = Depth, y = V1)) + 
  geom_bar(
    stat = "identity",          # Use the values as is (no additional counting)
    color = "black",            # Outline color of the bars
    fill = "lightblue"          # Fill color of the bars
  ) +
  theme_bw() +                  # Use the black-and-white theme
  theme(
    text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  labs(
    x = "Read Depth", 
    y = "Number of Sites"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(dgt$Depth), by = 1000)  # Adjust breaks as needed
  )
hist_plot
hist_plot <- hist_plot + scale_y_log10(                                  # Apply log scale to y-axis
    breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),  # Custom breaks
    labels = scales::comma                        # Format labels with commas
  ) +
  labs(
    y = "Number of Sites (log scale)"
  )
hist_plot

# Assuming you want to fit a normal distribution around the peak (depth ~361)
peak_data <- dgt[dgt$Depth >= 138 & dgt$Depth <= 582,]
# peak_data <- dgt[dgt$Depth >= 0 & dgt$Depth <= 720,]

# Assuming 'peak_mean' and 'peak_sd' are the mean and standard deviation of the peak
peak_mean <- mean(peak_data$Depth)  # Replace with your calculated mean
peak_sd <- sd(peak_data$Depth)      # Replace with your calculated standard deviation

# Calculate depths 2 standard deviations away from the mean
low_2sd <- peak_mean - 1 * peak_sd
high_2sd <- peak_mean + 1 * peak_sd

# Display the results
cat("Depth 1 standard deviations to the left of the mean:", low_2sd, "\n")
cat("Depth 1 standard deviations to the right of the mean:", high_2sd, "\n")

# Generate the x-values (depth range) for plotting the normal distribution
x_values <- seq(min(dgt$Depth), max(dgt$Depth), length.out = 440)

# Calculate the corresponding y-values using dnorm
y_values <- dnorm(x_values, mean = peak_mean, sd = peak_sd)

# Plot the normal distribution
normal_plot <- ggplot(data = data.frame(x_values, y_values), aes(x = x_values, y = y_values)) +
  geom_line(color = "red", linewidth = 1) +
  theme_minimal() +
  labs(x = "Read Depth", y = "Density") +
  scale_x_continuous(breaks = c(0, 500, 1000), limits = c(0,1000))

# Display the plot
print(normal_plot)
```

## x.check_files

### 3.samtools_coverage.sh

#### Produces a histogram or table of coverage per chromosome

- Takes 4 hours
- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

#SBATCH -J samtools_coverage.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-04:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_samtools_coverage_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_samtools_coverage_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/3.samtools_coverage/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in ${input_dir}*_clipped_indel_realigner.bam; do
    samtools coverage -o ${output_dir}$(basename "$file" _clipped_indel_realigner.bam)_coverage.out $file
done

# Iterate over the files in the input directory
    samtools coverage -o ${output_dir}all_coverage.out ${input_dir}*_clipped_indel_realigner.bam
```

- Download file

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/samtools_coverage .
```

## x.check_files

### 4.samtools_depth.sh

#### Computes read depth at each position or region

- Takes 8 hours
- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

#SBATCH -J samtools_depth.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-09:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_samtools_depth_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_samtools_depth_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/4.samtools_depth/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in ${input_dir}*_clipped_indel_realigner.bam; do
    samtools depth -aa $file | cut -f 3 | gzip > ${output_dir}$(basename "$file" _clipped_indel_realigner.bam)_depth.out.gz
done
```

- Download file

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/samtools_depth /Users/bailey/Documents/research/S_orbicularis/DNA/results
```

#### Create depth plots

``` r
## Install tidyverse and ggplot2 if you don't have it installed yet
library(tidyverse)
library(ggplot2)

basedir <- "/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/"
bam_list <- read_lines(paste0(basedir, "sample_list.txt"))

for (i in 1:length(bam_list)){
  bamfile = bam_list[i]
  # Compute depth stats
  depth <- read_tsv(paste0(basedir, "/4.samtools_depth/", bamfile, "depth"), col_names = F)$X1
  mean_depth <- mean(depth)
  sd_depth <- sd(depth)
  mean_depth_nonzero <- mean(depth[depth > 0])
  mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
  median <- median(depth)
  presence <- as.logical(depth)
  proportion_of_reference_covered <- mean(presence)
  output_temp <- tibble(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)

  # Bind stats into dataframe and store sample-specific per base depth and presence data
  if (i==1){
    output <- output_temp
    total_depth <- depth
    total_presence <- presence
  } else {
    output <- bind_rows(output, output_temp)
    total_depth <- total_depth + depth
    total_presence <- total_presence + presence
  }
}

output <- output %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

summary(total_depth)

# Plot the depth distribution (this may take a few minutes to run)
depth_distribution <- tibble(total_depth = total_depth, position = 1:length(total_depth))
depth_distribution

depth_distribution_plot <- ggplot(depth_distribution, aes(x = position, y = total_depth)) +
  geom_point(size = 0.1) +
  theme_bw()
ggsave("/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/4.samtools_depth/depth_distribution.png", depth_distribution_plot, width = 8, height = 6)

# Total depth per site across all individuals 
total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
str(total_depth_summary)
total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)
str(total_presence_summary)

total_depth_summary_plot1 <- ggplot(total_depth_summary, aes(x = total_depth, y = n)) +
  geom_point()
ggsave("/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/4.samtools_depth/total_depth_summary1.png", total_depth_summary_plot1, width = 8, height = 6)

total_depth_summary_plot2 <- ggplot(total_depth_summary, aes(x = total_depth, y = n)) +
  geom_point() +
  coord_cartesian(xlim=c(NA, 20))
ggsave("/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/4.samtools_depth/total_depth_summary2.png", total_depth_summary_plot2, width = 8, height = 6)

total_presence_summary_plot <- ggplot(total_presence_summary, aes(x = total_presence, y = n)) +
  geom_col()
ggsave("/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.check_files/4.samtools_depth/total_presence_summary.png", total_presence_summary_plot, width = 8, height = 6)
```

## 4.determine_allele_frequencies

### 1.angsd_gl_snp

- Takes 1.5 days
- angsd v0.94, samtools v1.13, gcc, gsl
- Replace bam_list.txt with a file containing the paths to your aligned
  and deduplicated BAM files. -b is a list of your bam files -ref is the
  reference sequence -out is the output file destination with name you
  want to give output files -uniqueOnly retains only uniquely mapped
  reads -remove_bads whether to keep only reads that were not tagged as
  bad, we removed bad reads -only_proper_pairs whether to include only
  properly paired sequences, we kept only proper paired reads -trim
  whether to trim, we chose not to trim -C reduces the effect of
  excessive mismatches -baq computes base alignment quality -minMapQ the
  minimum mapping quality score reads must have to be kept, we will
  filter after ngsParalog -minQ the minimum quality score reads must
  have to be kept, we will filter after ngsParalog -minInd only keep
  sites with a minimum number of individuals, used ~90% of individuals
  -setMinDepth filter out sites with less than this number of reads,
  used info from 2.angsd_dodepth for this -setMaxDepth filter out sites
  with more than this number of reads, did not set a max filter because
  we use ngsParalog -doCounts calculate various counts statistics -doMaf
  estimate allele frequencies -minMaf remove sites with MAF below the
  specified number, will not filter until 5.angsd_gl_snp_maf_nPc
  -SNP_pval remove sites with a pvalue larger than the specified number
  -GL estimate genotype likelihoods, we chose to use the samtools method
  -doGlf create likelihood file, we chose a beagle likelihood file
  -doMajorMinor specify how to assign the major and minor alleles, we
  used the reference allele as major (4) but later we will use the -site
  file to use SNPs (3) -doPost calculate the posterior probability of
  genotypes, we used a uniform prior -nthreads the number of threads to
  use for the processing, max is 10

``` bash
#!/bin/bash

#SBATCH -J 1.angsd_gl_snp.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 1-18:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_1.angsd_gl_snp_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_1.angsd_gl_snp_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/gatk_indel_realigner_bam.list"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/1.angsd_gl_snp/S_orbicularis_gl_snp_140_230/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

      # Call allele frequencies, genotype likelihoods, and SNPs
      $angsd -b $file_list -ref $ref -anc $ref -out ${output_dir}S_orbicularis_gl_snp_140_230 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -doCounts 1  -minInd 140 -setMinDepth 230 -GL 1 -doGlf 2 -doMajorMinor 4 \
        -doMaf 1 -doHWE 1-doSaf 1 -doIBS 1 -doCov 1 -makeMatrix 1 -nthreads 10

zcat ${output_dir}S_orbicularis_gl_snp_140_230.mafs.gz | cut -f 1-4 | tail -n+2 > ${output_dir}S_orbicularis_gl_snp_140_230.snps

$angsd sites index ${output_dir}S_orbicularis_gl_snp_140_230.snps
```

- Download file

``` bash
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/1.angsd_gl_snp/S_orbicularis_gl_snp_140_230.covMat /Users/bailey/Documents/research/S_orbicularis/DNA/results/gl_SNP/
```

#### Run a PCA based on your output SNPs

``` r
library(tidyverse)
# Define basedir
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/results/gl_SNP/"

# Load the covariance matrix
cov <- as.matrix(read.table(paste0(basedir, "S_orbicularis_gl_snp_140_230.covMat"), header = F))

# We will also add a column with population assingments
pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/S_orbicularis_154_samples.csv")
     
pca <- eigen(cov) # perform the pca using the eigen function. 
summary(pca)

eigenvectors = pca$vectors # extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) 
# combine with our population assignments

pca.eigenval.sum = sum(pca$values) #sum of eigenvalues
varPC1 <- (pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

# Biplot
plot(eigenvectors[,1:2])

# Run PCA
pca <- prcomp(cov, center = TRUE, scale. = TRUE)
summary(pca)
    
# Biplot
biplot(pca)

# Scree plot
plot(pca)

# Pull out principal component data frame
pca_coords <- pca$x

pca_coords <- cbind(pca_coords, pop)

# Loadings of variables on PC1
pca_loadings <- as.data.frame(pca$rotation)

# Type
pca_coords$Type <- factor(pca_coords$Type, levels = c("Ocean", "Mixed", "Stratified"))

custom_colors <- c("Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

# Site
pca_coords$Site <- factor(pca_coords$Site, levels = c("BAL", "CCM", "CCN", "CRF", "LMU", "NCO", "OCK", "OCM", "TKC", "FLK", "ULN", "CLM", "NLK", "OTM", "TLN"))

custom_colors <- c("BAL" = "#BE3428FF", "CCM" = "#8E2322FF", "CCN" = "#D9565CFF", "CRF" = "#FF3200FF", "LMU" = "#E9A17CFF", "NCO" = "#881C00FF", "OCK" = "#D6604DFF", "OCM" = "#ED3F39FF", "TKC" = "#803233FF", "FLK" = "#69D2E7FF", "ULN" = "#639CA4FF", "CLM" = "#204035FF", "NLK" = "#6D8325FF", "OTM" = "#417839FF", "TLN" = "palegreen3")

pca_plot <- ggplot(data = pca_coords, aes(x=PC1, y=PC2, colour = Site, fill = Site)) + 
  geom_jitter(size = 5, alpha = 0.75) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  scale_y_continuous(breaks = c(-10,-5,0,5,10), limits = c(-12,12)) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10), limits = c(-12,12)) +
  # scale_y_continuous(breaks = c(-5,0,5,10), limits = c(-7,13)) +
  # scale_x_continuous(breaks = c(-5,0,5,10,15,20,25), limits = c(-7,25)) +
labs(x= "PC1 2%", y= "PC2 2%", colour = "Site type", fill = "Site type")
pca_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/SNPs_PCA/pca_S_orbicularis_gl_snp_140_230_plot.png", pca_plot, width = 8, height = 6)
```

#### Split the SNPs by chromosome for ngsParalog

``` bash
#!/bin/bash

#SBATCH -J split_by_chromosome.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_split_by_chromosome_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_split_by_chromosome_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Directory to store files
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/1.angsd_gl_snp/S_orbicularis_gl_snp_140_230/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Path to the input file
SNP_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/1.angsd_gl_snp/S_orbicularis_gl_snp_140_230/S_orbicularis_gl_snp_140_230.snps"

# Define the valid chromosome numbers
for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24 mtDNA; do

    # Pad the chromosome number with leading zero if necessary
    output_SNP_file="S_orbicularis_gl_snp_140_230_${chr}.snps"

    # Use awk to filter lines for the current chromosome
    awk -v chr="${chr}" -v out="${output_SNP_file}" '$1 == chr { print > out }' $SNP_file

    # Count the number of lines in the created file
    line_count=$(wc -l < "$output_SNP_file")
    echo "Number of lines in ${output_SNP_file}: ${line_count}"
done
```

## 4.determine_allele_frequencies

### 2.ngsP_prune

#### Identifying regions like paralogy and repetitive sequences that confound short read mapping

- Takes 6 hours
- <https://github.com/tplinderoth/ngsParalog>

``` bash
#!/bin/bash

#SBATCH -J 2.ngsP_prune_140_230.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=1 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_2.ngsP_prune_140_230_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_2.ngsP_prune_140_230_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/gatk_indel_realigner_bam.list"
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/1.angsd_gl_snp/S_orbicularis_gl_snp_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nP_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ngsParalog and ANGSD
ngsParalog="/home/bcarlson4/local/ngsTools/ngsParalog/ngsParalog"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# convert site file to bed for ngsParalog
awk '{print $1"\t"$2-1"\t"$2}' ${input_dir}S_orbicularis_gl_snp_140_230.snps > ${output_dir}S_orbicularis_gl_snp_140_230.bed

### Run mpileup and ngsParalog without intermediate files
samtools mpileup -b $file_list -l ${output_dir}S_orbicularis_gl_snp_140_230.bed -q 0 -Q 0 --ff UNMAP,DUP |
$ngsParalog calcLR \
    -infile - \
    -outfile ${output_dir}S_orbicularis_gl_snp_nP_140_230.snps \
    -minQ 20 -minind 140 -allow_overwrite 1

### Convert ngsparalog output in list of canonical and deviant SNPs based on p-value threshold
module load anaconda3
source activate my-R
pval=0.001
Rscript convert_ngsparalog_to_sitelist.R \
    ${output_dir}S_orbicularis_gl_snp_nP_140_230.snps \
    ${input_dir}S_orbicularis_gl_snp_140_230.snps $pval
conda deactivate

# Copy files for later use
cp ${output_dir}S_orbicularis_gl_snp_nP_140_230.snps_deviant ${output_dir}S_orbicularis_gl_snp_nPd_140_230.snps
cp ${output_dir}S_orbicularis_gl_snp_nP_140_230.snps_canonical ${output_dir}S_orbicularis_gl_snp_nPc_140_230.snps

# Index files for ANGSD
$angsd sites index ${output_dir}S_orbicularis_gl_snp_nPd_140_230.snps
$angsd sites index ${output_dir}S_orbicularis_gl_snp_nPc_140_230.snps
```

#### Combine the out files from ngsParalog

- Only if you separated SNPs out by chromosome.

``` bash
#!/bin/bash

#SBATCH -J cat.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=1 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_cat_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_cat_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Directory
dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nP_140_230/"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Concatenate the files to get all SNPs together
cat ${dir}S_orbicularis_gl_snp_nP_140_230_*.snps_deviant > ${dir}S_orbicularis_gl_snp_nPd_140_230.snps

# Concatenate the files to get all SNPs together
cat ${dir}S_orbicularis_gl_snp_nP_140_230_*.snps_canonical > ${dir}S_orbicularis_gl_snp_nPc_140_230.snps

# Index file for ANGSD
$angsd sites index ${dir}S_orbicularis_gl_snp_nPc_140_230.snps
```

## 4.determine_allele_frequencies

### 3.angsd_saf

#### Sample allele frequency

- Takes 2 days and 5 hours
- angsd v0.94, samtools v1.13, gcc, gsl -doSaf 1: perform multisample GL
  estimation \#### Site

``` bash
#!/bin/bash

#SBATCH -J 3.angsd_saf_site.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.angsd_saf_site_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.angsd_saf_site_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nP_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/3.angsd_saf/S_orbicularis_gl_snp_nPc_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN;
do
        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_site_bam.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 \
            -doSaf 1 -doMaf 1 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 10

zcat ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230'.mafs.gz | cut -f 1-4 | tail -n+2 > ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230'.snps

$angsd sites index ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230'.snps
done
```

#### Type

``` bash
#!/bin/bash

#SBATCH -J 3.angsd_saf_type.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.angsd_saf_type_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.angsd_saf_type_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nP_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/3.angsd_saf/S_orbicularis_gl_snp_nPc_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

for POP in Ocean Mixed Stratified;
do
        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_type_bam.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 \
            -doSaf 1 -doMaf 1 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 5

zcat ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230'.mafs.gz | cut -f 1-4 | tail -n+2 > ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230'.snps

$angsd sites index ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230'.snps
done
```

#### Subsample

``` bash
#!/bin/bash

#SBATCH -J 3.angsd_saf_type_subsampled.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.angsd_saf_type_subsampled_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.angsd_saf_type_subsampled_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nP_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/3.angsd_saf/S_orbicularis_gl_snp_nPc_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

for POP in Ocean Mixed;
do

# Subsample; do not recreate bam list if you already have done so
# shuf -n 12 ${input_dir}$POP/$POP'_angsd_type_bam.list' > ${input_dir}$POP/$POP'_angsd_type_bam_subsampled.list'

        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_type_bam_subsampled.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 \
            -doSaf 1 -doMaf 1 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 5

zcat ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled'.mafs.gz | cut -f 1-4 | tail -n+2 > ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled'.snps

$angsd sites index ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled'.snps
done
```

## 4.determine_allele_frequencies

### 4.winsfs

### Site frequency spectrum

- Takes days
- winsfs, samtools v1.13, gcc, gsl \#### Site/Type/Subsample

``` bash
#!/bin/bash

#SBATCH -J 4.winsfs.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_4.winsfs_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_4.winsfs_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/3.angsd_saf/S_orbicularis_gl_snp_nPc_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/4.winsfs/S_orbicularis_gl_snp_nPc_140_230/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Site
# Run a for loop to calculate site frequency spectrums (SFS)
for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN; do
      
      #SFS
      winsfs -v ${input_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' -t 10 > ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${output_dir}$POP'_fold_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${output_dir}$POP'_fold_S_orbicularis_gl_snp_nPc_140_230.sfs' ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs'
done

populations=("BAL" "CCM" "CCN" "CLM" "CRF" "FLK" "LMU" "NCO" "NLK" "OCK" "OCM" "OTM" "TKC" "TLN" "ULN")

# Loop over populations. Arrays start at 0, which is why j will still be less than 15 when it gets to ULN
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # SFS
      winsfs ${input_dir}${pop1}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${input_dir}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' -v -t 10 \
        > ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'
    done
done

## Type
# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in Ocean Mixed Stratified; do
      
      #SFS
      winsfs -v ${input_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' -t 10 > ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${output_dir}$POP'_fold_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${output_dir}$POP'_fold_S_orbicularis_gl_snp_nPc_140_230.sfs' ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs'
done

populations=("Ocean" "Mixed" "Stratified")

# Loop over populations, arrays start at 0 which is why j will still be less than 3 when it gets to Ocean
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # SFS
      winsfs ${input_dir}${pop1}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${input_dir}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' -v -t 10 \
        > ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'
    done
done

# Subsample
# Run a for loop to calculate site frequency spectrums (SFS)
for POP in Ocean Mixed; do

      # SFS
      winsfs -v ${input_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.saf.idx' -t 10 > ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs'

      winsfs view -v -f ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs'

      tail -n +2 ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs' > ${output_dir}$POP'_fold_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs'
      mv ${output_dir}$POP'_fold_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs' ${output_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs'
done

# Ocean vs. Marine Lakes
# Loop over populations
for POP in CLM NLK OTM TLN FLK ULN; do

      # SFS
      winsfs ${input_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${input_dir}'Ocean_S_orbicularis_gl_snp_nPc_140_230_subsampled.saf.idx' -v -t 10 \
        > ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${output_dir}$POP'Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${output_dir}$POP'Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs' ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs'
done

# Mixed vs. Stratified Lakes
# Loop over populations
for POP in CLM NLK OTM TLN; do

      # SFS
      winsfs ${input_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${input_dir}'Mixed_S_orbicularis_gl_snp_nPc_140_230_subsampled.saf.idx' -v -t 10 \
        > ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${output_dir}$POP'Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${output_dir}$POP'Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs' ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs'
done
```

## 4.determine_allele_frequencies

### 5.angsd_gl_snp_maf_nPc

- Takes 1 day & 12 hours

``` bash
#!/bin/bash

#SBATCH -J 5.angsd_gl_snp_maf_nPc.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=4 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 1-18:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_gl_snp_maf_nPc_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_gl_snp_maf_nPc_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Assuming your files are in the current directory
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/gatk_indel_realigner_bam.list"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nPc_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

      # Call allele frequencies, genotype likelihoods, and SNPs
      $angsd -b $file_list -ref $ref -anc $ref -out ${output_dir}S_orbicularis_gl_snp_maf_nPc_140_230 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 140 -setMinDepth 230 -doCounts 1 -sites ${SNPs} \
        -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -GL 1 -doGlf 2 -doMajorMinor 3 -doPost 2 \
        -doIBS 1 -doCov 1 -makeMatrix 1 -doHWE 1 -nthreads 4

## Remove header from beagle file output and create new beagle file for ngsParalog
zcat ${output_dir}S_orbicularis_gl_snp_maf_nPc_140_230.beagle.gz | cut -f 4- | gzip  > ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle.gz

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle.gz > ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle

# Count how many lines are in the file
wc -l ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle

# Print how many columns are in the file
awk '{print NF}' S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle | sort -nu | tail -n 1

## Remove header from maf file output, replace semicolon with underscore, and create new maf file for ngsParalog
zcat ${output_dir}S_orbicularis_gl_snp_maf_nPc_140_230.mafs.gz | cut -f 1,2 | sed 's/:/_/g'| sed '1d' | gzip > ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos.gz

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos.gz > ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos

# Count how many lines are in the file
wc -l ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos
```

#### Site

``` bash
#!/bin/bash

#SBATCH -J 5.angsd_gl_snp_maf_nPc_site.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:01:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_gl_snp_maf_nPc_site_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_gl_snp_maf_nPc_site_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nPc_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call allele frequencies, genotype likelihoods, and SNPs
for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN;
do
        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_site_bam.list' \
            -ref $ref -anc $ref -sites $SNPs -r chr01:chr24 \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 -doPost 2 \
            -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 10

## Remove header from beagle file output and create new beagle file for ngsParalog
zcat ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.beagle.gz' | cut -f 4- | gzip  > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle.gz'

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle.gz' > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle'

# Count how many lines are in the file
wc -l ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle'

# Print how many columns are in the file
awk '{print NF}' $POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle' | sort -nu | tail -n 1

## Remove header from maf file output, replace semicolon with underscore, and create new maf file for ngsParalog
zcat ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.mafs.gz' | cut -f 1,2 | sed 's/:/_/g'| sed '1d' | gzip > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos.gz'

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos.gz' > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos'

# Count how many lines are in the file
wc -l ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos'
done
```

#### Type

``` bash
#!/bin/bash

#SBATCH -J 5.angsd_gl_snp_maf_nPc_type.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:01:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_gl_snp_maf_nPc_type_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_gl_snp_maf_nPc_type_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nPc_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for POP in Ocean Mixed Stratified;
do
        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_type_bam.list' \
            -ref $ref -anc $ref -sites $SNPs -r chr01:chr24 \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 -doPost 2 \
            -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 10

## Remove header from beagle file output and create new beagle file for ngsParalog
zcat ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.beagle.gz' | cut -f 4- | gzip  > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle.gz'

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle.gz' > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle'

# Count how many lines are in the file
wc -l ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle'

# Print how many columns are in the file
awk '{print NF}' $POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle' | sort -nu | tail -n 1

## Remove header from maf file output, replace semicolon with underscore, and create new maf file for ngsParalog
zcat ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.mafs.gz' | cut -f 1,2 | sed 's/:/_/g'| sed '1d' | gzip > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos.gz'

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos.gz' > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos'

# Count how many lines are in the file
wc -l ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos'
done
```

#### Subsample

``` bash
#!/bin/bash

#SBATCH -J 5.angsd_gl_snp_maf_nPc_type_subsampled.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:01:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_gl_snp_maf_nPc_type_subsampled_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_gl_snp_maf_nPc_type_subsampled_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/2.ngsP_prune/S_orbicularis_gl_snp_nPc_140_230/S_orbicularis_gl_snp_nPc_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for POP in Ocean Mixed; do

# Subsample
# shuf -n 12 ${input_dir}$POP/$POP'_angsd_type_bam.list' > ${input_dir}$POP/$POP'_angsd_type_bam_subsampled.list'

        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_type_bam_subsampled.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 -doPost 2 \
            -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 10

## Remove header from beagle file output and create new beagle file for ngsParalog
zcat ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.beagle.gz' | cut -f 4- | gzip  > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.beagle.gz'

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.beagle.gz' > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.beagle'

# Count how many lines are in the file
wc -l ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.beagle'

# Print how many columns are in the file
awk '{print NF}' $POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.beagle' | sort -nu | tail -n 1

## Remove header from maf file output, replace semicolon with underscore, and create new maf file for ngsParalog
zcat ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.mafs.gz' | cut -f 1,2 | sed 's/:/_/g'| sed '1d' | gzip > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.pos.gz'

# Unzip the file but keep the original gz copy of the file
gunzip -c ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.pos.gz' > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.pos'

# Count how many lines are in the file
wc -l ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_subsampled.pos'
done
```

#### Download cov file

``` bash
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_nP_140_230/S_orbicularis_gl_snp_maf_nPc_140_230.covMat /Users/bailey/Documents/research/S_orbicularis/DNA/results/gl_SNP/
```

#### Run a PCA based on your output SNPs

``` r
library(tidyverse)

# Define basedir
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/results/gl_snp/"

# Load the covariance matrix
cov <- as.matrix(read.table(paste0(basedir, "S_orbicularis_gl_snp_maf_nPc_140_230.covMat"), header = F))

# We will also add a column with population assingments
pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/S_orbicularis_154_samples.csv")
     
pca <- eigen(cov) # perform the pca using the eigen function. 
summary(pca)

eigenvectors = pca$vectors # extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) 
# combine with our population assignments

pca.eigenval.sum = sum(pca$values) #sum of eigenvalues
varPC1 <- (pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

# Biplot
plot(eigenvectors[,1:2])

# Run PCA
pca <- prcomp(cov, center = TRUE, scale. = TRUE)
summary(pca)
    
# Biplot
biplot(pca)

# Scree plot
plot(pca)

# Pull out principal component data frame
pca_coords <- pca$x

pca_coords <- cbind(pca_coords, pop)

# Loadings of variables on PC1
pca_loadings <- as.data.frame(pca$rotation)

# Type
pca_coords$Type <- factor(pca_coords$Type, levels = c("Ocean", "Mixed", "Stratified"))

custom_colors <- c("Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

# Site
pca_coords$Site <- factor(pca_coords$Site, levels = c("BAL", "CCM", "CCN", "CRF", "LMU", "NCO", "OCK", "OCM", "TKC", "FLK", "ULN", "CLM", "NLK", "OTM", "TLN"))

custom_colors <- c("BAL" = "#BE3428FF", "CCM" = "#8E2322FF", "CCN" = "#D9565CFF", "CRF" = "#FF3200FF", "LMU" = "#E9A17CFF", "NCO" = "#881C00FF", "OCK" = "#D6604DFF", "OCM" = "#ED3F39FF", "TKC" = "#803233FF", "FLK" = "#69D2E7FF", "ULN" = "#639CA4FF", "CLM" = "#204035FF", "NLK" = "#6D8325FF", "OTM" = "#417839FF", "TLN" = "palegreen3")

pca_plot <- ggplot(data = pca_coords, aes(x=PC1, y=PC2, colour = Site, fill = Site)) + 
  geom_jitter(size = 5, alpha = 0.75) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  # scale_y_continuous(breaks = c(-10,-5,0,5), limits = c(-15,7)) +
  # scale_x_continuous(breaks = c(-5,0,5,10,15), limits = c(-5,15)) +
  scale_y_continuous(breaks = c(-5,0,5,10), limits = c(-7,13)) +
  scale_x_continuous(breaks = c(-5,0,5,10,15,20,25), limits = c(-7,25)) +
labs(x= "PC1 62%", y= "PC2 15%", colour = "Site type", fill = "Site type")
pca_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/SNPs_PCA/pca_S_orbicularis_gl_snp_maf_nPc_140_230_plot.png", pca_plot, width = 8, height = 7)
```

#### Split the SNPs by chromosome for ngsLD

``` bash
#!/bin/bash

#SBATCH -J split_by_chromosome.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=1 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_split_by_chromosome_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_split_by_chromosome_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Path to the input file
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
pos_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos"
beagle_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle"

# Initialize the starting line for the beagle file
start_line=2  # since the first line is the header

# Define the valid chromosome numbers
for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24 mtDNA; do

    # Pad the chromosome number with leading zero if necessary
    output_pos_file=${output_dir}"S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.pos"

    echo "Processing chromosome: ${chr}, output file: ${output_pos_file}"

    # Use awk to filter lines for the current chromosome
    awk -v chr="${chr}" -v out="${output_pos_file}" '$1 == chr { print > out }' $pos_file

    # Count the number of lines in the created file
    line_count=$(wc -l < "$output_pos_file")
    echo "Number of lines in ${output_pos_file}: ${line_count}"
    gzip "$output_pos_file"

    # Extract corresponding lines from the beagle file, including the header
    output_beagle_file=${output_dir}"S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.beagle"
    head -n 1 "$beagle_file" > "$output_beagle_file"

    end_line=$((start_line + line_count - 1))
    sed -n "${start_line},${end_line}p" "$beagle_file" >> "$output_beagle_file"
    gzip "$output_beagle_file"

    # Update the start line for the next iteration
    start_line=$((end_line + 1))

    echo "Created ${output_beagle_file} with the corresponding lines."
done
```

#### Split pop files by chromosome

``` bash
#!/bin/bash

#SBATCH -J split_by_chromosome_pop.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=1 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_split_by_chromosome_pop_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_split_by_chromosome_pop_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Path to the input file
for pop in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN Ocean Mixed Stratified; do
    output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
    pos_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/${pop}_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.pos"
    beagle_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/${pop}_S_orbicularis_gl_snp_maf_nPc_nLD_140_230.beagle"

    # Initialize the starting line for the beagle file
    start_line=2  # since the first line is the header

    # Define the valid chromosome numbers
    for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24 mtDNA; do

        # Pad the chromosome number with leading zero if necessary
        output_pos_file=${output_dir}${pop}"_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.pos"

        echo "Processing chromosome: ${pop} ${chr}, output file: ${output_pos_file}"

        # Use awk to filter lines for the current chromosome
        awk -v chr="${chr}" -v out="${output_pos_file}" '$1 == chr { print > out }' $pos_file

        # Count the number of lines in the created file
        line_count=$(wc -l < "$output_pos_file")
        echo "Number of lines in ${output_pos_file}: ${line_count}"
        gzip "$output_pos_file"

        # Extract corresponding lines from the beagle file, including the header
        output_beagle_file=${output_dir}${pop}"_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.beagle"
        head -n 1 "$beagle_file" > "$output_beagle_file"

        end_line=$((start_line + line_count - 1))
        sed -n "${start_line},${end_line}p" "$beagle_file" >> "$output_beagle_file"
        gzip "$output_beagle_file"

        # Update the start line for the next iteration
        start_line=$((end_line + 1))

        echo "Created ${output_beagle_file} with the corresponding lines."
    done
done
```

## 4.determine_allele_frequencies

### 6.ngsLD_prune

#### Identify linked sites

- ngsLD v1.2.0, samtools v1.13, gcc, gsl
- Determine your number of sites for n_sites
- gzip -dc S_orbicularis.beagle.gz \| wc -l –geno –pos –probs –n_ind the
  number of individuals you have –n_sites the number of sites found in
  your beagle or glf file –max_kb_dist maximum distance between SNPs (in
  Kb) to calculate LD; set to 0(zero) to disable filter –n_threads the
  number of threads to use –out the output file name

``` bash
#!/bin/bash

#SBATCH -J 6.ngsLD_prune_chr*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=20 #number of cores per node
#SBATCH --mem-per-cpu 12G # the amount of memory per core to request in MB.
#SBATCH -t 0-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_6.ngsLD_prune_chr*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_6.ngsLD_prune_chr*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/"
ngsLD="/home/bcarlson4/local/ngsTools/ngsLD/ngsLD"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24; do

        # Calculate the number of sites
        sites=$(zcat "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.pos.gz" | wc -l)

        # Determine linkage groups
        $ngsLD \
          --geno "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.beagle.gz" \
          --pos "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.pos.gz" \
          --probs \
          --n_ind 154 \
          --n_sites 36073 \
          --max_kb_dist 0 \
          --n_threads 10 \
          --out "${output_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps"
done
```

#### Site/Type

``` bash
#!/bin/bash

#SBATCH -J 6.ngsLD_prune_site_type_chr*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=20 #number of cores per node
#SBATCH --mem-per-cpu 12G # the amount of memory per core to request in MB.
#SBATCH -t 0-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_6.ngsLD_prune_site_type_chr*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_6.ngsLD_prune_site_type_chr*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7


# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/"
ngsLD="/home/bcarlson4/local/ngsTools/ngsLD/ngsLD"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for pop in CLM FLK NLK OTM TLN ULN Ocean Mixed; do
    for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24; do

        # Calculate the number of sites
        sites=$(zcat "${input_dir}${pop}_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.pos.gz" | wc -l)

        # Determine linkage groups
        $ngsLD \
          --geno "${input_dir}${pop}_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.beagle.gz" \
          --pos "${input_dir}${pop}_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.pos.gz" \
          --probs \
          --n_ind 12 \
          --n_sites "$sites" \
          --max_kb_dist 0 \
          --n_threads 10 \
          --out "${output_dir}${pop}_S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps"
    done
done
```

#### Fit an LDdecay model to your SNPs

``` bash
#!/bin/bash

#SBATCH -J fit_LDdecay.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_fit_LDdecay_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_fit_LDdecay_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Give the paths to your directories
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/"
ngsLD="/home/bcarlson4/local/ngsTools/ngsLD/scripts/"

for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24; do

    echo -e "site1\tsite2\tdist\tr2_ExpG\tD\tDp\tr2" | cat - "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps" > "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}_filt.snps"
    
    mv "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}_filt.snps" "${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps"

# Load conda, activate R, run R script
module load anaconda3
source activate my-R
ls ${input_dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps | Rscript --vanilla --slave ${ngsLD}fit_LDdecay.R \
    --ld r2 --n_ind 154 \
    --max_kb_dist 500 \
    --fit_level 100 --fit_boot 100 \
    --out ${input_dir}ld_decay_140_230_${chr}.pdf
conda deactivate
done
```

#### Plot LD blocks

``` bash
#!/bin/bash

#SBATCH -J LD_blocks.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_LD_blocks_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_LD_blocks_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Give the paths to your directories
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"
ngsLD="/home/bcarlson4/local/ngsTools/ngsLD/scripts/"

# Plot LD blocks
cat "${input_dir}"S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_chr01.snps | bash ${ngsLD}LD_blocks.sh \
chr01 26674063 56506954

# chr01:26674063  chr01:56506954
```

#### Prune your SNPs to get rid of linked SNPs.

- <https://github.com/fgvieira/prune_graph>

``` bash
#!/bin/bash

#SBATCH -J LD_prune.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 8G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_LD_prune_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_LD_prune_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
module load anaconda3
source activate prune-graph

dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/"
prune_graph="/home/bcarlson4/local/ngsTools/ngsLD/prune_graph/target/release/prune_graph"

for chr in chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr24; do

# Use for SNP files without a header
#$prune_graph --in ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps --weight-field "column_7" --weight-filter "column_7 >= 0.5" --out ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.prunedSNPs

# Use for SNP files with a header
#$prune_graph --in ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.snps --header --weight-field "r2" --weight-filter "r2 >= 0.5" --out ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230_${chr}.prunedSNPs

# --weight-filter "column_3 <= 50000 && column_7 >= 0.5"

done
conda deactivate
```

#### Create your new SNP list that has been LDpruned

- Check to see if you can just use angsd index to create SNP list
- R code comes from here:
  <https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial3_ld_popstructure/markdowns/ld.md>

``` bash
#!/bin/bash

#SBATCH -J create_SNP_list.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_create_SNP_list_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_create_SNP_list_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/"

# Load conda, activate R, run R script
module load anaconda3
source activate my-R
Rscript ${dir}create_SNP_list.R
conda deactiavte 

cat ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_chr*.snps > ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Index SNPs
$angsd sites index ${dir}S_orbicularis_gl_snp_maf_nPc_nLD_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps
```

``` r
# Define the base directory
basedir <- "/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/"

# Define chromosomes, skipping chr23
chromosomes <- c(paste0("chr", sprintf("%02d", 1:22)), "chr24")

# Loop over chromosomes
for (chr in chromosomes) {
  
  # Read the pruned positions for the current chromosome
  pruned_position <- as.integer(gsub(paste0(chr, ":"), "", 
                                     readLines(paste0(basedir, "6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/S_orbicularis_gl_snp_maf_nPc_nLD_140_230_", chr, ".prunedSNPs"))))
  
  # Read the SNP list (assuming the SNP list file contains data for all chromosomes)
  snp_list <- read.table(paste0(basedir, "5.angsd_gl_snp_maf_nPc/S_orbicularis_gl_snp_maf_nPc_140_230/S_orbicularis_gl_snp_maf_nPc_140_230.mafs.gz"), 
                         stringsAsFactors = F, header = T)[,1:4]
  
  # Filter SNPs for the current chromosome and matching positions
  snp_list_chr <- snp_list[snp_list$chromo == chr & snp_list$position %in% pruned_position, ]
  
  # Write the filtered SNP list to a file
  write.table(snp_list_chr, 
              paste0(basedir, "6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_", chr, ".snps"), 
              col.names = F, row.names = F, quote = F, sep = "\t")
  
  # Print progress
  print(paste("Processed", chr, "and wrote output to file"))
}
```

## 4.determine_allele_frequencies

### 7.angsd_gl_snp_maf_nPc_nLDp

#### Generate genotype likelihoods and SNPs based on linkage disequilibrium

- angsd v0.94, samtools v1.13, gcc, gsl
- Make sure the ‘site’ file contains at least two columns (chromosome
  and position) plus you can include 2 additional columns (major and
  minor).
- Need to index ‘sites’ file before hand, angsd sites index filename

``` bash
#!/bin/bash

#SBATCH -J 7.angsd_gl_snp_maf_nPc_nLDp.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_7.angsd_gl_snp_maf_nPc_nLDp_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_7.angsd_gl_snp_maf_nPc_nLDp_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/gatk_indel_realigner_bam.list"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/6.ngsLD_prune/S_orbicularis_gl_snp_maf_nPc_nLD_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps "
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

      # Call variant with informed linkage
      $angsd -b $file_list -ref $ref -anc $ref \
        -out ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 140 -setMinDepth 230 -doCounts 1 \
        -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -GL 1 -doGlf 2 -doMajorMinor 3 \
        -doPost 2 -doIBS 1 -doCov 1 -makeMatrix 1 -nthreads 10 -sites ${SNPs} \
        -doHWE 1  -HWE_pval_F 1 -dobcf 1 --ignore-RG 0 -dogeno 1

zcat ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.mafs.gz | cut -f1-4 | tail -n+2 > ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps

$angsd sites index ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps
```

- Download file

``` bash
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.covMat /Users/bailey/Documents/research/S_orbicularis/DNA/results/gl_SNP/
```

#### Run a PCA based on your output SNPs

``` r
library(tidyverse)
library(dplyr)

# Define basedir
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/results/gl_SNP/"

# Load the covariance matrix
cov <- as.matrix(read.table(paste0(basedir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.covMat"), header = F))

# We will also add a column with population assingments
pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/S_orbicularis_154_samples.csv")

# Run PCA
pca <- prcomp(cov, center = TRUE, scale. = TRUE)
summary(pca)
    
# Biplot
biplot(pca)

# Scree plot
plot(pca)

# Pull out principal component data frame
pca_coords <- pca$x

pca_coords <- cbind(pca_coords, pop)

# Loadings of variables on PC1
pca_loadings <- as.data.frame(pca$rotation)

# Step 1: Select only 12 random rows with "Ocean"
ocean_sample <- pca_coords %>%
  filter(Site == "Ocean") %>%
  sample_n(48)

# Step 2: Filter out all rows with "Ocean" from the original dataframe
pca_coords_no_ocean <- pca_coords %>%
  filter(Site != "Ocean")

# Step 3: Combine the rows without "Ocean" with the 12 sampled "Ocean" rows
pca_coords <- bind_rows(pca_coords_no_ocean, ocean_sample)

# Type
pca_coords$Type <- factor(pca_coords$Type, levels = c("Ocean", "Mixed", "Stratified"))

custom_colors <- c("Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

# Population
#pca_coords$Population <- factor(pca_coords$Site, levels = c("BAL", "CCM", "CCN", "CRF", "LMU", "NCO", "OCK", "OCM", "TKC", "FLK", "ULN", "CLM", "NLK", "OTM", "TLN"))

#custom_colors <- c("BAL" = "#BE3428FF", "CCM" = "#8E2322FF", "CCN" = "#D9565CFF", "CRF" = "#FF3200FF", "LMU" = "#E9A17CFF", "NCO" = "#881C00FF", "OCK" = "#D6604DFF", "OCM" = "#ED3F39FF", "TKC" = "#803233FF", "FLK" = "#69D2E7FF", "ULN" = "#639CA4FF", "CLM" = "#204035FF", "NLK" = "#6D8325FF", "OTM" = "#417839FF", "TLN" = "palegreen3")

# Site
pca_coords$Site <- factor(pca_coords$Site, levels = c("NLK", "CLM", "OTM", "TLN", "FLK", "ULN", "Ocean"))

custom_colors <- c("Ocean" = "#EE6363", "FLK" = "#69D2E7FF", "ULN" = "#639CA4FF", "CLM" = "#204035FF", "NLK" = "#6D8325FF", "OTM" = "#417839FF", "TLN" = "palegreen3")

# Separate "Ocean" points from the other points
ocean_points <- pca_coords %>% filter(Site == "Ocean")
other_points <- pca_coords %>% filter(Site != "Ocean")

pca_plot <- ggplot(data = pca_coords, aes(x=PC1, y=PC2, colour = Site, shape = Type)) + 
  geom_jitter(data = ocean_points, aes(x = PC1, y = PC2, color = Site), size = 5, alpha = 0.75, width = 0.5, height = 0.5, fill = NA, stroke = 1) +
  geom_jitter(data = other_points, aes(x = PC1, y = PC2, color = Site), size = 5, alpha = 0.75, width = 0.5, height = 0.5, fill = NA, stroke = 1) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(22, 21, 24)) + 
#  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  scale_y_continuous(breaks = c(-10,-5,0,5,10), limits = c(-15,11.5)) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10), limits = c(-15,11.5)) +
guides(color = "none", fill = "none", shape = "none") +
labs(x= "PC1 10%", y= "PC2 8%", colour = "Site type", fill = "Site type")
pca_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/SNPs_PCA/pca_PC1_PC2_plot.png", pca_plot, width = 6, height = 5)

pca_plot <- ggplot(data = pca_coords, aes(x=PC3, y=PC4, colour = Site, shape = Type)) + 
  geom_jitter(data = ocean_points, aes(x = PC3, y = PC4, color = Site), size = 5, alpha = 0.75, width = 0.5, height = 0.5, fill = NA, stroke = 1) +
  geom_jitter(data = other_points, aes(x = PC3, y = PC4, color = Site), size = 5, alpha = 0.75, width = 0.5, height = 0.5, fill = NA, stroke = 1) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c(22, 21, 24)) + 
#  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  scale_y_continuous(breaks = c(-10,-5,0,5,10), limits = c(-15,11.5)) +
  scale_x_continuous(breaks = c(-10,-5,0,5,10), limits = c(-15,11.5)) +
guides(color = "none", fill = "none", shape = "none") +
labs(x= "PC3 8%", y= "PC4 6%", colour = "Site type", fill = "Site type")
pca_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/SNPs_PCA/pca_PC3_PC4_plot.png", pca_plot, width = 6, height = 5)
```

## 4.determine_allele_frequencies

### 8.angsd_saf_maf

#### Sample allele frequency

- Takes 2 days and 5 hours
- angsd v0.94, samtools v1.13, gcc, gsl -doSaf 1: perform multisample GL
  estimation \#### Site

``` bash
#!/bin/bash

#SBATCH -J 8.angsd_saf_maf_site.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:01:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_8.angsd_saf_maf_site_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_8.angsd_saf_maf_site_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/8.angsd_saf_maf/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Call allele frequencies, genotype likelihoods, and SNPs
for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN;
do
        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_site_bam.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 140 -setMinDepth 230 -doCounts 1 \
            -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -GL 1 -doGlf 2 -doMajorMinor 3 \
            -doPost 2 -doIBS 1 -doCov 1 -makeMatrix 1 -nthreads 10 \
            -doHWE 1  -HWE_pval_F 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -doSaf 1

zcat ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.mafs.gz | cut -f1-4 | tail -n+2 > ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps

$angsd sites index ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps
done
```

#### Type

``` bash
#!/bin/bash

#SBATCH -J 8.angsd_saf_maf_type.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:01:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_gl_snp_maf_nPc_type_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_gl_snp_maf_nPc_type_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/8.angsd_saf_maf/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

for POP in Ocean Mixed Stratified; do

        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_type_bam.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20 -minInd 140 -setMinDepth 230 -doCounts 1 \
            -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -GL 1 -doGlf 2 -doMajorMinor 3 \
            -doPost 2 -doIBS 1 -doCov 1 -makeMatrix 1 -nthreads 10 \
            -doHWE 1  -HWE_pval_F 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -doSaf 1

zcat ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.mafs.gz | cut -f1-4 | tail -n+2 > ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps

$angsd sites index ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps
done
```

#### Subsample

``` bash
#!/bin/bash

#SBATCH -J 8.angsd_saf_maf_type_subsampled.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:01:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_8.angsd_saf_maf_type_subsampled_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_8.angsd_saf_maf_type_subsampled_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/4.gatk_indel_realigner/"
SNPs="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.snps"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/8.angsd_saf_maf/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.reference/ref.fna"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Call ANGSD
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

for POP in Ocean Mixed; do

# Subsample; don't subsample again if you already have
# shuf -n 12 ${input_dir}$POP/$POP'_angsd_type_bam.list' > ${input_dir}$POP/$POP'_angsd_type_bam_subsampled.list'

        echo $POP
        $angsd -b ${input_dir}$POP/$POP'_angsd_type_bam_subsampled.list' \
            -ref $ref -anc $ref -sites $SNPs \
            -out ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_subsampled' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -minQ 20  -doCounts 1 -GL 1 -doGlf 2 -doMajorMinor 3 -doPost 2 \
            -doMaf 1 -minMaf 0.05 -SNP_pval 1e-6 -doHWE 1 -doIBS 1 -doCov 1 -makeMatrix 1 -P 10

zcat ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_subsampled.mafs.gz | cut -f1-4 | tail -n+2 > ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_subsampled.snps

$angsd sites index ${output_dir}$POP_S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_subsampled.snps
done
```

## 4.determine_allele_frequencies

### 9.winsfs_maf

#### Site frequency spectrum

- Takes 6 hours
- angsd v0.94, samtools v1.13, gcc, gsl
- Fst weighted is sum(a)/sum(a+b), while unweighted is
  average(a/(a+b))\] \#### Site/Type/Subsample

``` bash
#!/bin/bash

#SBATCH -J 9.winsfs_maf_site_type.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_9.winsfs_maf_site_type_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_9.winsfs_maf_site_type_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/8.angsd_saf_maf/S_orbicularis_gl_snp_maf_nPc_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/9.winsfs_maf/S_orbicularis_gl_snp_maf_nPc_140_230/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Run a for loop to calculate site frequency spectrums (SFS)
for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN; do
      
      #SFS
      winsfs -v ${input_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' -t 10 > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' > ${output_dir}$POP'_fold_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
      mv ${output_dir}$POP'_fold_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
done

populations=("BAL" "CCM" "CCN" "CLM" "CRF" "FLK" "LMU" "NCO" "NLK" "OCK" "OCM" "OTM" "TKC" "TLN" "ULN")

# Loop over populations. Arrays start at 0, which is why j will still be less than 15 when it gets to ULN
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # SFS
      winsfs ${input_dir}${pop1}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${input_dir}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' -v -t 10 \
        > ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      tail -n +2 ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' > ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
      mv ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
    done
done

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in Ocean Mixed Stratified; do
      
      #SFS
      winsfs -v ${input_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' -t 10 > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' > ${output_dir}$POP'_fold_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
      mv ${output_dir}$POP'_fold_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
done

populations=("Ocean" "Mixed" "Stratified")

# Loop over populations, arrays start at 0 which is why j will still be less than 3 when it gets to Ocean
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # SFS
      winsfs ${input_dir}${pop1}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${input_dir}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' -v -t 10 \
        > ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      tail -n +2 ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' > ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
      mv ${output_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' ${output_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
    done
done

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in Ocean Mixed; do

      #SFS
      winsfs -v ${input_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.saf.idx' -t 10 > ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs'

      winsfs view -v -f ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs'

      tail -n +2 ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs' > ${output_dir}$POP'_fold_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs'
      mv ${output_dir}$POP'_fold_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs' ${output_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs'
done

# Ocean vs. Marine Lakes
# Loop over populations.
for POP in CLM NLK OTM TLN FLK ULN; do

      # SFS
      winsfs ${input_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${input_dir}'Ocean_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.saf.idx' -v -t 10 \
        > ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' > ${output_dir}$POP'Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
      mv ${output_dir}$POP'Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' ${output_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
done

# Mixed vs. Stratified Lakes
# Loop over populations.
for POP in CLM NLK OTM TLN; do

      # SFS
      winsfs ${input_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${input_dir}'Mixed_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.saf.idx' -v -t 10 \
        > ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      winsfs view -v -f ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'

      tail -n +2 ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' > ${output_dir}$POP'Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
      mv ${output_dir}$POP'Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' ${output_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs'
done
```

## 5.analyze

### 1.angsd_theta_fst

#### 1.angsd_theta, 2.angsd_fst,

### Fst, Tajima’s D

- Takes days
- angsd v0.94, samtools v1.13, gcc, gsl
- Fst weighted is sum(a)/sum(a+b), while unweighted is
  average(a/(a+b))\] \### Site/Type/Subsample

``` bash
#!/bin/bash

#SBATCH -J 1.angsd_theta_fst.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_1.angsd_theta_fst_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_1.angsd_theta_fst_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
saf_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/3.angsd_saf/S_orbicularis_gl_snp_nPc_140_230/"
sfs_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/4.winsfs/S_orbicularis_gl_snp_nPc_140_230/"
theta_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.analyze/1.angsd_theta_fst/1.angsd_theta/S_orbicularis_gl_snp_nPc_140_230/"
fst_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.analyze/1.angsd_theta_fst/2.angsd_fst/S_orbicularis_gl_snp_nPc_140_230/"
realSFS="/home/bcarlson4/local/ngsTools/angsd/misc/realSFS"
thetaStat="/home/bcarlson4/local/ngsTools/angsd/misc/thetaStat"

# Create the output directory if it doesn't exist
mkdir -p "$theta_dir"
mkdir -p "$fst_dir"

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN; do

      # Thetas
      $realSFS saf2theta ${saf_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        -outname ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230' \
        -sfs ${sfs_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs' -fold 1

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.idx' \
        -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.txt'

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.idx' \
        -win 100000 -step 10000 -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.windows.txt'
done

populations=("BAL" "CCM" "CCN" "CLM" "CRF" "FLK" "LMU" "NCO" "NLK" "OCK" "OCM" "OTM" "TKC" "TLN" "ULN")

# Loop over populations. Arrays start at 0, which is why j will still be less than 15 when it gets to ULN
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # Fst
      $realSFS fst index \
        ${saf_dir}${pop1}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${saf_dir}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
       -sfs ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' \
       -fstout ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
       > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt'
    done
done

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in Ocean Mixed Stratified; do
      $realSFS saf2theta ${saf_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        -outname ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230' \
        -sfs ${sfs_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.sfs' -fold 1

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.idx' \
        -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.txt'

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.idx' \
        -win 100000 -step 10000 -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.thetas.windows.txt'
done

populations=("Ocean" "Mixed" "Stratified")

# Loop over populations, arrays start at 0 which is why j will still be less than 3 when it gets to Ocean
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # SFS
      winsfs ${saf_dir}${pop1}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${saf_dir}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' -v -t 10 \
        > ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      winsfs view -v -f ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      tail -n +2 ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' > ${sfs_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'
      mv ${sfs_dir}${pop1}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs'

      # Fst
      $realSFS fst index \
        ${saf_dir}${pop1}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${saf_dir}${pop2}'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
       -sfs ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.sfs' \
       -fstout ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
       > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt'
    done
done

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in Ocean Mixed; do

      # Theta
      $realSFS saf2theta ${saf_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.saf.idx' \
        -outname ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled' \
        -sfs ${sfs_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.sfs' -fold 1

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.thetas.idx' \
        -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.thetas.txt'

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.thetas.idx' \
        -win 100000 -step 10000 -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230_subsampled.thetas.windows.txt'
done

# Ocean vs. Marine Lakes
# Loop over populations
for POP in CLM NLK OTM TLN FLK ULN; do

      # Fst
      $realSFS fst index \
        ${saf_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${saf_dir}'Ocean_S_orbicularis_gl_snp_nPc_140_230_subsampled.saf.idx' \
       -sfs ${sfs_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.sfs' \
       -fstout ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
       > ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt'
done

# Mixed vs. Stratified Lakes
# Loop over populations
for POP in CLM NLK OTM TLN; do

      # Fst
      $realSFS fst index \
        ${saf_dir}$POP'_S_orbicularis_gl_snp_nPc_140_230.saf.idx' \
        ${saf_dir}'Mixed_S_orbicularis_gl_snp_nPc_140_230_subsampled.saf.idx' \
       -sfs ${sfs_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.sfs' \
       -fstout ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
       > ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt'
done
```

#### Watterson’s theta & Tajima’s D Plots

``` r
library(ggplot2)
library(dplyr)

# Read a text file with space or tab-delimited data
NLK_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/NLK_S_orbicularis_gl_snp_nPc_140_230.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

NLK_data <- NLK_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(NLK_data) <- columns

# Remove "chr" prefix and convert to numeric
NLK_data$chr<- as.factor(gsub("chr", "", NLK_data$chr))

NLK_data$Site <- "Strat1"

# Read a text file with space or tab-delimited data
CLM_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/CLM_S_orbicularis_gl_snp_nPc_140_230.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

CLM_data <- CLM_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(CLM_data) <- columns

# Remove "chr" prefix and convert to numeric
CLM_data$chr<- as.factor(gsub("chr", "", CLM_data$chr))

CLM_data$Site <- "Strat2"

# Read a text file with space or tab-delimited data
OTM_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/OTM_S_orbicularis_gl_snp_nPc_140_230.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

OTM_data <- OTM_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(OTM_data) <- columns

# Remove "chr" prefix and convert to numeric
OTM_data$chr<- as.factor(gsub("chr", "", OTM_data$chr))

OTM_data$Site <- "Strat3"

# Read a text file with space or tab-delimited data
TLN_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/TLN_S_orbicularis_gl_snp_nPc_140_230.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

TLN_data <- TLN_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(TLN_data) <- columns

# Remove "chr" prefix and convert to numeric
TLN_data$chr<- as.factor(gsub("chr", "", TLN_data$chr))

TLN_data$Site <- "Strat4"

# Read a text file with space or tab-delimited data
FLK_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/FLK_S_orbicularis_gl_snp_nPc_140_230.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

FLK_data <- FLK_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(FLK_data) <- columns

# Remove "chr" prefix and convert to numeric
FLK_data$chr<- as.factor(gsub("chr", "", FLK_data$chr))

FLK_data$Site <- "Mix1"

# Read a text file with space or tab-delimited data
ULN_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/ULN_S_orbicularis_gl_snp_nPc_140_230.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

ULN_data <- ULN_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(ULN_data) <- columns

# Remove "chr" prefix and convert to numeric
ULN_data$chr<- as.factor(gsub("chr", "", ULN_data$chr))

ULN_data$Site <- "Mix2"

# Read a text file with space or tab-delimited data
Ocean_data <- read.table("/Users/bailey/Documents/research/S_orbicularis/DNA/results/theta/Ocean_S_orbicularis_gl_snp_nPc_140_230_subsampled.thetas.txt.pestPG", header = F, sep = "\t", stringsAsFactors = FALSE)

Ocean_data <- Ocean_data[,c(1:4,9,14)]

columns <- c("region", "chr", "WinCenter", "tW", "Tajima", "nSites")

colnames(Ocean_data) <- columns

# Remove "chr" prefix and convert to numeric
Ocean_data$chr<- as.factor(gsub("chr", "", Ocean_data$chr))

Ocean_data$Site <- "Ocean"

site_data <- rbind(NLK_data,CLM_data,OTM_data,TLN_data,FLK_data,ULN_data,Ocean_data)

# Site
site_data$Site <- factor(site_data$Site, levels = c("Strat2", "Strat3", "Strat1", "Strat4", "Mix1", "Mix2", "Ocean"))

custom_colors <- c("Strat2" = "#204035FF", "Strat3" = "#417839FF", "Strat1" = "#6D8325FF", "Strat4" = "palegreen3", "Mix1" = "#69D2E7FF", "Mix2" = "#639CA4FF", "Ocean" = "#EE6363")

site_plot <- ggplot(site_data, aes(x = chr, y = tW, color = Site, fill = Site)) +
  geom_point(size = 5, alpha = 0.75) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_blank(), 
    legend.text = element_blank(),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.key = element_blank()) +
  labs(x= "Chromosome", y= expression("Watterson" * "'" * "s" ~ theta ~ (k)), colour = "Site type", fill = "Site type") +
  scale_y_continuous(breaks = seq(0, 20000, by = 5000), labels = function(x) x / 1000) +
  guides(color = "none", fill = "none")
site_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/theta/wattersons_site_chr_plot.png", site_plot, width = 14, height = 5.5)

site_plot <- ggplot(site_data, aes(x = chr, y = Tajima, color = Site, fill = Site)) +
  geom_point(size = 5, alpha = 0.75) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  labs(x= "Chromosome", y= "Tajima's D", colour = "Site type", fill = "Site type") +
  guides(color = "none", fill = "none")
site_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/theta/tajimaD_site_chr_plot.png", site_plot, width = 14, height = 5.5)
```

#### Fst Heatmap

``` r
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)

# Set the working directory to the folder containing your text files
{setwd("/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/")

# List all text files in the directory
files <- list.files(pattern = "*_S_orbicularis_gl_snp_nPc_140_230.fst.txt")

# Create a data frame to store the pairwise comparison values
comparison_data <- data.frame()

# Loop through each file, extract the second value, and store it
for (file in files) {
  # Load the matrix
  matrix_data <- as.matrix(read.table(file, header = FALSE))
  
  # Assuming the second value you need is in the (1,2) position
  value <- matrix_data[1, 2]
  
  # Extract pair names from the filename (customize based on your file naming scheme)
  pair_names <- strsplit(file, "_")[[1]]
  
  # Store the extracted value along with pair names
  comparison_data <- rbind(comparison_data, data.frame(First = pair_names[1], Second = pair_names[2], Value = value))
}}

# write.csv(comparison_data, file = "/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/population_fst.csv")

#Switch sites around to change the structure of the plot

comparison_data_edited <- read.csv("/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/population_fst.csv")

# Plot the heatmap
fst_heatmap <- ggplot(data = comparison_data_edited, 
                      aes(x = factor(First, levels = c("NLK", "CLM", "OTM", "TLN", "FLK", "ULN", "Ocean")),
                          y = factor(Second, levels = c("Ocean", "ULN", "FLK", "TLN", "OTM", "CLM")),
                          fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), color = "white", size = 5) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(0,1), limits = c(0,1)) +
  scale_x_discrete(labels = c("CLM" = "Strat2", "FLK" = "Mix1", "NLK" = "Strat1", 
                              "OTM" = "Strat3", "TLN" = "Strat4", "ULN" = "Mix2"), expand = c(0,0)) +
  scale_y_discrete(labels = c("Ocean" = "Ocean", "ULN" = "Mix2", "TLN" = "Strat4", 
                              "OTM" = "Strat3","FLK" = "Mix1", "CLM" = "Strat2"), expand = c(0,0)) +
  labs(x = "Sites", y = "Sites", fill = expression(F[ST])) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())
fst_heatmap
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/fst_heatmap.png", fst_heatmap, width = 8, height = 6)
```

#### Fst Manhattan Plots

``` r
# Load the library
library(qqman)
library(ggplot2)
library(dplyr)

# Set the working directory to the folder containing your text files
setwd("/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/")

# Before loading in txt file make sure that "Fst" is added to the header for the last column

# Read a text file with space or tab-delimited data
CLM_Ocean_data <- read.table("CLM_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
CLM_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", CLM_Ocean_data$chr))

# Replace negative numbers with 0
CLM_Ocean_data$Fst_adjusted <- ifelse(CLM_Ocean_data$Fst < 0, 0, CLM_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/CLM_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(CLM_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="CLM-Ocean Fst", col = c("#204035FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
NLK_Ocean_data <- read.table("NLK_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
NLK_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", NLK_Ocean_data$chr))

# Replace negative numbers with 0
NLK_Ocean_data$Fst_adjusted <- ifelse(NLK_Ocean_data$Fst < 0, 0, NLK_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/NLK_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(NLK_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="NLK-Ocean Fst", col = c("#6D8325FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
OTM_Ocean_data <- read.table("OTM_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
OTM_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", OTM_Ocean_data$chr))

# Replace negative numbers with 0
OTM_Ocean_data$Fst_adjusted <- ifelse(OTM_Ocean_data$Fst < 0, 0, OTM_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/OTM_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(OTM_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="OTM-Ocean Fst", col = c("#417839FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
TLN_Ocean_data <- read.table("TLN_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
TLN_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", TLN_Ocean_data$chr))

# Replace negative numbers with 0
TLN_Ocean_data$Fst_adjusted <- ifelse(TLN_Ocean_data$Fst < 0, 0, TLN_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/TLN_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(TLN_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="TLN-Ocean Fst", col = c("palegreen3", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
FLK_Ocean_data <- read.table("FLK_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
FLK_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", FLK_Ocean_data$chr))

# Replace negative numbers with 0
FLK_Ocean_data$Fst_adjusted <- ifelse(FLK_Ocean_data$Fst < 0, 0, FLK_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/FLK_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(FLK_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="FLK-Ocean Fst", col = c("#69D2E7FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
ULN_Ocean_data <- read.table("ULN_Ocean_S_orbicularis_gl_snp_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
ULN_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", ULN_Ocean_data$chr))

# Replace negative numbers with 0
ULN_Ocean_data$Fst_adjusted <- ifelse(ULN_Ocean_data$Fst < 0, 0, ULN_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/ULN_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(ULN_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="ULN-Ocean Fst", col = c("#639CA4FF", "#EE6363"))
dev.off()

# Let's highlight SNPs of interest, with a bit of customization on the plot
#manhattan(subset(data, chr_numeric == 3), highlight = snpsOfInterest, xlim = c(200, 500), main = "Chr 3")
```

## 5.analyze

### 2.angsd_admix

### Determine population structure estimating number of populations

- Takes 20 hours
- angsd v0.94, samtools v1.13, gcc, gsl

``` bash
#!/bin/bash

#SBATCH -J 2.angsd_admix.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_2.angsd_admix_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_2.angsd_admix_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.analyze/2.angsd_admix/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/"
NGSadmix="/home/bcarlson4/local/ngsTools/angsd/misc/NGSadmix"
evalAdmix="/home/bcarlson4/local/evalAdmix"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

      # Admixture
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
        $NGSadmix -likes ${input_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.beagle.gz -K ${i} \
        -P 10 -o ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K${i}

        $evalAdmix -beagle ${input_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.beagle.gz \
        -fname ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K${i}.fopt.gz \
        -qname ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K${i}.qopt \
        -P 10 -o ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K${i}.corres.txt
done
```

- Download folder

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.analyze/2.angsd_admix/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/ /Users/bailey/Documents/research/S_orbicularis/DNA/results/admixture/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230
```

#### Determine which K is best supported

``` r
source("/Users/bailey/git/evalAdmix/visFuns.R")

# Define the base directory for the files
base_dir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/results/admixture/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230"
fig_dir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/figures/admixture/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230"

pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/S_orbicularis_154_samples.csv")
pop <- pop$Stratification

# Site
pop <- factor(pop, levels = c("Strat1", "Strat2", "Strat3", "Strat4", "Mix1", "Mix2", "Ocean"))

custom_colors <- c("Strat2" = "#204035FF", "Strat3" = "#417839FF", "Strat1" = "#6D8325FF", "Strat4" = "palegreen3", "Mix1" = "#69D2E7FF", "Mix2" = "#639CA4FF", "Ocean" = "#EE6363")

# row.names(pop) <- pop$Type
# q_file <- as.matrix(read.table(paste0(base_dir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K5.qopt")))
# r_file <- as.matrix(read.table(paste0(base_dir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K5.corres.txt")))

# Loop through K values from 2 to 15
for (k in 2:15) {
  # Define the filenames based on the current K value
  q_file <- paste0(base_dir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K", k, ".qopt")
  r_file <- paste0(base_dir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K", k, ".corres.txt")
  
  # Read in q and r matrices
  q <- as.matrix(read.table(q_file))
  r <- as.matrix(read.table(r_file))
  
  # Order individuals for consistent plotting
  ord <- orderInds(pop = pop, q = q)
  
  # Define file paths for the output images
  q_plot_file <- paste0(fig_dir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K", k, ".qopt.png")
  r_plot_file <- paste0(fig_dir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K", k, ".corres.png")
  
  # Generate the admixture plot and save as PNG
  png(q_plot_file, width = 1200, height = 800, res = 100)
  plotAdmix(q = q, pop = pop, ord = ord)
  dev.off()
  
  # Generate the correlation of residuals plot and save as PNG
  png(r_plot_file, width = 1400, height = 1200, res = 100)
  plotCorRes(cor_mat = r, pop = pop, ord = ord, title = "", max_z = 0.1, min_z = -0.1)
  dev.off()
  
  # Print progress message
  print(paste("Completed processing for K =", k))
}
```

#### Generate specific admixture plot based on best K

``` r
# Load necessary library
library(tidyverse)

# Load the data (replace with actual path if needed)
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/results/admixture/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230/" 
admix <- read_table(paste0(basedir, "S_orbicularis_gl_snp_maf_nPc_nLDp_140_230_admix_K5.qopt"), col_names = FALSE)

# Load the population data
pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/S_orbicularis_154_samples.csv")
pop <- pop[,c(1:4,10:11)]

# Combine pop and admix data
admix.id <- as.data.frame(cbind(pop, admix[,-6]))
names(admix.id) <- c("Sample", "Site", "Type", "Island", "Population", "Stratification", "q1","q2","q3","q4","q5")

# Sort dataframe by Stratification to help with visual grouping
sorted_admix.id <- admix.id[order(admix.id$Stratification), ]

# Reshape the data to long format for ggplot
admix_long <- sorted_admix.id %>%
  select(Sample, Stratification, q1:q5) %>%
  pivot_longer(cols = q1:q5, names_to = "Component", values_to = "Proportion")

# Set the levels for Stratification
admix_long$Stratification <- factor(admix_long$Stratification, levels = c("Strat1", "Strat2", "Strat3", "Strat4", "Mix1", "Mix2", "Ocean"))

# Set Sample as a factor with levels sorted by Stratification, so it plots in the desired order
admix_long$Sample <- factor(admix_long$Sample, levels = sorted_admix.id$Sample)

# Plot using ggplot, with Sample on the x-axis
admix_plot <- ggplot(admix_long, aes(x = Sample, y = Proportion, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6D8325FF", "palegreen3", "#204035FF", "#417839FF", "#EE6363")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    plot.margin = unit(c(1, 0.5, 0, 0.5), "lines"),
    strip.text = element_text(size = 13),  
    strip.placement = "outside" # Place facet labels outside plot area
  ) +
  labs(y = "Admixture Proportion") +
  scale_y_continuous(breaks = c(0,0.5,1), expand = c(0,0)) +
  guides(color = "none", fill = "none") +
  facet_grid(~Stratification, scales = "free_x", space = "free_x", switch = "x") # Use switch to move labels to the bottom
admix_plot
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/admixture/S_orbicularis_gl_snp_maf_nPc_nLDp_140_230final_admix_K5_plot.png", admix_plot, width = 7.25, height = 4)
```

## 5.analyze

### 3.angsd_fst_maf

#### Fst

- Takes days
- angsd v0.94, samtools v1.13, gcc, gsl
- Fst weighted is sum(a)/sum(a+b), while unweighted is
  average(a/(a+b))\] \#### Site/Type/Subsample

``` bash
#!/bin/bash

#SBATCH -J 3.angsd_fst_maf.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=10 #number of cores per node
#SBATCH --mem-per-cpu 16G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.angsd_fst_maf_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.angsd_fst_maf_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Assuming your files are in the current directory
saf_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/8.angsd_saf_maf/S_orbicularis_gl_snp_maf_nPc_140_230/"
sfs_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/9.winsfs_maf/S_orbicularis_gl_snp_maf_nPc_140_230/"
theta_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.analyze/3.angsd_fst_maf/1.angsd_theta/S_orbicularis_gl_snp_maf_nPc_140_230/"
fst_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.analyze/3.angsd_fst_maf/2.angsd_fst/S_orbicularis_gl_snp_maf_nPc_140_230/"
realSFS="/home/bcarlson4/local/ngsTools/angsd/misc/realSFS"
thetaStat="/home/bcarlson4/local/ngsTools/angsd/misc/thetaStat"

# Create the output directory if it doesn't exist
mkdir -p "$theta_dir"
mkdir -p "$fst_dir"

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in BAL CCM CCN CLM CRF FLK LMU NCO NLK OCK OCM OTM TKC TLN ULN; do

      # Theta
      $realSFS saf2theta ${saf_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        -outname ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230' \
        -sfs ${sfs_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' -fold 1

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.idx' \
        -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.txt'

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.idx' \
        -win 100000 -step 10000 -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.windows.txt'
done

populations=("BAL" "CCM" "CCN" "CLM" "CRF" "FLK" "LMU" "NCO" "NLK" "OCK" "OCM" "OTM" "TKC" "TLN" "ULN")

# Loop over populations. Arrays start at 0, which is why j will still be less than 15 when it gets to ULN
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # Fst
      $realSFS fst index \
        ${saf_dir}${pop1}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${saf_dir}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
       -sfs ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' \
       -fstout ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
       > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt'
    done
done

# Run a for loop to calculate thetas
for POP in Ocean Mixed Stratified; do

      # Theta
      $realSFS saf2theta ${saf_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        -outname ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230' \
        -sfs ${sfs_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' -fold 1

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.idx' \
        -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.txt'

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.idx' \
        -win 100000 -step 10000 -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.thetas.windows.txt'
done

populations=("Ocean" "Mixed" "Stratified")

# Loop over populations, arrays start at 0 which is why j will still be less than 3 when it gets to Ocean
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}

      # Fst
      $realSFS fst index \
        ${saf_dir}${pop1}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${saf_dir}${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
       -sfs ${sfs_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' \
       -fstout ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
       > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}${pop1}_${pop2}'_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt'
    done
done

# Run a for loop to calculate site frequency spectrums (SFS) and thetas
for POP in Ocean Mixed; do
      
      # Theta
      $realSFS saf2theta ${saf_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.saf.idx' \
        -outname ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled' \
        -sfs ${sfs_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs' -fold 1

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.thetas.idx' \
        -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.thetas.txt'

      $thetaStat do_stat ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.thetas.idx' \
        -win 100000 -step 10000 -outnames ${theta_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.thetas.windows.txt'
done

# Ocean vs. Marine Lakes
# Loop over populations.
for POP in CLM NLK OTM TLN FLK ULN; do

      # Fst
      $realSFS fst index \
        ${saf_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${saf_dir}'Ocean_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.saf.idx' \
       -sfs ${sfs_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' \
       -fstout ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
       > ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}$POP'_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt'
done

# Mixed vs. Stratified Lakes
# Loop over populations.
for POP in CLM NLK OTM TLN; do

      # Fst
      $realSFS fst index \
        ${saf_dir}$POP'_S_orbicularis_gl_snp_maf_nPc_140_230.saf.idx' \
        ${saf_dir}'Mixed_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.saf.idx' \
       -sfs ${sfs_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.sfs' \
       -fstout ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230' -whichFst 1

      $realSFS fst stats ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
       > ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.fst.txt'

      $realSFS fst stats2 ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.fst.idx' \
        -win 100000 -step 10000 > ${fst_dir}$POP'_Mixed_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt'
done
```

#### Fst_maf Heatmap

``` r
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)

# Set the working directory to the folder containing your text files
{setwd("/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/")

# List all text files in the directory
files <- list.files(pattern = "*_S_orbicularis_gl_snp_maf_nPc_140_230.fst.txt")

# Create a data frame to store the pairwise comparison values
comparison_data <- data.frame()

# Loop through each file, extract the second value, and store it
for (file in files) {
  # Load the matrix
  matrix_data <- as.matrix(read.table(file, header = FALSE))
  
  # Assuming the second value you need is in the (1,2) position
  value <- matrix_data[1, 2]
  
  # Extract pair names from the filename (customize based on your file naming scheme)
  pair_names <- strsplit(file, "_")[[1]]
  
  # Store the extracted value along with pair names
  comparison_data <- rbind(comparison_data, data.frame(First = pair_names[1], Second = pair_names[2], Value = value))
}}

# write.csv(comparison_data, file = "/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/population_fst.csv")

#Switch sites around to change the structure of the plot

comparison_data_edited <- read.csv("/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/population_fst.csv")

# Plot the heatmap
fst_heatmap <- ggplot(data = comparison_data_edited, 
                      aes(x = factor(First, levels = c("NLK", "CLM", "OTM", "TLN", "FLK", "ULN", "Ocean")),
                          y = factor(Second, levels = c("Ocean", "ULN", "FLK", "TLN", "OTM", "CLM")),
                          fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), color = "white", size = 5) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(0,1), limits = c(0,1)) +
  scale_x_discrete(labels = c("CLM" = "Strat2", "FLK" = "Mix1", "NLK" = "Strat1", 
                              "OTM" = "Strat3", "TLN" = "Strat4", "ULN" = "Mix2"), expand = c(0,0)) +
  scale_y_discrete(labels = c("Ocean" = "Ocean", "ULN" = "Mix2", "TLN" = "Strat4", 
                              "OTM" = "Strat3","FLK" = "Mix1", "CLM" = "Strat2"), expand = c(0,0)) +
  labs(x = "Sites", y = "Sites", fill = expression(F[ST])) +
  theme_bw() +
  theme(text = element_text(size = 24), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.text = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())
fst_heatmap
ggsave("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/fst_heatmap.png", fst_heatmap, width = 8, height = 6)
```

#### Fst Manhattan Plot

``` r
# Load the library
library(qqman)
library(ggplot2)
library(dplyr)

# Set the working directory to the folder containing your text files
setwd("/Users/bailey/Documents/research/S_orbicularis/DNA/results/fst/")

# Before loading in txt file make sure that "Fst" is added to the header for the last column

# Read a text file with space or tab-delimited data
CLM_Ocean_data <- read.table("CLM_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
CLM_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", CLM_Ocean_data$chr))

# Replace negative numbers with 0
CLM_Ocean_data$Fst_adjusted <- ifelse(CLM_Ocean_data$Fst < 0, 0, CLM_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/CLM_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(CLM_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="CLM-Ocean Fst", col = c("#204035FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
NLK_Ocean_data <- read.table("NLK_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
NLK_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", NLK_Ocean_data$chr))

# Replace negative numbers with 0
NLK_Ocean_data$Fst_adjusted <- ifelse(NLK_Ocean_data$Fst < 0, 0, NLK_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/NLK_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(NLK_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="NLK-Ocean Fst", col = c("#6D8325FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
OTM_Ocean_data <- read.table("OTM_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
OTM_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", OTM_Ocean_data$chr))

# Replace negative numbers with 0
OTM_Ocean_data$Fst_adjusted <- ifelse(OTM_Ocean_data$Fst < 0, 0, OTM_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/OTM_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(OTM_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="OTM-Ocean Fst", col = c("#417839FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
TLN_Ocean_data <- read.table("TLN_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
TLN_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", TLN_Ocean_data$chr))

# Replace negative numbers with 0
TLN_Ocean_data$Fst_adjusted <- ifelse(TLN_Ocean_data$Fst < 0, 0, TLN_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/TLN_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(TLN_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="TLN-Ocean Fst", col = c("palegreen3", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
FLK_Ocean_data <- read.table("FLK_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
FLK_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", FLK_Ocean_data$chr))

# Replace negative numbers with 0
FLK_Ocean_data$Fst_adjusted <- ifelse(FLK_Ocean_data$Fst < 0, 0, FLK_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/FLK_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(FLK_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="FLK-Ocean Fst", col = c("#69D2E7FF", "#EE6363"))
dev.off()

# Read a text file with space or tab-delimited data
ULN_Ocean_data <- read.table("ULN_Ocean_S_orbicularis_gl_snp_maf_nPc_140_230.fst.windows.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove "chr" prefix and convert to numeric
ULN_Ocean_data$chr_numeric <- as.numeric(gsub("chr", "", ULN_Ocean_data$chr))

# Replace negative numbers with 0
ULN_Ocean_data$Fst_adjusted <- ifelse(ULN_Ocean_data$Fst < 0, 0, ULN_Ocean_data$Fst)

# Make the Manhattan plot on the gwasResults dataset
png("/Users/bailey/Documents/research/S_orbicularis/DNA/figures/fst/ULN_ocean_fst_windows_manhattan_plot.png", width = 1200, height = 800, res = 100)
manhattan(ULN_Ocean_data, chr="chr_numeric", bp="midPos", snp="Nsites", p="Fst_adjusted", logp=FALSE, ylab="ULN-Ocean Fst", col = c("#639CA4FF", "#EE6363"))
dev.off()

# Let's highlight SNPs of interest, with a bit of customization on the plot
#manhattan(subset(data, chr_numeric == 3), highlight = snpsOfInterest, xlim = c(200, 500), main = "Chr 3")
```

#### Heterozygosity

``` r
# Set the working directory to the folder containing your text files
setwd("/Users/bailey/Documents/research/S_orbicularis/DNA/results/sfs/")

# Read a text file with space or tab-delimited data
CLM_data <- scan("CLM_S_orbicularis_gl_snp_maf_nPc_140_230.sfs")
CLM_data[2]/sum(CLM_data)

FLK_data <- scan("FLK_S_orbicularis_gl_snp_maf_nPc_140_230.sfs")
FLK_data[2]/sum(FLK_data)

Ocean_data <- scan("Ocean_S_orbicularis_gl_snp_maf_nPc_140_230_subsampled.sfs")
Ocean_data[2]/sum(Ocean_data)

a<-scan("est.ml")
a[2]/sum(a)
```

## 5.analyze

### 4.ngsF_ngsRelate

#### 

``` bash
#!/bin/bash

#SBATCH -J ngsF.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=20 #number of cores per node
#SBATCH --mem-per-cpu 12G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_ngsF_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_ngsF_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_nPc_nLDp_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/ngsF/"
ngsF="/home/bcarlson4/local/ngsTools/ngsF/ngsF"

zcat ${input_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.beagle.gz | tail -n +2 | perl -an -e 'for($i=3; $i<=$#F; $i++){print(pack("d",($F[$i]==0 ? -inf : log($F[$i]))))}' > ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.glf

$ngsF --n_ind 154 --n_sites 449,419 --glf ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.glf --out ${output_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.ib
```

## 5.analyze

### 7.pcangsd_analyses

#### 

``` bash
#!/bin/bash

#SBATCH -J pcangsd.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --ntasks=20 #number of cores per node
#SBATCH --mem-per-cpu 12G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_pcangsd_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_pcangsd_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

source activate pcangsd

input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/7.angsd_gl_snp_maf_nPc_nLDp/S_orbicularis_gl_nPc_nLDp_140_230/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.determine_allele_frequencies/pcangsd/"

pcangsd -b ${input_dir}S_orbicularis_gl_snp_maf_nPc_nLDp_140_230.beagle.gz -t 20 \
-o ${output_dir}S_orbicularis_gl_nPc_nLDp_140_230/S_orbicularis_gl_nPc_nLDp_140_230 \
--filter ${output_dir}individual.txt --filter_sites ${output_dir}site.txt \
--selection --pcadapt --inbreed_samples --inbreed_sites --admix --admix_K 7 --tree 
```

## Check file qualities

### picard validate summary

#### 

- Takes ? hours

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/aligned/merge/"

# Iterate over the files in the input directory
for file in "${input_dir}"*_merge.bam; do

    # Remove duplicates
    picard ValidateSamFile \
        I="$file" \
        MODE=SUMMARY \
        MAX_OUTPUT=10000
done
```

### picard validate verbose

#### 

- Takes ? hours

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/aligned/merge/"

# Iterate over the files in the input directory
for file in "${input_dir}"*_merge.bam; do

    # Remove duplicates
    picard ValidateSamFile \
        I="$file" \
        IGNORE_WARNINGS=true \
        MODE=VERBOSE \
        MAX_OUTPUT=10000
done
```

## Extra/unused code

### create_mtDNA

#### Create the mitochondrial genome

- Takes -j \[seqid\]  
  -a \[assembly.fasta\]  
  -r \[genbank_reference.gb\]  
  -o \[genetic_code\] -p \[threads\] -m \[memory\]
  <https://www.ncbi.nlm.nih.gov/nuccore/?term=%22mitochondrion%22+and+%22complete+genome%22+and+%22Apogoninae%22>

``` bash
#!/bin/bash

#SBATCH -J mitofinder_mtDNA.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 7-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p dept.les # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_mitofinder_mtDNA_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_mitofinder_mtDNA_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load singularity/3.8.3

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.created_mtDNA/mitofinder_mtDNA.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.created_mtDNA/"
mitofinder="/home/bcarlson4/local/mitofinder/mitofinder_latest.sif"

# Execute your processing script
      singularity run --home $dir \
        $mitofinder \
        -j S_orb \
        -1 fSphaOr1_27320_1%231_R1_paired.fq.gz \
        -2 fSphaOr1_27320_1%231_R2_paired.fq.gz \
        -r apogoninae_ncbi_mt_genomes.gb \
        -o 2 -p 28 -m 240
```

### samtools fixmate

#### Fix the mate pairing information in your sequences, file name samtools_fixmate.sh and corresponds to samtools_fixmate.slurm

- Takes 2 hours
- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/aligned/sort/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/aligned/fixmate/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in "${input_dir}"*_pos.srt.bam; do

    # Fix mate information
    samtools fixmate \
        -O bam,level=1 \
        -m "$file" \
        -@28 \
        "${output_dir}$(basename "${file}" pos.srt.bam)fixmate.bam"
done
```

### samtools sort

#### samtools_sort.sh

- Takes 2.5 hours
- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/1.viewed/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.edit_alignments/2.sorted/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in "${input_dir}"*_hqr.bam; do

      # Sort the aligned reads
      samtools sort \
            -O bam,level=1 \
            -@64 \
            -o "${output_dir}$(basename "${file}" hqr.bam)pos.srt.bam" \
            -T $output_dir \
            "$file"
done
```

### samtools merge

#### Merge lane files together, file name samtools_merge.sh and corresponds to samtools_merge.slurm

- Takes 21 hours
- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/aligned/fixmate/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/aligned/merge/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file_L001 in "${input_dir}"*_L001_fixmate.bam; do
    # Generate file_L002
    file_L002="${file_L001/_L001_/_L002_}"
    # Generate file_L003
    file_L003="${file_L001/_L001_/_L003_}"
    # Generate file_L004
    file_L004="${file_L001/_L001_/_L004_}"
    
      # Fix mate information
      samtools merge \
        "${output_dir}$(basename "${file_L001}" L001_fixmate.bam)merge.bam" \
        "$file_L001" \
        "$file_L002" \
        "$file_L003" \
        "$file_L004" \
        -@28
done
```

### gatk br bqsr

#### 

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/clean_align/libraries/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/clean_align/recalibrate/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in "${input_dir}"*_library.bam; do

    # 
    gatk BaseRecalibrator \
        -I "$file" \
        -R ref.fa \
        -O "${output_dir}$(basename "${file}" library.bam)recal_report.table"

    #
    gatk ApplyBQSR \
        -I "$file" \
        -O "${output_dir}$(basename "${file}" library.bam)recalibrated.bam" \
        -bqsr "${output_dir}$(basename "${file}" library.bam)recal_report.table"

done
```
