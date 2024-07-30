S. orbicularis lcwgs processing
================
Bailey Carlson

#### R Markdown

## Quick notes/code to remember

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
```

## 1.load_raw_sequences

### Trim file names

### Takes less than 10 minutes

``` bash
#!/bin/bash

#SBATCH -J load_raw_sequences.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 128G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:10:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_load_raw_sequences_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_load_raw_sequences_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.raw_sequences/"

tar -xf "$input_dir"/L001/fastq_23203Daw_N23157_L001.tar -C "$input_dir"/L001/
tar -xf "$input_dir"/L002/fastq_23203Daw_N23157_L002.tar -C "$input_dir"/L002/
tar -xf "$input_dir"/L003/fastq_23203Daw_N23157_L003.tar -C "$input_dir"/L003/
tar -xf "$input_dir"/L004/fastq_23203Daw_N23157_L004.tar -C "$input_dir"/L004/
```

#### Use mv 23203Daw_R23157_M0D0\* ../. to move all the sequence files up into the main lane folder L00\*

## 2.rename_sequences

### Trim file names

### Takes less than 10 minutes

- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J rename_sequences_L001.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 24G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_rename_sequences_L001_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_rename_sequences_L001_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/1.raw_sequences/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.rename_sequences/L00*/"

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

## 3.trimmomatic_trim_filter

### Trim sequences according to quality and adapters

### Takes 16.5 hours

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
- I used the TruSeq3-PE-2.fa for adapters found on the
  <https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE-2.fa>
- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J trimmomatic_trim_filter_L00*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 0-16:30:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_trimmomatic_trim_filter_L00*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_trimmomatic_trim_filter_L00*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load trimmomatic/0.39

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.trimmomatic_trim_filter/trimmomatic_trim_filter_L00*.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/2.rename_sequences/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.trimmomatic_trim_filter/L00*/"

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

#### Before working with the reference nuclear genome (GCF_902148855.1) all single digit chromosome numbers need to be changed to have 0 before the number (i.e. chr01). Additionally, you have to add in the complete mtDNA genome (AP018927.1) to the reference nuclear genome fasta file.

## 4.index_reference

### Index your reference genome

### Takes less than 7 minutes

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
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
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

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/index_reference.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
ref_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/"

      # Indexing the reference sequence
        bwa-mem2 index \
          -p ${ref_dir}ref \
          ${ref_dir}ref.fna

      gatk CreateSequenceDictionary \
          -R ${ref_dir}ref.fna

        samtools faidx \
          ${ref_dir}ref.fna
```

## 5.bwamem2_alignment

### Align your sequences to your reference genome

### Takes 3 days, less time for certain lanes

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
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_bwamem2_alignment_L00*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_bwamem2_alignment_L00*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load bwa-mem2/2.2.1

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.bwamem2_alignment/bwamem2_alignment_L00*.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/3.trimmomatic_trim_filter/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.bwamem2_alignment/L00*/"
ref_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/"

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

## 6.edit_alignments

### 1.samtools_view_sort

### Transform SAM to BAM and sort sequences

### Takes 6 hours

- samtools v1.13
- samtools manual <https://www.htslib.org/doc/samtools.html>
- Replace \* with lane number

``` bash
#!/bin/bash

#SBATCH -J 1.samtools_view_sort_L00*.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 0-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_1.samtools_view_sort_L00*_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_1.samtools_view_sort_L00*_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/1.samtools_view_sort/1.samtools_view_sort_L00*.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/5.bwamem2_alignment/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/1.samtools_view_sort/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in "${input_dir}"*_alignment.sam; do

      # Filter out reads with low mapping quality
      samtools view \
        -h \
        -@64 \
        -q 20 \
        -b "$file" |
      samtools view -buS - |
      samtools sort \
        -O bam,level=1 \
        -@64 \
        -o "${output_dir}$(basename "${file}" _bwamem2_alignment.sam)view_sort.bam" \
        -T $output_dir
done
```

## x.check_files

### samtools_flagstat.sh

### Takes 4 hour

- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/1.samtools_view_sort/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.checked_files/flagstat/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in ${input_dir}*_view_sort.bam; do
    samtools flagstat \
        -O tsv $file > ${output_dir}$file.flagstat.out
done
```

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.checked_files/flagstat /Users/bailey/Documents/research/S_orbicularis/DNA/analyses/
```

## 6.edit_alignment

### 2.picard_markdups

### Mark duplicate sequences

### Takes 1 day

- picard v2.26.2
- Input lane files
- Output library file (required)
- Metrics file (required)

``` bash
#!/bin/bash

#SBATCH -J 2.picard_markdups.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 1-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_2.picard_markdups_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_2.picard_markdups_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load  module
module load picard/2.26.2

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/2.picard_markdups/2.picard_markdups.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/1.samtools_view_sort/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/2.picard_markdups/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

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
        -O "${output_dir}$(basename "$file_L001" L001_pos.srt.bam)markdups.bam" \
          -M "${output_dir}$(basename "$file_L001" L001_pos.srt.bam)markdups_metrics.txt" \
        -VALIDATION_STRINGENCY SILENT \
        -REMOVE_DUPLICATES true
done
```

## 6.edit_alignment

### 3.bamutils_clip

### Clip overlapping sequences

### Takes 13 hours

- bamutil v1.0.15

``` bash
#!/bin/bash

#SBATCH -J 3.bamutils_clip.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 0-13:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p medium # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.bamutils_clip_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.bamutils_clip_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load bamutil/1.0.15

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/3.bamutils_clip/3.bamutils_clip.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/2.picard_markdups/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over the files in the input directory
for file in "${input_dir}"*_markdups.bam; do

      # Clip overlapping ends of read pairs
      bam clipOverlap \
            --in "$file" \
            --out "${output_dir}$(basename "${file}" markdups.bam)clip.bam" \
            --stats
done
```

## 6.edit_alignment

### 4a.gatk_target_creator

### Identify where realignment around indels will occur

### Takes 1 day and 6 hours

- gatk v3.7, java v1.8, samtools v1.13
- <https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md>

``` bash
#!/bin/bash

#SBATCH -J 4.gatk_target_creator.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 1-06:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_4.gatk_target_creator_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_4.gatk_target_creator_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/4.gatk_indel_realigner/4.gatk_target_creator.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/4.gatk_indel_realigner/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/ref.fna"
gatk="/home/bcarlson4/local/gatk3.7.0/GenomeAnalysisTK.jar"
java="/home/bcarlson4/local/jdk1.8.0_121/bin/java"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Make the list file
rm -r /home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list
touch ${input_dir}"bamutils_clip_bam.list"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list"

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

## 6.edit_alignment

### 4.gatk_indel_realigner

#### 4b.gatk_indel_realigner

#### Realign sequences around indels

#### Takes 6 days

- gatk v3.7, java v1.8
- <https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md>

``` bash
#!/bin/bash

#SBATCH -J 5.gatk_indel_realigner.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 6-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p dept.les # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.gatk_indel_realigner_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.gatk_indel_realigner_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/4.gatk_indel_realigner/4b.gatk_indel_realigner.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/4.gatk_indel_realigner/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/ref.fna"
gatk="/home/bcarlson4/local/gatk3.7.0/GenomeAnalysisTK.jar"
java="/home/bcarlson4/local/jdk1.8.0_121/bin/java"

# Use the list file
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list"

cd $output_dir

      ## Run the indel realigner tool
      $java -Xmx40g -jar $gatk \
        -T IndelRealigner \
        -R $ref \
        -I $file_list \
        -targetIntervals ${output_dir}indel_realigner.intervals \
        --consensusDeterminationModel USE_READS  \
        --nWayOut _indel_realigner.bam
```

``` bash
#!/bin/bash

#SBATCH -J scp_files.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 0-01:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p short # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_scp_files_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_scp_files_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/5.angsd_saf/"

mkdir ${dir}BAL ${dir}CCM ${dir}CCN ${dir}CLM ${dir}CRF ${dir}FLK ${dir}LMU ${dir}NL1 ${dir}NL2 ${dir}NLK ${dir}OCK ${dir}OCM ${dir}OTM ${dir}TKC ${dir}TLN ${dir}ULN

scp M0D039797C* M0D039798D* M0D039920V* M0D039921W* M0D039922X* M0D039923Y* M0D039924Z* M0D039925A* ${dir}BAL
scp M0D040019Q* M0D040020R* M0D040021S* M0D040022T* M0D040023U* M0D040024V* M0D040025W* M0D040026X* M0D040028Z* M0D040031C* ${dir}CCM
scp M0D039950Z* M0D039951A* M0D039952B* M0D039953C* M0D039954D* M0D039955E* M0D039956F* M0D039957G* ${dir}CCN
scp M0D039846Z* M0D039847A* M0D039848B* M0D039852F* M0D039855I* M0D039856J* M0D039970T* M0D040011I* M0D040012J* M0D040014L* M0D040017O* M0D040018P* ${dir}CLM
scp M0D039800F* M0D039859M* M0D039860N* M0D039861O* M0D039862P* M0D039913O* M0D039914P* M0D039915Q* M0D039916R* M0D039917S* M0D039937M* M0D039872Z* ${dir}CRF
scp M0D039864R* M0D039865S* M0D039866T* M0D039867U* M0D039868V* M0D039926B* M0D039927C* M0D039928D* M0D039929E* M0D039931G* M0D039932H* M0D039933I* ${dir}FLK
scp M0D039958H* M0D039959I* M0D039960J* M0D039961K* M0D039962L* M0D039963M* M0D039964N* M0D039965O* M0D039966P* M0D039967Q* ${dir}LMU
scp M0D040047S* M0D040048T* M0D040049U* M0D040050V* M0D040051W* ${dir}NL1
scp M0D040052X* M0D040053Y* ${dir}NL2
scp M0D039884L* M0D039886N* M0D039888P* M0D039889Q* M0D039891S* M0D039892T* M0D039894V* M0D039895W* M0D039898Z* M0D039900B* M0D039904F* M0D039905G* ${dir}NLK
scp M0D039858L* M0D039907I* M0D039908J* M0D039909K* M0D039910L* M0D039911M* M0D039912N* M0D039934J* M0D039935K* M0D039936L* ${dir}OCK
scp M0D039878F* M0D039879G* M0D039880H* M0D039881I* M0D040034F* M0D040035G* M0D040037I* M0D040038J* M0D040039K* M0D040041M* M0D040042N* M0D040043O* ${dir}OCM
scp M0D039802H* M0D039803I* M0D039804J* M0D039805K* M0D039810P* M0D039811Q* M0D039812R* M0D039813S* M0D039820Z* M0D039873A* M0D039875C* M0D039877E* ${dir}OTM
scp M0D039944T* M0D039946V* M0D039947W* M0D039948X* M0D039949Y* ${dir}TKC
scp M0D039828H* M0D039829I* M0D039830J* M0D039831K* M0D039832L* M0D039833M* M0D039835O* M0D039837Q* M0D039838R* M0D039841U* M0D039842V* M0D039844X* ${dir}TLN
scp M0D039918T* M0D039919U* M0D039938N* M0D039939O* M0D039940P* M0D039941Q* M0D039943S* M0D040027Y* M0D040029A* M0D040030B* M0D040033E* M0D040045Q* ${dir}ULN

mkdir ${dir}Ocean ${dir}Meromictic ${dir}Holomictic

scp M0D039944T* M0D039946V* M0D039947W* M0D039948X* M0D039949Y* M0D039878F* M0D039879G* M0D039880H* M0D039881I* M0D040034F* M0D040035G* M0D040037I* M0D040038J* M0D040039K* M0D040041M* M0D040042N* M0D040043O* M0D039858L* M0D039907I* M0D039908J* M0D039909K* M0D039910L* M0D039911M* M0D039912N* M0D039934J* M0D039935K* M0D039936L* M0D040052X* M0D040053Y* M0D040047S* M0D040048T* M0D040049U* M0D040050V* M0D040051W* M0D039958H* M0D039959I* M0D039960J* M0D039961K* M0D039962L* M0D039963M* M0D039964N* M0D039965O* M0D039966P* M0D039967Q* M0D039800F* M0D039859M* M0D039860N* M0D039861O* M0D039862P* M0D039913O* M0D039914P* M0D039915Q* M0D039916R* M0D039917S* M0D039937M* M0D039872Z* M0D039950Z* M0D039951A* M0D039952B* M0D039953C* M0D039954D* M0D039955E* M0D039956F* M0D039957G* M0D040019Q* M0D040020R* M0D040021S* M0D040022T* M0D040023U* M0D040024V* M0D040025W* M0D040026X* M0D040028Z* M0D040031C* M0D039797C* M0D039798D* M0D039920V* M0D039921W* M0D039922X* M0D039923Y* M0D039924Z* M0D039925A* ${dir}Ocean
scp M0D039828H* M0D039829I* M0D039830J* M0D039831K* M0D039832L* M0D039833M* M0D039835O* M0D039837Q* M0D039838R* M0D039841U* M0D039842V* M0D039844X* M0D039802H* M0D039803I* M0D039804J* M0D039805K* M0D039810P* M0D039811Q* M0D039812R* M0D039813S* M0D039820Z* M0D039873A* M0D039875C* M0D039877E* M0D039884L* M0D039886N* M0D039888P* M0D039889Q* M0D039891S* M0D039892T* M0D039894V* M0D039895W* M0D039898Z* M0D039900B* M0D039904F* M0D039905G* M0D039846Z* M0D039847A* M0D039848B* M0D039852F* M0D039855I* M0D039856J* M0D039970T* M0D040011I* M0D040012J* M0D040014L* M0D040017O* M0D040018P* ${dir}Meromictic
scp M0D039918T* M0D039919U* M0D039938N* M0D039939O* M0D039940P* M0D039941Q* M0D039943S* M0D040027Y* M0D040029A* M0D040030B* M0D040033E* M0D040045Q* M0D039864R* M0D039865S* M0D039866T* M0D039867U* M0D039868V* M0D039926B* M0D039927C* M0D039928D* M0D039929E* M0D039931G* M0D039932H* M0D039933I* ${dir}Holomictic
```

## x.check_files

### samtools_coverage.sh

### Takes 8 hours

- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Use the list file
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list"

# Determine the average coverage
for file in `cat $file_list`; do
        samtools coverage -o $file.coverage.out $file
done

cd /home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/
mv *clip.bam.coverage.out ../../x.checked_files/coverage/
```

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.checked_files/coverage .
```

## x.check_files

### samtools_depth.sh

### Takes 9 hours

- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Use the list file
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/bamutils_clip_bam.list"

for file in `cat $file_list`; do
    ## Count per position depth per sample
    samtools depth -aa $file | cut -f 3 | gzip > $file.depth.gz
done

cd /home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/3.bamutils_clip/
mv *clip.bam.depth.gz ../../x.checked_files/depth/
```

``` bash
scp -r bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/x.checked_files/depth /Users/bailey/Documents/research/S_orbicularis/DNA/analyses/
```

``` r
## install.packages("tidyverse") ## Install tidyverse if you don't have it installed yet
library(tidyverse)

basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/analyses/"
bam_list <- read_lines(paste0(basedir, "sample_list.txt"))

for (i in 1:length(bam_list)){
  bamfile = bam_list[i]
  # Compute depth stats
  depth <- read_tsv(paste0(basedir, "/depth/", bamfile, "clipped_realigned.bam.depth.gz"), col_names = F)$X1
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

output %>%
  mutate(across(where(is.numeric), round, 3))

# Plot the depth distribution (this may take a few minutes to run)
tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
  ggplot(aes(x = position, y = total_depth)) +
  geom_point(size = 0.1)

# Total depth per site across all individuals 
total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)
total_depth_summary %>%
  ggplot(aes(x = total_depth, y = n)) +
  geom_point()
total_depth_summary %>%
  ggplot(aes(x = total_depth, y = n)) +
  geom_point() +
  coord_cartesian(xlim=c(NA, 20))
total_presence_summary %>%
  ggplot(aes(x = total_presence, y = n)) +
  geom_col()
```

## 7.analyze

### 1.angsd_gl_SNP

### Takes 2 days

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
  minimum mapping quality score reads must have to be kept -minQ the
  minimum quality score reads must have to be kept -minInd only keep
  sites with a minimum number of individuals -setMinDepth filter out
  sites with less than this number of reads -setMaxDepth filter out
  sites with more than this number of reads -doCounts calculate various
  counts statistics -doMaf estimate allele frequencies -minMaf remove
  sites with MAF below the specified number -SNP_pval remove sites with
  a pvalue larger than the specified number -GL estimate genotype
  likelihoods, we chose to use the GATK method -doGlf create likelihood
  file, we chose a beagle likelihood file -doMajorMinor specify how to
  assign the major and minor alleles, we used the reference allele as
  major -doPost calculate the posterior probability of genotypes, we
  used a uniform prior -nthreads the number of threads to use for the
  processing, max is 8

``` bash
#!/bin/bash

#SBATCH -J 1.angsd_gl_SNP.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_1.angsd_gl_SNP_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_1.angsd_gl_SNP_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/1.angsd_gl_SNP/1.angsd_gl_SNP.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/1.angsd_gl_SNP/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Make the list file
rm -r /home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/1.angsd_gl_SNP/gatk_indel_realigner_bam.list
touch ${input_dir}"gatk_indel_realigner_bam.list"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/1.angsd_gl_SNP/gatk_indelrealigner_bam.list"

# Iterate over the files in the input directory
for sample in "${input_dir}"*indel_realign.bam; do
    # Append the full path to the file list
    echo $sample >> $file_list
done

      # Call allele frequencies, genotype likelihoods, and SNPs
      $angsd -b $file_list -ref $ref -out ${output_dir}S_orbicularis_gl_SNP \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 -GL 2 -doGlf 2 -doMajorMinor 4 \
        -doPost 2 -doIBS 1 -doCov 1 -makeMatrix 1 -nthreads 24

## Remove header from beagle file output
zcat ${output_dir}S_orbicularis_gl_SNP.beagle.gz | cut -f 4- | gzip  > ${output_dir}S_orbicularis_gl_SNP_LD.beagle.gz

gunzip -c ${output_dir}S_orbicularis_gl_SNP_LD.beagle.gz > ${output_dir}S_orbicularis_gl_SNP_LD.beagle
wc -l ${output_dir}S_orbicularis_gl_SNP_LD.beagle
awk '{print NF}' S_orbicularis_gl_SNP_LD.beagle | sort -nu | tail -n 1

## Remove header from maf file output and replace semicolon with underscore
zcat ${output_dir}S_orbicularis_gl_SNP.mafs.gz | cut -f 1,2 | sed 's/:/_/g'| sed '1d' | gzip > ${output_dir}S_orbicularis_gl_SNP_LD.pos.gz

gunzip -c ${output_dir}S_orbicularis_gl_SNP_LD.pos.gz > ${output_dir}S_orbicularis_gl_SNP_LD.pos
wc -l ${output_dir}S_orbicularis_gl_SNP_LD.pos
```

``` bash
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/1.angsd_gl_SNP/S_orbicularis_angsd_gl_SNP.covMat /Users/bailey/Documents/research/S_orbicularis/DNA/analyses/gl_SNP/
```

``` r
library(tidyverse)
#Define basedir
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/analyses/gl_SNP/"

#Load the covariance matrix
cov <- as.matrix(read.table(paste0(basedir, "S_orbicularis_angsd_gl.covMat"), header = F))

#We will also add a column with population assingments
pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/DNA/data/S_orbicularis_samples.csv")
     
pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors = pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments

pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = Site)) + geom_point()

ggsave(filename = paste0(basedir, "pca_allSNPs_plot.pdf"), plot = pca) #change file path if data on your own computer

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
```

``` bash
#!/bin/bash

#SBATCH -J split_by_chromosome.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_split_by_chromosome_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_split_by_chromosome_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Path to the input file
pos_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/1.angsd_gl_SNP/S_orbicularis_gl_SNP_LD.pos"
beagle_file="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/1.angsd_gl_SNP/S_orbicularis_gl_SNP_LD.beagle"

# Define the valid chromosome numbers
valid_chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 24)

# Initialize the starting line for the beagle file
start_line=2  # since the first line is the header

for chr in "${valid_chromosomes[@]}"; do
    # Pad the chromosome number with leading zero if necessary
    padded_chr=$(printf "%02d" $chr)
    output_pos_file="S_orbicularis_gl_SNP_LD_chr${padded_chr}.pos"

    echo "Processing chromosome: chr${padded_chr}, output file: ${output_pos_file}"

    # Use awk to filter lines for the current chromosome
    awk -v chr="chr${padded_chr}" -v out="${output_pos_file}" '$1 == chr { print > out }' $pos_file

    # Count the number of lines in the created file
    line_count=$(wc -l < "$output_pos_file")
    echo "Number of lines in ${output_pos_file}: ${line_count}"

    # Extract corresponding lines from the beagle file, including the header
    output_beagle_file="S_orbicularis_gl_SNP_LD_chr${padded_chr}.beagle"
    head -n 1 "$beagle_file" > "$output_beagle_file"
    
    end_line=$((start_line + line_count - 1))
    sed -n "${start_line},${end_line}p" "$beagle_file" >> "$output_beagle_file"

    # Update the start line for the next iteration
    start_line=$((end_line + 1))

    echo "Created ${output_beagle_file} with the corresponding lines."
done
```

#### Check that all files with chromosome numbers have the 0 before the single digit chromosomes (i.e. chr01)

## 7.analyze

### 2.ngsLD_prune

### ngsLD/

- ngsLD v1.2.0, samtools v1.13, gcc, gsl
- Determine your number of sites for n_sites
- gzip -dc S_orbicularis_PWKO.beagle.gz \| wc -l –geno –pos –probs
  –n_ind the number of individuals you have –n_sites the number of sites
  found in your beagle or glf file –max_kb_dist maximum distance between
  SNPs (in Kb) to calculate LD; set to 0(zero) to disable filter
  –n_threads the number of threads to use –out the output file name

``` bash
#!/bin/bash

#SBATCH -J 2.ngsLD_prune.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 7-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p dept.les # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_2.ngsLD_prune_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_2.ngsLD_prune_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load modules
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/2.ngsLD_prune/2.ngsLD_prune.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/1.angsd_gl_SNP/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/2.ngsLD_prune/"
ngsLD="/home/bcarlson4/local/ngsTools/ngsLD/ngsLD"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

      # Determine linkage groups
      $ngsLD \
        --geno ${input_dir}S_orbicularis_gl_SNP_LD_chr01.beagle.gz \
        --pos ${input_dir}S_orbicularis_gl_SNP_LD_chr01.pos.gz \
        --probs \
        --n_ind 154 \
        --n_sites 6658 \
        --max_kb_dist 0 \
        --n_threads 32 \
        --out ${output_dir}S_orbicularis_chr01_prune.out
```

## 7.analyze

### 3.angsd_gl_SNP_LDp

### Generate genotype likelihoods and SNPs based on linkage disequilibrium

- angsd v0.94, samtools v1.13, gcc, gsl

``` bash
#!/bin/bash

#SBATCH -J 3.angsd_gl_SNP_LDp.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 0-00:05:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p test # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_3.angsd_gl_SNP_LDp_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_3.angsd_gl_SNP_LDp_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/3.angsd_gl_SNP_LDp/3.angsd_gl_SNP_LDp.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/2.ngsLD_prune/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/3.angsd_gl_SNP_LDp/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"
file_list="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/4.gatk_indel_realigner/gatk_indel_realigner_bam.list"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Input SNP list made from ngsLD
SNP_list=${input_dir}LDprune_snps.list

# Index SNP list
$angsd sites index $SNP_list

      # Call variant with informed linkage
      $angsd -b $file_list -ref $ref -out ${output_dir}S_orbicularis_gl_SNP_LDp \
        -doCounts 1 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
        -GL 2 -doGlf 2 -doMajorMinor 4 -doPost 2 -P 8 \
        -doIBS 1  -doCov 1 -makeMatrix 1 -sites $SNP_list
```

``` bash
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/3.angsd_gl_LDp/S_orbicularis_gl_LDp.covMat /Users/bailey/Documents/research/S_orbicularis/DNA/analyses/gl_SNP_LDp/
```

``` r
library(tidyverse)
#Define basedir
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/analyses/gl_SNP_LDp/"

#Load the covariance matrix
LDp_cov <- as.matrix(read.table(paste0(basedir, "S_orbicularis_gl_LDp.covMat"), header = F))

#We will also add a column with population assingments
pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/DNA/data/S_orbicularis_samples.csv")

LDp_pca <- eigen(LDp_cov) #perform the pca using the eigen function. 

LDp_eigenvectors = LDp_pca$vectors #extract eigenvectors 
LDp_pca.vectors = as_tibble(cbind(pop, data.frame(LDp_eigenvectors))) #combine with our population assignments

LDp_pca = ggplot(data = LDp_pca.vectors, aes(x=X1, y=X2, colour = Site)) + geom_point()

ggsave(filename = paste0(basedir, "LDp_pca_allSNPs_plot.pdf"), plot = pca) #change file path if data on your own computer

LDp_pca.eigenval.sum = sum(LDp_pca$values) #sum of eigenvalues
varPC1 <- (LDp_pca$values[1]/LDp_pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (LDp_pca$values[2]/LDp_pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (LDp_pca$values[3]/LDp_pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (LDp_pca$values[4]/LDp_pca.eigenval.sum)*100 #Variance explained by PC4
```

## 7.analyze

### 4.angsd_admix

### Determine population structure estimating number of populations

### Takes 20 hours

- angsd v0.94, samtools v1.13, gcc, gsl

``` bash
#!/bin/bash

#SBATCH -J 4.angsd_admix.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 1-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_4.angsd_admix_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_4.angsd_admix_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/4.angsd_admix/4.angsd_admix.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/3.angsd_gl_SNP_LDp/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/4.angsd_admix/"
ngs_admix="/home/bcarlson4/local/ngsTools/angsd/misc/NGSadmix"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

      # Admixture S_orbicularis_gl_LDp_vcf.beagle.gz
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 2 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K2_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 3 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K3_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 4 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K4_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 5 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K5_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 6 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K6_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 7 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K7_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 8 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K8_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 9 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K9_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 10 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K10_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 11 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K11_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 12 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K12_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 13 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K13_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 14 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K14_out

      # Admixture
        $ngs_admix -likes ${input_dir}S_orbicularis_gl_SNP_LDp.beagle.gz -K 15 \
        -P 8 -o ${output_dir}S_orbicularis_gl_SNP_LDp_admix_K15_out
```

``` bash
scp bcarlson4@login.rc.ucmerced.edu:/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/4.angsd_admix/S_orbicularis_gl_SNP_LDp_admix_K* /Users/bailey/Documents/research/S_orbicularis/DNA/analyses/admixture/
```

``` r
basedir <- "/Users/bailey/Documents/research/S_orbicularis/DNA/analyses/admixture/" 
library(tidyverse) #load the tidyverse package for formatting and plotting
library(viridis)

#Load the covariance matrix
admix = read_table(paste0(basedir, "S_orbicularis_gl_SNP_LDp_admix_K3_out.qopt"), col_names = F)

pop <- read.csv("/Users/bailey/Documents/research/S_orbicularis/DNA/data/S_orbicularis_samples.csv")

admix.id = as.data.frame(cbind(pop, admix[,-4]))
names(admix.id) = c("Sample", "Site", "Type", "Collection", "Details", "Island", "q1","q2","q3")
#,"q4","q5","q6","q7","q8","q9","q10","q11","q12","q13","q14","q15"

# Sort dataframe by multiple columns
sorted_admix.id <- admix.id[order(admix.id$Type, admix.id$Site), ]

# Generate a color palette with three colors
palette <- viridis_pal(alpha = 1, begin = 0.15, end = 0.85, option = "A")(15)
palette

row.names(sorted_admix.id) <- sorted_admix.id$Type

# Get the row names
row_names <- row.names(subset(sorted_admix.id, select = q1:q3))

# Create the bar plot
plot = barplot(t(as.matrix(subset(sorted_admix.id, select = q1:q3))),
               col=c("#251256FF", "#A1307EFF", "#FB8861FF", "#5C167FFF", "#E44F64FF", "#38106CFF","#F7735CFF", "#6D1D81FF", "#D6456CFF", "#4A1079FF", "#F05F5EFF", "#7E2482FF", "#C53C75FF", "#902A81FF", "#B3367AFF"),
               border=NA,
               names.arg=row_names,
               las=2, # Rotate labels vertically
               cex.axis = 0.5,
               cex.names = 0.5)  
```

## 7.analyze

### 5.angsd_saf

### Sample allele frequency

### Takes 2 days and 5 hours

- angsd v0.94, samtools v1.13, gcc, gsl -doSaf 1: perform multisample GL
  estimation Site

``` bash
#!/bin/bash

#SBATCH -J 5.angsd_saf_site.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_saf_site_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_saf_site_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/5.angsd_saf_site.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/5.angsd_saf/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Make the list file
rm -r ${output_dir}*_angsd_saf_site_bam.list
for POP in BAL CCM CCN CLM CRF FLK LMU NL1 NL2 NLK OCK OCM OTM TKC TLN ULN
do
echo $POP
touch ${output_dir}$POP'angsd_saf_site_bam.list'
    # Iterate over the files in the input directory
  for sample in ${output_dir}$POP/*_indel_realigner.bam; do
      # Append the full path to the file list
      echo $sample >> ${output_dir}$POP'angsd_saf_site_bam.list'
  done
done

for POP in BAL CCM CCN CLM CRF FLK LMU NL1 NL2 NLK OCK OCM OTM TKC TLN ULN
do
        echo $POP
        $angsd -b ${output_dir}$POP'angsd_saf_site_bam.list' -ref $ref -anc $ref -out ${output_dir}$POP'_S_orbicularis' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
            -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
            -GL 2 -doSaf 1 -P 8
done
```

Type

``` bash
#!/bin/bash

#SBATCH -J 5.angsd_saf_type.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_5.angsd_saf_type_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_5.angsd_saf_type_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/5.angsd_saf_type.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edit_alignments/4.gatk_indel_realigner/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyze/5.angsd_saf/"
ref="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Make the list file
rm -r ${output_dir}*angsd_saf_type_bam.list
for POP in Ocean Holomictic Meromictic
do
echo $POP
touch ${output_dir}$POP'angsd_saf_type_bam.list'
    # Iterate over the files in the input directory
  for sample in ${output_dir}$POP/*_indel_realigner.bam; do
      # Append the full path to the file list
      echo $sample >> ${output_dir}$POP'_angsd_saf_type_bam.list'
  done
done

for POP in Ocean Holomictic Meromictic
do
        echo $POP
        $angsd -b ${output_dir}$POP'_type_bam.list' -ref $ref -anc $ref -out ${output_dir}$POP'_S_orbicularis' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
            -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
            -GL 2 -doSaf 1 -P 8
done
```

## 7.analyze

### 6.angsd_sfs

### Site frequency spectrum

### Takes days

- angsd v0.94, samtools v1.13, gcc, gsl

``` bash
#!/bin/bash

#SBATCH -J 6.angsd_sfs.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
#SBATCH -t 3-00:00:00 # --time= Runtime in minutes. Default is 10 minutes.
#SBATCH -p long # --partition Partition to submit to the standard compute node.
#SBATCH -o slurm_6.angsd_sfs_%j.out # --output= Standard out goes to this file.
#SBATCH -e slurm_6.angsd_sfs_%j.err # Standard err goes to this file.
#SBATCH --mail-user bcarlson4@ucmerced.edu # this is the email for notifications.
#SBATCH --mail-type ALL # this specifies what events will be emailed.

# Load module
module load samtools/1.13

module load gcc/12.2.0

module load gsl/2.7

# Execute your processing script
/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/6.angsd_sfs.sh
```

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/4.realigned/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/5.angsd_saf/"
anc="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference_indexed/ref.fna"
realSFS="/home/bcarlson4/local/ngsTools/angsd/misc/realSFS"

for POP in BAL CCM CCN CLM CRF FLK LMU NL1 NL2 NLK OCK OCM OTM TKC TLN ULN
do
        echo $POP
        $realSFS ${output_dir}$POP.saf.idx -fold 1 > ${output_dir}$POP.sfs
done

POP2=CCM
POP3=CCN
POP4=CRF
POP5=LMU
POP6=NL1
POP7=NL2
POP8=OCK
POP9=OCM
POP10=TKC
for POP in BAL
do
        echo $POP
        $realSFS ${output_dir}$POP.saf.idx ${output_dir}$POP2.saf.idx > $BASEDIR/results/ocean.sfs
done

POP2=ULN
for POP in FLK
do
        echo $POP
        $realSFS ${output_dir}$POP.saf.idx ${output_dir}$POP2.saf.idx > ${output_dir}holomictic.sfs
done

POP2=NLK
POP3=OTM
POP4=TLN
for POP in CLM
do
        echo $POP
        $realSFS ${output_dir}$POP.saf.idx ${output_dir}$POP2.saf.idx ${output_dir}$POP3.saf.idx ${output_dir}$POP4.saf.idx > ${output_dir}meromictic.sfs
done
```

## 7.analyze

### 7.angsd_fst

### 

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/4.realigned/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/5.angsd_sfs/"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"
realSFS="/home/bcarlson4/local/ngsTools/angsd/misc/realSFS"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

      # Fst
      realSFS fst index Results/MAQU.saf.idx Results/MBNS.saf.idx Results/JIGA.saf.idx -sfs Results/MAQU.MBNS.sfs -sfs Results/MAQU.JIGA.sfs -sfs Results/MBNS.JIGA.sfs -fstout Results/JIGA.pbs -whichFst 1
```

## 7.analyze

### 8.angsd_nd

### 

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/4.realigned/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/7.analyzed/5.angsd_sfs/"
anc="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/4.reference_indexed/ref.fna"
angsd="/home/bcarlson4/local/ngsTools/angsd/angsd"
thetaStat="/home/bcarlson4/local/ngsTools/angsd/misc/thetaStat"

for POP in BAL CCM CCN CLM CRF FLK LMU NL1 NL2 NLK OCK OCM OTM TKC TLN ULN
do
        $angsd -b ${input_dir}$POP'_bam.list' -anc $anc -out ${output_dir}$POP'_S_orbicularis' \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
            -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 5 -setMaxDepth 60 -doCounts 1 \
            -GL 2 -doSaf 1 -P 8 -doThetas 1 -pest ${output_dir}$POP'_S_orbicularis.sfs'
done


for POP in BAL CCM CCN CLM CRF FLK LMU NL1 NL2 NLK OCK OCM OTM TKC TLN ULN
do
      # estimate for the whole region
      $thetaStat do_stat ${output_dir}$POP'_S_orbicularis.thetas.idx'
      # perform a sliding-window analysis
      $thetaStat do_stat ${output_dir}$POP'_S_orbicularis.thetas.idx' \
      -win 10000 -step 1000 -outnames ${output_dir}$POP'_S_orbicularis.thetas.windows'
done
```

## 7.analyze

### 9.ngsRelate_ib

### 

``` bash
git clone https://github.com/ANGSD/NgsRelate.git
```

#### Check file qualities

## picard validate summary

### picard_validate_summary.sh

### Takes ? hours

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

## picard validate verbose

### picard_validate_verbose.sh

### Takes ? hours

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

# Extra/unused code

## create_mtDNA

### Create the mitochondrial genome

### Takes

-j \[seqid\]  
-a \[assembly.fasta\]  
-r \[genbank_reference.gb\]  
-o \[genetic_code\] -p \[threads\] -m \[memory\]
<https://www.ncbi.nlm.nih.gov/nuccore/?term=%22mitochondrion%22+and+%22complete+genome%22+and+%22Apogoninae%22>

``` bash
#!/bin/bash

#SBATCH -J mitofinder_mtDNA.batch # --job-name= Name for your job.
#SBATCH -N 1 # --nodes= Number of nodes to spread cores across - default is 1.
#SBATCH --mem-per-cpu 240G # the amount of memory per core to request in MB.
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

## samtools fixmate

### Fix the mate pairing information in your sequences, file name samtools_fixmate.sh and corresponds to samtools_fixmate.slurm

### Takes 2 hours

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

## samtools sort

### samtools_sort.sh

### Takes 2.5 hours

- samtools manual <https://www.htslib.org/doc/samtools.html>

``` bash
#!/bin/bash

# Assuming your files are in the current directory
input_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/1.viewed/L00*/"
output_dir="/home/bcarlson4/borgstore/bcarlson4/Sphaeramia/DNA/6.edited_alignments/2.sorted/"

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

## samtools merge

### Merge lane files together, file name samtools_merge.sh and corresponds to samtools_merge.slurm

### Takes 21 hours

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

## gatk br bqsr

### gatk_br_bqsr.sh

### 

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
