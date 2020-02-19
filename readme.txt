Assignemnt 3
Author: Christopher Nguyen

How to run:
1. Must have following files in the same directory:
    pipeline.py
    dgorgon_reference.fa
    hawkins_pooled_sequences.fastq
    harrington_clinical_data.txt
2. Have samtools, venv, and BWA installed
3. Change the absolute pathway for BWA and Samtools within pipeline.py. By default the pathway are set as if /BWA and /samstool directories are within the directory containing the files listed in step 1
4. When steps above are met, run pipeline.py

Pipeline.py takes the subjects who developed mold in their ear listed in the harrington_clinical_data.txt and create individual .fastq for each patient within folder fastqs. The sequences are added into their individual .fastq files based on the barcode indicator within the clinical_data.txt.

After the .fastq files are created for each subject, the .fastq files are then converted to .bam and indexed through using BWA and the dgorgon_reference.fa. The .bam files are stored in the bams folder.

Dgorgon_reference.fa contains the wildtype sequence. Once .bam files are sorted and indexed, pysam is utilized to compare the sequences of each subject to the wildtype sequence.

The sequences and bases that do not match the wild type sequenced is reported. The findings of each subject is reported in the text file report.txt.
