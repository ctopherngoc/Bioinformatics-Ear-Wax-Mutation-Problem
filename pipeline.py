import os
import pysam
text_file = open("report.txt", "w")
Guide = "AGCGGTCATAAGTGGTACATTACGAGATTCGGAGTACCATAGATTCGCATGAATCCCTGTGGATACGAGAGTGTGAGATATATGTACGCCAATCCAGTGTGATACCCATGAGATTTAGGACCGATGATGGTTGAGGACCAAGGATTGACCCGATGGATGCAGATTTGACCCCAGATAGAATAAATGCGATGAGATGATTTGGCCGATAGATAGATAG"

#Create directories and create absolute paths. REMEMBER TO CHANGE BWA AND SAMTOOLS PATHWAYS
os.mkdir("fastqs")
os.mkdir("sams")
os.mkdir("bams")
script_dir = os.path.dirname(__file__)
rel_path_sams = "sams/"
rel_path_fastqs = "fastqs/"
rel_path_bams = "bams/"
abs_file_path_sams = os.path.join(script_dir, rel_path_sams)
abs_file_path_fastqs = os.path.join(script_dir, rel_path_fastqs)
abs_file_path_BWA = "/Users/chris/Desktop/week6/BWA/"  ######### change to absolute path for BWA for local machine
abs_file_path_samtools = "/Users/chris/Desktop/week6/samtools/" ##########change to absolute path for samtools local machine
abs_file_path_bams = os.path.join(script_dir, rel_path_bams)

# Create list of patient names, barcode, and color
f = open("harrington_clinical_data.txt", "r")
lines = f.readlines()
name = []
code = []
color = []

# Transfer hawken_pooled_sequence.fastq into list
for x in lines:
    name.append(x.split('\t')[0])
    color.append(x.split('\t')[1])
    code.append(x.split('\t')[2])
f.close()

# Remove first entry in each list
code = [x.strip('\n') for x in code]
name.remove('Name')
code.remove('Barcode')
color.remove('Color')

# create counters
count = 0
blueprint = []  # list containing hawkin_pooled.fastq data

#append hackins_pooled_Sequence into blueprint
g = open("hawkins_pooled_sequences.fastq", 'r')
for x in g:
    blueprint.append(x.split('\n')[0])
g.close()

# Create list for sequences that match patient barcode
for (i, j) in zip(name, code):
    i_list = []
    count = 1
    while count <= len(blueprint):
        if blueprint[count].startswith(j):
            i_list.append(blueprint[count-1])
            i_list.append(blueprint[count])
            i_list.append(blueprint[count+1])
            i_list.append(blueprint[count+2])
            count += 4
        else:
            count += 4

        # count = 1
        while count < (len(i_list)):
            if i_list[count].startswith(j):
                i_list[count] = (i_list[count].lstrip(j))
                count += 4
            else:
                count += 4

    count = 1
    while count < (len(i_list)):
        code_len = (len(j))
        i_list[count] = i_list[count][code_len:]
        i_list[count+2] = i_list[count+2][code_len:]
        temp = i_list[count+2].split('DD', 1)[0]
        temp = (temp.split('FF', 1)[0])
        temp = (temp.split('FD', 1)[0])
        temp = (temp.split('DF', 1)[0])
        temp = temp.rstrip('D')
        i_list[count+2] = temp.rstrip('F')
        temp_var = len(i_list[count+2])
        i_list[count] = i_list[count][:temp_var]
        count += 4
    file = open(abs_file_path_fastqs + i + "_trimmed.fastq", "w")
    for line in i_list:
        file.write(line + "\n")
    file.close()

#run BWA index and mem through loop
os.system(abs_file_path_BWA + "bwa index " + script_dir + "/dgorgon_reference.fa")

#create bam and bam.bai files
for i in name:
    os.system(abs_file_path_BWA + "bwa mem dgorgon_reference.fa " + abs_file_path_fastqs + i + "_trimmed.fastq > " + abs_file_path_sams + i +".sam")
    os.system(abs_file_path_samtools + "samtools view -bS " + abs_file_path_sams + i + ".sam > " + abs_file_path_bams + i + ".bam")
    os.system(abs_file_path_samtools + "samtools sort " + abs_file_path_bams + i + ".bam -o " + abs_file_path_bams + i + ".sorted.bam")
    os.system(abs_file_path_samtools + "samtools index " + abs_file_path_bams + i + ".sorted.bam")
    os.system("rm -r " + abs_file_path_bams + i + (".bam"))

#remove not needed files and directories
os.system("rm -r " + abs_file_path_sams)
os.system("rm -r " + script_dir + "/dgorgon_reference.fa.amb")
os.system("rm -r " + script_dir + "/dgorgon_reference.fa.ann")
os.system("rm -r " + script_dir + "/dgorgon_reference.fa.bwt")
os.system("rm -r " + script_dir + "/dgorgon_reference.fa.pac")
os.system("rm -r " + script_dir + "/dgorgon_reference.fa.sa")


def pileup() -> object:
    count = 0
    text_file.write(
        'The green mold was caused by a mutation in position 134. The wildtype base was G and the mutation was C.\n')
    text_file.write(
        'The black mold was caused by a mutation in position 13. The wildtype base was G and the mutation was C.\n')
    text_file.write(
        'The orange mold was caused by a mutation in position 83. The wildtype base was G and the mutation was C.\n')
    text_file.write(
        'The yellow mold was caused by a mutation in position 120. The wildtype base was C and the mutation was G.\n\n')
    for i in name:
        # test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
        samfile = pysam.AlignmentFile(abs_file_path_bams + i + ".sorted.bam", "rb")

        for pileupcolumn in samfile.pileup(): ###create base counts for each subject
            base_A = 0
            base_T = 0
            base_C = 0
            base_G = 0

            # use a dictionary to count up the bases at each position
            ntdict = {}
            ntdict.update({'name': i})
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:

                    base = pileupread.alignment.query_sequence[pileupread.query_position]

                    if base == "A":
                        base_A += 1
                    elif base == "T":
                        base_T += 1
                    elif base == "C":
                        base_C += 1
                    else:
                        base_G += 1
            ###if any of base is not equal to 0 or sequence length
            if (base_A != pileupcolumn.n and base_A != 0) or (base_T != pileupcolumn.n and base_T != 0) or (
                    base_C != pileupcolumn.n and base_C != 0) or (base_G != pileupcolumn.n and base_G != 0):
                ntdict.update({"Position": pileupcolumn.pos})
                ntdict.update({"Count": pileupcolumn.n})

                ##compares refence strand to base count and appends keys and values to dictionary
                if str(Guide[pileupcolumn.pos]) == "A":
                    if base_T > 0:
                        ntdict.update({"M count": (base_T)})
                        ntdict.update({"Freq": ((base_T / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'T'})
                    if base_C > 0:
                        ntdict.update({"M count": (base_C)})
                        ntdict.update({"Freq": ((base_C / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'C'})
                    else:  # base G
                        ntdict.update({"M count": (base_G)})
                        ntdict.update({"Freq": ((base_G / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'G'})
                elif str(Guide[pileupcolumn.pos]) == "T":
                    if base_A > 0:
                        ntdict.update({"M count": (base_A)})
                        ntdict.update({"Freq": ((base_A / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'A'})
                    if base_C > 0:
                        tdict.update({"M count": (base_C)})
                        ntdict.update({"Freq": ((base_C / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'C'})
                    else:  # base G
                        ntdict.update({"M count": (base_G)})
                        ntdict.update({"Freq": ((base_G / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'G'})
                elif str(Guide[pileupcolumn.pos]) == "C":
                    if base_T > 0:
                        ntdict.update({"M count": (base_T)})
                        ntdict.update({"Freq": ((base_T / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'T'})
                    if base_A > 0:
                        tdict.update({"M count": (base_A)})
                        ntdict.update({"Freq": ((base_A / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'A'})
                    else:  # base G
                        ntdict.update({"M count": (base_G)})
                        ntdict.update({"Freq": ((base_G / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'G'})
                else:  # G
                    if base_T > 0:
                        ntdict.update({"M count": (base_T)})
                        ntdict.update({"Freq": ((base_T / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'T'})
                    if base_C > 0:
                        ntdict.update({"M count": (base_C)})
                        ntdict.update({"Freq": ((base_C / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'C'})
                    else:  # base A
                        ntdict.update({"M count": (base_A)})
                        ntdict.update({"Freq": ((base_A / pileupcolumn.n) * 100)})
                        ntdict.update({'Mutation': 'A'})

                ##appends report into text file based on color.
                if color[count] == "Orange":
                    haha = ('Subject ' + str(ntdict["name"]) + ' has an ' + color[count] + ' mold, ' + str(ntdict["Count"]) + ' reads, and had ' + "{:.0f}".format(ntdict['Freq']) + "% of the reads at position " + str(ntdict['Position']) + ' had the mutation ' + str(ntdict['Mutation']) + '.')
                    text_file.write(haha + '\n')
                    count += 1
                else:
                    haha = ('Subject ' + str(ntdict["name"]) + ' has a ' + color[count] + ' mold, ' + str(ntdict["Count"]) + ' reads, and had ' + "{:.0f}".format(ntdict['Freq']) + "% of the reads at position " + str(ntdict['Position']) + ' had the mutation ' + str(ntdict['Mutation']) + '.')
                    text_file.write(haha + '\n')
                    count += 1

    samfile.close()
    text_file.close()


if __name__ == "__main__":
    pileup()
