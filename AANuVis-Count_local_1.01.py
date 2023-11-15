#this is a local version of AANuVis, dependency: seqkit, prodigal, and of course python
#the code imports all sequences from the fasta file, trains and translates with prodigal
#_sorted.csv are usage counts of each genome; _perced.csv are percentaged count of each genome, suited for concecutive visualization unless otherwise needed
#AA count order: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
#DiNucl count order: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT
#Tetra count order: AAAA,AAAC,AAAG,AAAT,AACA,AACC,AACG,AACT,AAGA,AAGC,AAGG,AAGT,AATA,AATC,AATG,AATT,ACAA,ACAC,ACAG,ACAT,ACCA,ACCC,ACCG,ACCT,ACGA,ACGC,ACGG,ACGT,ACTA,ACTC,ACTG,ACTT,AGAA,AGAC,AGAG,AGAT,AGCA,AGCC,AGCG,AGCT,AGGA,AGGC,AGGG,AGGT,AGTA,AGTC,AGTG,AGTT,ATAA,ATAC,ATAG,ATAT,ATCA,ATCC,ATCG,ATCT,ATGA,ATGC,ATGG,ATGT,ATTA,ATTC,ATTG,ATTT,CAAA,CAAC,CAAG,CAAT,CACA,CACC,CACG,CACT,CAGA,CAGC,CAGG,CAGT,CATA,CATC,CATG,CATT,CCAA,CCAC,CCAG,CCAT,CCCA,CCCC,CCCG,CCCT,CCGA,CCGC,CCGG,CCGT,CCTA,CCTC,CCTG,CCTT,CGAA,CGAC,CGAG,CGAT,CGCA,CGCC,CGCG,CGCT,CGGA,CGGC,CGGG,CGGT,CGTA,CGTC,CGTG,CGTT,CTAA,CTAC,CTAG,CTAT,CTCA,CTCC,CTCG,CTCT,CTGA,CTGC,CTGG,CTGT,CTTA,CTTC,CTTG,CTTT,GAAA,GAAC,GAAG,GAAT,GACA,GACC,GACG,GACT,GAGA,GAGC,GAGG,GAGT,GATA,GATC,GATG,GATT,GCAA,GCAC,GCAG,GCAT,GCCA,GCCC,GCCG,GCCT,GCGA,GCGC,GCGG,GCGT,GCTA,GCTC,GCTG,GCTT,GGAA,GGAC,GGAG,GGAT,GGCA,GGCC,GGCG,GGCT,GGGA,GGGC,GGGG,GGGT,GGTA,GGTC,GGTG,GGTT,GTAA,GTAC,GTAG,GTAT,GTCA,GTCC,GTCG,GTCT,GTGA,GTGC,GTGG,GTGT,GTTA,GTTC,GTTG,GTTT,TAAA,TAAC,TAAG,TAAT,TACA,TACC,TACG,TACT,TAGA,TAGC,TAGG,TAGT,TATA,TATC,TATG,TATT,TCAA,TCAC,TCAG,TCAT,TCCA,TCCC,TCCG,TCCT,TCGA,TCGC,TCGG,TCGT,TCTA,TCTC,TCTG,TCTT,TGAA,TGAC,TGAG,TGAT,TGCA,TGCC,TGCG,TGCT,TGGA,TGGC,TGGG,TGGT,TGTA,TGTC,TGTG,TGTT,TTAA,TTAC,TTAG,TTAT,TTCA,TTCC,TTCG,TTCT,TTGA,TTGC,TTGG,TTGT,TTTA,TTTC,TTTG,TTTT
#Usage: AANuVis-count_local.py import.fasta , no none-verbose option is provided, so recommended nohup (-) & for huge files 
#1.01 version fixed the issue if there are special characters in seq names

import os
import re
import csv
import sys
def run_prodigal(seq_id):
    with open(f'{seq_id}.fasta', 'r') as file:
        sequence = file.read()
    with open(f'{seq_id}_modified.fasta', 'w') as file:
        file.write(f'>{seq_id}\n')
        file.write(sequence)
    os.system(f"prodigal -i '{seq_id}_modified.fasta' -a '{seq_id}.faa' -t '{fasta_file}_train'")
    os.system(f"rm '{seq_id}_modified.fasta'")
    
def delete_files(seq_id):
    os.system(f"rm '{seq_id}.fasta'")
    os.system(f"rm '{seq_id}.faa'")

def calculate_dinucleotide_count(seq_id):
    nucleotides = ['A', 'C', 'G', 'T']
    dinucleotides = [a+b for a in nucleotides for b in nucleotides]
    dinucleotide_counts = {dinucleotide: 0 for dinucleotide in dinucleotides}
    with open(f'{seq_id}.fasta', 'r') as file:
        sequence = file.read().replace('\n', '')
    for i in range(len(sequence)-1):
        dinucleotide = sequence[i:i+2]
        if dinucleotide in dinucleotide_counts:
            dinucleotide_counts[dinucleotide] += 1
    with open('di-nucl.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        for dinucleotide, count in dinucleotide_counts.items():
            writer.writerow([dinucleotide, count])
def calculate_tetranucleotide_count(seq_id):
    nucleotides = ['A', 'C', 'G', 'T']
    tetranucleotides = [a+b+c+d for a in nucleotides for b in nucleotides for c in nucleotides for d in nucleotides]
    tetranucleotide_counts = {tetranucleotide: 0 for tetranucleotide in tetranucleotides}
    with open(f'{seq_id}.fasta', 'r') as file:
        sequence = file.read().replace('\n', '')
    for i in range(len(sequence)-3):
        tetranucleotide = sequence[i:i+4]
        if tetranucleotide in tetranucleotide_counts:
            tetranucleotide_counts[tetranucleotide] += 1
    with open('tetra-nucl.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        for tetranucleotide, count in tetranucleotide_counts.items():
            writer.writerow([tetranucleotide, count])
def calculate_amino_acid_usage(seq_id):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_acid_counts = {aa: 0 for aa in amino_acids}
    with open(f'{seq_id}.faa', 'r') as file:
        lines = file.readlines()
    for line in lines:
        if line.startswith('>'):
            continue
        sequence = line.strip()
        for aa in amino_acids:
            amino_acid_counts[aa] += sequence.count(aa)
    with open('AA-usage.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        for aa, count in amino_acid_counts.items():
            writer.writerow([aa, count])
def process_fasta_file(fasta_file):
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
    seq_id = ''
    seq = ''
    for line in lines:
        if line.startswith('>'):
            if seq_id != '':
                with open(f'{seq_id}.fasta', 'w') as file:
                    file.write(seq)
                run_prodigal(seq_id)
                calculate_dinucleotide_count(seq_id)
                calculate_tetranucleotide_count(seq_id)
                calculate_amino_acid_usage(seq_id)
                delete_files(seq_id)
                seq = ''
            seq_id = re.search(r'>(\S+)', line).group(1)
        else:
            seq += line.strip()
    if seq_id != '':
        with open(f'{seq_id}.fasta', 'w') as file:
            file.write(seq)
        run_prodigal(seq_id)
        calculate_dinucleotide_count(seq_id)
        calculate_tetranucleotide_count(seq_id)
        calculate_amino_acid_usage(seq_id)
        delete_files(seq_id)
fasta_file = sys.argv[1]
with open('di-nucl.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Dinucleotide', 'Count'])
with open('tetra-nucl.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Tetranucleotide', 'Count'])
with open('AA-usage.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Amino Acid', 'Count'])
prodigal_train = 'prodigal -i {} -t {}_train'.format(fasta_file,fasta_file)
os.system(prodigal_train)
process_fasta_file(fasta_file)

seqkit_command = 'seqkit seq -n {} -o {}_string.txt'.format(fasta_file, fasta_file)
os.system(seqkit_command)

def collect_nc(line):
    return "".join(line.strip())
def split_data(data, num_rows):
    new_data = []
    for i in range(0, len(data), num_rows):
        new_data.append(data[i:i+num_rows])
    return new_data
with open('AA-usage.csv', 'r') as file:
    csv_data = list(csv.reader(file))
    data = [row[1] for row in csv_data[1:]]
    new_data = split_data(data, 20)
with open('{}_string.txt'.format(fasta_file), 'r') as file:
    lines = file.readlines()
    nc_lines = [collect_nc(line) for line in lines]
for i, line in enumerate(new_data):
    line.insert(0, nc_lines[i])
with open('AA-usage_sorted.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(new_data)
with open('di-nucl.csv', 'r') as file:
    csv_data = list(csv.reader(file))
    data = [row[1] for row in csv_data[1:]]
    new_data = split_data(data, 16)
for i, line in enumerate(new_data):
    line.insert(0, nc_lines[i])
with open('di-nucl_sorted.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(new_data)
with open('tetra-nucl.csv', 'r') as file:
    csv_data = list(csv.reader(file))
    data = [row[1] for row in csv_data[1:]]
    new_data = split_data(data, 256)
for i, line in enumerate(new_data):
    line.insert(0, nc_lines[i])
with open('tetra-nucl_sorted.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(new_data)

def calculate_percentage(row):
    total = sum(int(value) for value in row[1:])
    if total == 0:
        total = 1
    percentages = [round(float(value) / total, 9) for value in row[1:]]
    return [row[0]] + percentages
def process_csv(input_file, output_file):
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        rows = list(reader)
    rows = [calculate_percentage(row) for row in rows]
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(rows)
process_csv('AA-usage_sorted.csv', 'AA-usage_perced.csv')
process_csv('di-nucl_sorted.csv', 'di-nucl_perced.csv')
process_csv('tetra-nucl_sorted.csv', 'tetra-nucl_perced.csv')
