from Bio import SeqIO
Bed_Path=''
Bed1_path=''
Bed2_path=''
def remove_overlap():
    # translocation_list = []
    # with open('alltype.bed', 'r') as f:
    #     for line in f:
    #         seq = line.strip().split('\t')
    #         if seq[3] == 'reciprocal translocation':
    #             start = int(seq[4].split(':')[2])
    #             length = int(seq[2]) - int(seq[1])
    #             translocation_list.append((start - length, start + length))
    # # print(translocation_list)
    translocation_list = []
    with open('alltype.bed', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[3] == 'translocation cut-paste':
                start = int(seq[4].split(':')[2])
                length = int(seq[2]) - int(seq[1])
                translocation_list.append((start - length, start + length))
    # print(translocation_list)
    chrom = ''
    last_start = 0
    last_end = 0
    file = open('visor/del_hp2_rmo.bed', 'w')
    with open('visor/del.21_22.hp2.bed', 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[0] != chrom:
                last_start = int(seq[1])
                last_end = int(seq[2])
                chrom = seq[0]
                continue
            overlap = False
            if int(seq[1]) <= last_end:
                overlap = True
                print(line)
            # if seq[0] == 'chr2':
            #     for item in translocation_list:
            #         if (int(seq[1]) <= item[0] and int(seq[2]) >= item[0]) or (item[0] <= int(seq[1]) <= item[1]):
            #             overlap = True
            if not overlap:
                file.write(line)
                last_start = int(seq[1])
                last_end = int(seq[2])
    file.close()

def combine_bed():
    output = open(Bed_Path, 'w')
    h1_dict = dict()
    h2_dict = dict()
    with open(Bed1_path, 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[0] not in h1_dict:
                h1_dict[seq[0]] = list()
            h1_dict[seq[0]].append(seq)
    with open(Bed2_Path, 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            if seq[0] not in h2_dict:
                h2_dict[seq[0]] = list()
            h2_dict[seq[0]].append(seq)
    for chrom in h1_dict:
        h1_idx = 0
        h2_idx = 0
        while h1_idx < len(h1_dict[chrom]) and h2_idx < len(h2_dict[chrom]):
            if int(h1_dict[chrom][h1_idx][1]) < int(h2_dict[chrom][h2_idx][1]):
                output.write('\t'.join(h1_dict[chrom][h1_idx]))
                output.write('\t1|0\n')
                h1_idx += 1
            elif int(h2_dict[chrom][h2_idx][1]) < int(h1_dict[chrom][h1_idx][1]):
                output.write('\t'.join(h2_dict[chrom][h2_idx]))
                output.write('\t0|1\n')
                h2_idx += 1
            else:
                if int(h2_dict[chrom][h2_idx][2]) == int(h1_dict[chrom][h1_idx][2]) and len(h1_dict[chrom][h1_idx][4]) == len(h2_dict[chrom][h2_idx][4]):
                    # 1|1
                    output.write('\t'.join(h1_dict[chrom][h1_idx][:-1]))
                    output.write('\t1|1\n')
                else:
                    output.write('\t'.join(h1_dict[chrom][h1_idx]))
                    output.write('\t1|0\n')
                    output.write('\t'.join(h2_dict[chrom][h2_idx]))
                    output.write('\t0|1\n')
                h1_idx += 1
                h2_idx += 1
        while h1_idx < len(h1_dict[chrom]):
            output.write('\t'.join(h1_dict[chrom][h1_idx]))
            output.write('\t1|0\n')
            h1_idx += 1
        while h2_idx < len(h2_dict[chrom]):
            output.write('\t'.join(h2_dict[chrom][h2_idx]))
            output.write('\t0|1\n')
            h2_idx += 1
    output.close()

def Generation_VCF_header(file):
    # General header
    file.write("##fileformat=VCFv4.2\n")
    import time
    file.write("##fileDate=%s\n"%(time.strftime('%Y-%m-%d %H:%M:%S %w-%Z',time.localtime())))
    #file.write("##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr2,length=242193529>\n##contig=<ID=chr3,length=198295559>\n##contig=<ID=chr4,length=190214555>\n##contig=<ID=chr5,length=181538259>\n##contig=<ID=chr6,length=170805979>\n##contig=<ID=chr7,length=159345973>\n##contig=<ID=chr8,length=145138636>\n##contig=<ID=chr9,length=138394717>\n##contig=<ID=chr10,length=133797422>\n##contig=<ID=chr11,length=135086622>\n##contig=<ID=chr12,length=133275309>\n##contig=<ID=chr13,length=114364328>\n##contig=<ID=chr14,length=107043718>\n##contig=<ID=chr15,length=101991189>\n##contig=<ID=chr16,length=90338345>\n##contig=<ID=chr17,length=83257441>\n##contig=<ID=chr18,length=80373285>\n##contig=<ID=chr19,length=58617616>\n##contig=<ID=chr20,length=64444167>\n##contig=<ID=chr21,length=46709983>\n##contig=<ID=chr22,length=50818468>\n##contig=<ID=chrX,length=156040895>\n##contig=<ID=chrY,length=57227415>")
    # Specific header
    file.write("##contig=<ID=chr1,length=225280621>\n##contig=<ID=chr2,length=238204518>\n##contig=<ID=chr3,length=194797135>\n##contig=<ID=chr4,length=187661676>\n##contig=<ID=chr5,length=177695260>\n##contig=<ID=chr6,length=167395066>\n##contig=<ID=chr7,length=155353663>\n##contig=<ID=chr8,length=142888922>\n##contig=<ID=chr9,length=120143431>\n##contig=<ID=chr10,length=131314738>\n##contig=<ID=chr11,length=131129516>\n##contig=<ID=chr12,length=130481393>\n##contig=<ID=chr13,length=95589878>\n##contig=<ID=chr14,length=88289540>\n##contig=<ID=chr15,length=81694766>\n##contig=<ID=chr16,length=78884753>\n##contig=<ID=chr17,length=77498584>\n##contig=<ID=chr18,length=74657229>\n##contig=<ID=chr19,length=55808983>\n##contig=<ID=chr20,length=59505520>\n##contig=<ID=chr21,length=35106642>\n##contig=<ID=chr22,length=34894545>\n##contig=<ID=chrX,length=151100560>\n##contig=<ID=chrY,length=22984529>\n")
    # ALT
    file.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence relative to the reference\">\n")
    file.write("##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">\n")
    file.write("##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference\">\n")
    file.write("##ALT=<ID=INV,Description=\"Inversion of reference sequence\">\n")
    file.write("##ALT=<ID=BND,Description=\"Breakend of translocation\">\n")

    # INFO
    file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n")
    file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    # file.write("##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    file.write("##FILTER=<ID=q5,Description=\"Quality below 5\">\n")
    # FORMAT
    # file.write("\n")
    file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    # file.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n")
    # file.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n")
    # file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n")
    # file.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n")

def generate_output():
    
    '''
    Generation of VCF format file.
    VCF version: 4.2
    '''

    # genotype_trigger = TriggerGT[args.genotype]
    sv_trans = {'insertion':'INS', 'deletion':'DEL', 'inversion':'INV', 'tandem duplication':'DUP', 'reciprocal translocation':'BND'}
    svid = dict()
    svid["INS"] = 0
    svid["DEL"] = 0
    svid["BND"] = 0
    svid["DUP"] = 0
    svid["INV"] = 0

    file = open('reslut.vcf', 'w')
    ref_g = SeqIO.to_dict(SeqIO.parse('fullRef.fasta', "fasta"))
    Generation_VCF_header(file)
    file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
    with open(Bed_path, 'r') as f:
        for line in f:
            seq = line.strip().split('\t')
            chrom = seq[0]
            start = int(seq[1])
            end = int(seq[2])
            svtype = sv_trans[seq[3]]
            gt = seq[-1]
            if svtype == 'INS':
                file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}\n".format(
                    CHR = chrom, 
                    POS = str(start), 
                    ID = "%s.%d"%('INS', svid['INS']),
                    REF = str(ref_g[chrom].seq[max(start-1, 0)]),
                    ALT = str(ref_g[chrom].seq[max(start-1, 0)]) + seq[4], 
                    INFO = "SVTYPE=%s;SVLEN=%d;END=%d"%(svtype, len(seq[4]), start), 
                    FORMAT = "GT", 
                    GT = gt,
                    QUAL = '.',
                    PASS = 'PASS'))
            if svtype == 'DEL':
                file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}\n".format(
                    CHR = seq[0], 
                    POS = str(start), 
                    ID = "%s.%d"%('DEL', svid['DEL']),
                    REF = str(ref_g[chrom].seq[max(start-1, 0):end]),
                    ALT = str(ref_g[chrom].seq[max(start-1, 0)]), 
                    INFO = "SVTYPE=%s;SVLEN=%d;END=%d"%(svtype, start-end, end), 
                    FORMAT = "GT", 
                    GT = gt,
                    QUAL = '.',
                    PASS = 'PASS'))
            if svtype == 'INV' or svtype == 'DUP':
                file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}\n".format(
                    CHR = seq[0], 
                    POS = str(start), 
                    ID = "%s.%d"%('INV', svid['INV']) if svtype == 'INV' else "%s.%d"%('DUP', svid['DUP']),
                    REF = str(ref_g[chrom].seq[start]),
                    ALT = "<%s>"%(svtype), 
                    INFO = "SVTYPE=%s;SVLEN=%d;END=%d"%(svtype, end-start, end), 
                    FORMAT = "GT", 
                    GT = gt,
                    QUAL = '.',
                    PASS = 'PASS'))
            if svtype == 'BND':
                chr2_info = seq[4].split(':')
                chr2 = chr2_info[1]
                chr2_start = int(chr2_info[2])
                file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}\n".format(
                    CHR = seq[0], 
                    POS = str(start), 
                    ID = "%s.%d"%('BND', svid['BND']),
                    REF = str(ref_g[chrom].seq[start]),
                    ALT = "N[%s:%d["%(chr2, chr2_start), 
                    INFO = "SVTYPE=%s;SVLEN=%d;CHR2=%s;END=%d"%(svtype, end-start, chr2, chr2_start), 
                    FORMAT = "GT", 
                    GT = gt,
                    QUAL = '.',
                    PASS = 'PASS'))

            svid[svtype] += 1
    file.close()

if __name__ == '__main__':
    # remove_overlap()
    combine_bed()
    generate_output()
