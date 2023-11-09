dict = {}
gene_dict = {}
import argparse

print("start")


def calculate_heterozygosity(genotype_list):
    value = 0
    one_one_value = 0
    count = 0
    for element in genotype_list:
        genotype = element.split(":")[0]
        if genotype[0] == '.' or genotype[2] == '.':
            continue
        if genotype[0] == genotype[2] and int(genotype[0]) + int(genotype[2]) != 0:
            value += 1
        if genotype[0] == '1' and genotype[2] == '1':
            one_one_value += 1
        count += 1

    return value / count, one_one_value / count


with open('extractVars_onefam_DBVDC.vcf', "r") as file:
    for line in file:
        if not line.startswith('#'):
            #print(line)
            text = line.split("\t")
            selection = text[7].split(";")
            #print("nl\n", text[9:])
            het_val, pure_het_value = calculate_heterozygosity(text[9:])
            for el in selection:
                if el[0:2] == "AF":
                    af = el[3:]
            dict[text[1]] = af, het_val, pure_het_value
    print("extract end")

with open("All_chromosomes_one_fam.vcf", "r") as file:
    for line in file:
        text = line.split("\t")
        chrm = text[0]
        position = text[1]
        gene = text[19]
        impact = text[18]
        gene_dict[position] = gene, impact, chrm
    print("read chr end")
with open("af_filtered_All_chr_one_fam.vcf", "w") as file:
    file.write("chromosome\tposition\tallele_frequency\tgene_name\timpact\thomozygosity_value\tpure_homozygosity_value\n")
    for key in dict.keys():
        if key in gene_dict.keys():
            file.write(gene_dict[key][2] + "\t" + key + "\t" + dict[key][0] + "\t" + gene_dict[key][0] + "\t" + gene_dict[key][1] + "\t"
                       + str(dict[key][1]) + "\t" + str(dict[key][2]) + "\n")
