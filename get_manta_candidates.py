import argparse
import os

""" This script is used to parse the raw manta data and to add the corresponding gene to each location. It takes an 
"expected genotype" file, where the dogs ID and its genotype should be organized as follows :

Dog_id  expected_genotype_1 expected_genotype_2  ... (the header mustn't be in the file)
6958RE  1/1 0/1 1/0
BC0480	0/1	0/0	1/0
BC547	1/1

supported format are csv (comma or semicolon) and tsv.

command example :

python3 get_manta_candidate -g dog_genotype_path -i manta_file_path -c chrX -a annotation_file_path -o output_name

The second argument is the manta file containing the structural variants as well as the genotypes.

The third (facultative) argument is the specific chromosome(s), it can be all the chromosome separated by a space or it 
can be an interval. Example: chrX chr3 chr1-15

The fourth (facultative) argument is the file containing the path to the gtf annotation files that allows the affected 
genes to be added to the output 

The fifth (facultative) argument is the output name of the modified file

Another option is to use the -p or --parameter with the path to the parameter file that contains more complete options

python3 get_manta_candidate -p manta_parameter

The output is the original manta line with the gene_id and the location of the variant added. In the case of a 
translocation the gene_id is complemented with the gene_id of the gene the translocation is located  

CHROM	POS	GENE_ID	LOCATION	ID	REF	ALT	
chr1	17938	chr_start-ENPP1_3-(LOC111090445_3)	Intergenic	MantaBND:0:75537:75538:0:7:0:1	A	A[CHR17:53307608[

In this case, the translocation is located between the start of the chromosome and a gene (ENPP1_3) and the fragment
translocated went in another gene (LOC111090445_3) located on chromosome 17

"""

parser = argparse.ArgumentParser()
parser.add_argument('-g', "--genotype", help='path to the dog expected dog genotype file')
parser.add_argument("-i", "--input", help="path to the raw manta file")
parser.add_argument("-a", "--annotation",
                    help="path to the file containing the path of all the gtf annotation files")
parser.add_argument("-c", "--chromosome", action='append', nargs='+',
                    help="facultative, used when only one or more chromosome is wanted")
parser.add_argument("-gl", "--gene_list", action='append', nargs='+',
                    help="facultative, used when only one or more gene is wanted")
parser.add_argument("-s", "--start_end", action='append', nargs='+',
                    help="facultative, used when only one portion of a chromosome is wanted (-c option must be used)")
parser.add_argument("-l", "--location", action='append', nargs='+',
                    help="facultative, used when only one or more variant location is wanted")
parser.add_argument("-q", "--quality", action='append', nargs='+',
                    help="facultative, used when only one or more variant location is wanted")
parser.add_argument("-p", "--parameter", help="path to the parameter file")
parser.add_argument("-o", "--output", help="facultative, name of the output file")

args = parser.parse_args()
ann_dict = {}


def parse_parameter_file(path):
    with open(path, "r") as file:
        cpt = 0
        gene_list = []
        location_list = []
        for line in file:
            if not line.isspace() and '#' not in line and line != '\n':
                match cpt:
                    case 0:  # path to the manta file
                        m_path = line[11:-1]
                    case 1:  # path to the annotation file
                        a_path = line[16:-1]
                    case 2:  # path to the genotype file
                        g_path = line[5:-1]
                    case 3:  # name of the output
                        o_path = line[12:-1]
                    case 4:  # list of genes that will only be included in the output file
                        if line[7:-2] != '':
                            gene_list = line[7:-2].split(" ")
                        else:
                            gene_list = []
                    case 5:  # chromosome list that will be used as a filter in the output
                        if line[13:-2] != '':
                            c_list = line[13:-2].split(" ")
                            c_list = parse_chromosome(c_list)
                        else:
                            c_list = []
                    case 6 | 7:  # start-end filter
                        if line.startswith("start='"):
                            start = line[7:-2]
                        else:
                            end = line[5:-2]
                    case 8:  # Location filter
                        if line[10:-2] == 'ALL':
                            location_list = []  # 'Intronic', 'Exonic', 'Promoter' or 'Intergenic'
                        else:
                            location_list = line[10:-2].split(" ")
                    case 9:  # quality filter
                        if line[16:-2] == 'ALL':
                            quality_filter_list = []
                        else:
                            quality_filter_list = line[16:-2].split(" ")
                cpt += 1
        return m_path, a_path, g_path, o_path, gene_list, c_list, start, end, location_list, quality_filter_list


def parse_chromosome(chromosome_list):
    chr_list = chromosome_list
    for el in chr_list:
        if '-' in el:
            chr_pos = el[3:].split("-")
            chr_list.pop(chr_list.index(el))
            for i in range(int(chr_pos[0]), int(chr_pos[1]) + 1):
                pos = "chr" + str(i)
                if pos not in chr_list:
                    chr_list.append(pos)
    return chr_list


def define_separator(path):
    dog_dict = {}
    with open(path, "r") as file:
        text = file.read().split("\n")
        sep = ""
        for char in text[0]:
            match char:
                case "\t":
                    sep = "\t"
                    break
                case ";":
                    sep = ";"
                    break
                case ",":
                    sep = ","
                    break
                case _:
                    pass
        if sep == "":
            raise Exception("Invalid separator, please use csv or tsv files")

        for line in text:
            current_dog = line.split("\t")[0]
            if current_dog not in dog_dict.keys():
                dog_dict[current_dog] = []
            for genotype in line.split("\t")[1:]:
                dog_dict[current_dog].append(genotype.split("\n")[0])
    return sep, dog_dict


def read_manta(path, dog_dict, sep):
    """
    Used to read the manta file and to filter the variant that are coherent based on the dogs' genotype.

    :param path: path to the manta file
    :param dog_dict: dictionary containing the genotype of the dogs
    :param sep: separator of the manta file
    :return: body : the lines of the manta file, file_name : base name of the manta file
    """
    body = ""
    dog_list = []

    def line_analysis(li, cont):

        chromosome = cont[0]
        if 'chrY' in chromosome or 'chrUn' in chromosome:
            return 'unannotated_chr'
        position = cont[1]
        location, gene_id = get_corresponding_info(position, chromosome)

        if '[' in cont[4] or ']' in cont[4]:
            gene_id += '-' + get_corresponding_chromosome(cont[4])
        line_split = li.split(sep)
        line_split.insert(2, gene_id)
        line_split.insert(3, location)

        li = list_to_sep_line(line_split, sep)
        print(li)
        if args.chromosome and chromosome not in args.chromosome:
            # filter based on chromosome
            return ''
        if args.gene_list and gene_id[-2] not in args.gene_list:
            return ''
        if args.start_end and not args.start_end[0] < int(position) < args.start_end[1]:
            return ''
        if args.location and not test_location(location):
            return ''
        if args.quality and content[5] not in args.quality:
            return ''
        return li

    with open(path, "r") as file:

        text = file.read().split("\n")
        first_line = text[0]
        text = text[1:]  # deleting the header for the main part of the text
        first_line_split = first_line.split("\t")

        for j in range(14, len(first_line_split)):  # getting the dogs name from the manta file header
            dog_list.append(first_line_split[j])

        if args.annotation:  # adding new columns to the header if --annotation is specified
            first_line_split.insert(2, 'GENE_ID')
            first_line_split.insert(3, 'LOCATION')
            first_line = list_to_sep_line(first_line_split, sep)

        body += first_line + "\n"
        i = 0
        for line in text[0:-1]:
            content = line.split("\t")
            valid_line = ''
            if args.genotype:
                for j in range(14, len(content)):  # tests if the genotype for the variant correspond to those expected
                    if content[j] not in dog_dict[dog_list[j - 14]]:
                        break

                    if j == len(content) - 1:
                        i += 1
                        valid_line = line_analysis(line, content)
            else:
                valid_line = line_analysis(line, content)
            if valid_line != '':
                body += valid_line + "\n"
        print(i)
    return body


def test_location(loc):
    """

    :param loc: Str, correponds to the location of the variant
    :return: boolean, true if the variant location corresponds to the one specified in options
    """
    if 'exon' in loc:
        loc = 'Exonic'
    if loc in args.location:
        return True
    return False


def read_annotation_data(annotation_file):
    """
    reads the annotation data and store it into ann_dict dictionary.

    :param annotation_file: gtf format file that is parsed to get the annotation data
    """
    text_file = open(annotation_file, 'r')
    text = text_file.read().split('\n')
    for line in text:
        if line != '':
            with open(line, "r") as file:
                for lines in file:
                    get_info_from_gff(lines)


def get_info_from_gff(line):
    """
    used to extract the name and exon location of a gene, that will be stored in a dictionary.

    :param line: one line of annotation that will be treated
    """

    if "#" not in line:
        content = line.split('\t')
        description = content[8].split(";")
        gene_id = description[0][9:-1]
        chromosome_id = content[0]
        gene_start = content[3]
        gene_end = content[4]
        if content[2] == 'gene':

            if chromosome_id not in ann_dict:
                ann_dict[chromosome_id] = {}
                ann_dict[chromosome_id][gene_id] = {'coord': [gene_start, gene_end], 'exon': {}}

            else:
                ann_dict[chromosome_id][gene_id] = {'coord': [gene_start, gene_end], 'exon': {}}

        elif content[2] == 'exon':  # searching for the exon_number information that has a variable location
            exon_number_position = [j for j, word in enumerate(description) if word.startswith(' exon_number "')][0]
            exon_number = "exon " + str(description[exon_number_position][14:-1])
            ann_dict[chromosome_id][gene_id]['exon'][exon_number] = [gene_start, gene_end]


def list_to_sep_line(str_list, sep):
    """
    used to reconstitute the manta file lines after inserting new columns.

    :param str_list: any kind of list
    :param sep: separator of the file you will write the string to
    :return: sep_line : string with the elements of the list separated
    """
    sep_line = str_list[0]
    for element in str_list[1:]:
        sep_line += sep + str(element)
    return sep_line


def get_corresponding_chromosome(content):
    """
    Gets the chromosome and the location from a string, for example : "G]CHR29:12369845", in this case,
    chromosome will take the value chr29 and other_chromosome_location will take the value 12369845.

    :param content: manta generated string corresponding to the variant linked chromosome
    :return: the hypothetical gene_id corresponding to the location
    """

    content = content.split(':')

    if 'CHRUN' in content[0] or 'CHRY' in content[0]:
        return ''

    other_chromosome_location = content[1][0:[j for j, number in enumerate(content[1]) if number.isdigit()][-1] + 1]

    chromosome_id = [j for j, number in enumerate(content[0]) if number.isdigit() or number == 'X']
    chromosome = 'chr' + content[0][chromosome_id[0]:chromosome_id[-1] + 1]

    gene_id = get_corresponding_info(other_chromosome_location, chromosome)[1]
    if gene_id == 'none':
        return ''
    return "(" + gene_id + ")"


def get_corresponding_info(var_pos, chrom):
    """
    Used to search in ann_dict dictionary to know if the location of the variant lands in a gene and if it does, in what
    exon if any.

    :param var_pos: int of the position on a given chromosome to get the hypothetical corresponding gene
    :param chrom: name of a chromosome
    :return: exon_number: number of the exon if the mutation is in one, gene_id : id of the gene
    """
    previous_gene = ''
    previous_gene_location = [0, 0]
    cpt = 0
    for gene_id in ann_dict[chrom]:
        gene_location = ann_dict[chrom][gene_id]['coord']
        #  if the variant is located within a gene
        if int(gene_location[0]) < int(var_pos) < int(gene_location[1]) \
                or int(gene_location[1]) < int(var_pos) < int(gene_location[0]):

            for exon_number in ann_dict[chrom][gene_id]['exon']:
                exon_location = ann_dict[chrom][gene_id]['exon'][exon_number]
                if int(exon_location[0]) < int(var_pos) < int(exon_location[1]) \
                        or int(exon_location[1]) < int(var_pos) < int(exon_location[0]):
                    return exon_number, gene_id
            exon_number = 'Intronic'
            return exon_number, gene_id

        #  if the variant is located between two genes

        elif int(previous_gene_location[1]) < int(var_pos) < int(gene_location[0]) and cpt > 0:
            if int(gene_location[0]) - 1024 < int(var_pos) < int(gene_location[0]):
                return 'promoter', gene_id
            else:
                return 'Intergenic', previous_gene + '-' + gene_id

        previous_gene_location = gene_location
        previous_gene = gene_id
        cpt += 1

    return 'none', 'none'


def format_location(loc_list):
    """

    :param loc_list: list of str
    :return: list of str with the string reformated
    """
    for i in range(0, len(loc_list)):
        loc_list[i] = loc_list[i][0] + loc_list[i][1:].lower()

    return loc_list


def write_to_output(body, output_name):
    with open(output_name, "w") as file:
        file.write(body)


def main():
    if args.parameter:
        arg = parse_parameter_file(args.parameter)

        manta_path = arg[0]
        args.input = arg[0]

        annotation = arg[1]
        args.annotation = arg[1]

        dog_genotype_path = arg[2]
        if dog_genotype_path:
            args.genotype = dog_genotype_path

        output_name = arg[3]

        if output_name == '':
            file_name = os.path.basename(manta_path)
            file_name = os.path.splitext(file_name)[0]
            output_name = file_name + "_structural_candidates.txt"

        gene_list = arg[4]
        if gene_list:
            args.gene_list = gene_list

        chr_list = arg[5]
        if chr_list:
            args.chromosome = chr_list

        start = arg[6]
        end = arg[7]
        if start and end:
            args.start_end = int(start) * 1000000, int(end) * 1000000

        location_list = arg[8]
        if location_list:
            location_list = format_location(location_list)
            args.location = location_list

        quality_filter_list = arg[9]
        if quality_filter_list:
            args.quality = quality_filter_list

        if args.start_end and not args.chromosome:
            raise RuntimeError("No chromosome found, please specify a chromosome to search the region in")
        if args.start_end and len(args.chromosome) != 1:
            raise RuntimeError('Found ' + str(len(args.chromosome)) + " chromosome, expected 1. Please give only one"
                               + " chromosome when searching for a specific region")

    else:

        dog_genotype_path = args.genotype

        manta_path = args.input
        annotation = args.annotation
        if args.gene_list:
            gene_list = args.gene_list[0]
            args.gene_list = gene_list

        if args.start_end:
            start, end = args.start_end[0]
            args.start_end = [start * 1000000, end * 1000000]

        if args.location:
            location_list = args.location[0]
            args.location = location_list

        if args.quality:
            quality_filter_list = args.quality[0]
            args.quality = quality_filter_list

        file_name = os.path.basename(manta_path)
        file_name = os.path.splitext(file_name)[0]
        output_name = file_name + "_structural_candidates.txt"

        if args.output:
            output_name = args.output

        if args.chromosome:
            args.chromosome = args.chromosome[0]
            args.chromosome = parse_chromosome(args.chromosome)

        if args.start_end and not args.chromosome:
            raise Exception("No chromosome found, please specify a chromosome to search the region in")
        if args.start_end and len(args.chromosome) != 1:
            raise RuntimeError('Found ' + str(len(args.chromosome)) + " chromosome, expected 1.")
    separator, dog_genotype_dict = "\t", {}
    if args.genotype:
        separator, dog_genotype_dict = define_separator(dog_genotype_path)

    read_annotation_data(annotation)
    output_body = read_manta(manta_path, dog_genotype_dict, separator)
    write_to_output(output_body, output_name)


if __name__ == "__main__":
    main()
