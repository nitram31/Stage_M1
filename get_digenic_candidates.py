import argparse

"""
This script is used to parse the vcf file containing the variant in the population to find two genes that are only 
present as homozygous in affected subjects. 

First a .vcf file must be created, containing only the header and the lines that contains the variants of interest. 

Then "expected genotype" file has to be created, where the dogs ID and its 
genotype should be organized as follows, with one line for each expected genotype:

Note that the dogs id should be identical to the one contained in the header of the first .vcf file

Dog_id  expected_genotype_for_the_second_gene
6958RE  1/1
6958RE  1/0
BC548   1/1   
BC547   1/1



supported format is tsv.

The last file is the .vcf file containing the variant and the genotype of the dogs in the population, note that a header 
with the identical dog_id as the other two files is required. 

command example :

python3 get_digenic_candidates -g dog_genotype_path -k known_variant_path -p population_variant_path -o output_name


"""
# sys.tracebacklimit = -1

parser = argparse.ArgumentParser()
parser.add_argument('-g', "--genotype", help='path to the expected dog genotype file')
parser.add_argument('-k', '--known_variant', help='path to the file containing the manually filtered variant of the \
 first gene that will be matched with others ')
parser.add_argument('-p', '--population_variant', help='path to the file containing the genotype of the population')
parser.add_argument('-a', '--annotation', help='path to the file containing the variant for which the dogs are \
                                                homozygous')
parser.add_argument('-o', '--output', help='path to output file', default='digenic_candidates')
args = parser.parse_args()


def read_genotype(path):
    """

    :param path: path to the genotype file
    :return: sep ; separator used in the genotype file, dog_dict ; dictionary containing the dogs id and genotype
    """
    dog_dict = {}
    with open(path, "r") as file:
        text = file.read().split("\n")[1:]
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
            raise Exception("Invalid separator, please use tsv files")
        for line in text:
            current_dog = line.split(sep)[0]
            if current_dog not in dog_dict.keys():
                dog_dict[current_dog] = {'genotype': ''}

            if line != '':
                genotype = line.split(sep)[1]
            dog_dict[current_dog]['genotype'] = genotype
    return sep, dog_dict


def read_known_variants(path):
    """

    :param path: path to the vcf file containing the variant that are suspected to be linked with another gene
    :return: known_variant_dict, dictionary containing chromosome, location and the full line
    """
    known_variant_dict = {}
    with open(path, 'r') as file:
        text = file.read().split('\n')
        for line in text[1:]:
            if line and line[0] != '':
                line = line.split('\t')
                chr = line[0]
                position = line[1]
                if chr not in known_variant_dict.keys():
                    known_variant_dict[chr] = [position], [list_to_sep_line(line, '\t')]
                else:
                    known_variant_dict[chr][0].append(position)
                    known_variant_dict[chr][1].append(list_to_sep_line(line, '\t'))

    return known_variant_dict


def test_homozyg(string):
    """
    Used to test if the string contains a homozygous marker
    :param string: str '1/1' or '0|1'
    :return: boolean, true if the string contains a homozygous marker
    """
    genotype = string[0:3]
    if genotype == "1/1" or genotype == "1|1":
        return True
    return False


def read_population_variant(path):
    """
    Used to read the  homozygosity of the population for all the different variants
    :param path: path to the population variant file
    :return: a dictionnary containing the chromosome, location and list of homozygous dogs for this position and the
    size of the population for stats purposes
    """
    dog_id_list = []

    with open(path, 'r') as file:
        text = file.read().split('\n')
        header = text[0].split('\t')
        text = text[1:]

        for dog in header[9:]:
            dog_id_list.append(dog)

        variant_dict = {}

        for line in text:
            line = line.split('\t')

            if line and line[0] != '':
                chromosome = line[0]
                location = line[1]

                if chromosome not in list(variant_dict.keys()):
                    variant_dict[chromosome] = {location: {'dog_genotype': []}}

                else:
                    variant_dict[chromosome][location] = {'dog_genotype': []}

                dog_genotype = variant_dict[chromosome][location]['dog_genotype']
                cpt = 0
                population_size = len(line[9:])
                for genotype in line[9:]:
                    if test_homozyg(genotype):
                        dog_genotype.append(dog_id_list[cpt])
                    cpt += 1

    return variant_dict, population_size


def get_variant_of_interest(variant_dict, dog_dict, known_variant_dict):
    """

    :param variant_dict: dictionary of the variant known in the population
    :param dog_dict: dictionary of the dogs ids and their expected genotype
    :param known_variant_dict: dictionary of the variants known
    :return: variant_dict modified to include the variant that are linked to the known variants
    """
    for known_chromosome in known_variant_dict:

        for known_position in known_variant_dict[known_chromosome][0]:

            for chromosome in variant_dict:

                variant_position = variant_dict[chromosome]
                for variant in variant_position:
                    if 'number_digenic_homozygous' not in variant_position[variant]:
                        variant_position[variant][
                            'number_digenic_homozygous_' + known_chromosome + "_" + known_position] = 0

                    if set(list(dog_dict.keys())).issubset(variant_position[variant]['dog_genotype']):
                        for dog in variant_position[variant]['dog_genotype']:
                            if dog in variant_dict[known_chromosome][known_position]['dog_genotype']:
                                variant_position[variant]['number_digenic_homozygous_'
                                                          + known_chromosome
                                                          + "_"
                                                          + known_position] += 1
                    else:
                        variant_position[variant]['number_digenic_homozygous_'
                                                  + known_chromosome
                                                  + "_"
                                                  + known_position] = -1
    return variant_dict


def output_variants(annotation_path, variant_dict, known_variant_dict):
    with open(annotation_path, 'r') as file:
        text = file.read().split('\n')[:-1]
        header = text[0].split('\t')
        header.insert(5, '#HOMOZYGOUS_INDIVIDUALS')
        header = list_to_sep_line(header, '\t')
        body = header + '\n'
        for known_chromosome in known_variant_dict:
            for known_position in known_variant_dict[known_chromosome][0]:
                line_index = known_variant_dict[known_chromosome][0].index(known_position)
                body += ">" + known_variant_dict[known_chromosome][1][line_index] + "\n"
                for line in text[1:]:
                    line = line.split('\t')
                    chromosome = line[0]    
                    position = line[1]
                    if chromosome in variant_dict and position in variant_dict[chromosome]:
                        line.insert(5, variant_dict[chromosome][position]['number_digenic_homozygous_'
                                                                          + known_chromosome
                                                                          + "_"
                                                                          + known_position])
                        body += list_to_sep_line(line, '\t') + "\n"
    with open(args.output, 'w') as file:
        file.write(body)


def list_to_sep_line(str_list, sep):
    """
    used to reconstitute the file lines after inserting new columns.

    :param str_list: any kind of list
    :param sep: separator of the file you will write the string to
    :return: sep_line : string with the elements of the list separated
    """
    sep_line = str_list[0]
    for element in str_list[1:]:
        sep_line += sep + str(element)
    return sep_line


def main():
    annotation_path = args.annotation
    genotype_path = args.genotype
    known_variant_path = args.known_variant
    population_variant_path = args.population_variant
    sep, dog_dict = read_genotype(genotype_path)
    known_variant_dict = read_known_variants(known_variant_path)
    variant_dict = read_population_variant(population_variant_path)[0]
    variant_dict = get_variant_of_interest(variant_dict, dog_dict, known_variant_dict)
    output_variants(annotation_path, variant_dict, known_variant_dict)


if __name__ == '__main__':
    main()
