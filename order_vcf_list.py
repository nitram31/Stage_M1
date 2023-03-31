import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("vcf_list", help="put the file that you want to order here, note that it will only work on dogs")
args = parser.parse_args()

order = {'chr1': '', 'chr2': '', 'chr3': '', 'chr4': '', 'chr5': '', 'chr6': '', 'chr7': '', 'chr8': '', 'chr9': '',
         'chr10': '', 'chr11': '', 'chr12': '', 'chr13': '', 'chr14': '', 'chr15': '', 'chr16': '', 'chr17': '',
         'chr18': '', 'chr19': '', 'chr20': '', 'chr21': '', 'chr22': '', 'chr23': '', 'chr24': '', 'chr25': '',
         'chr26': '', 'chr27': '', 'chr28': '', 'chr29': '', 'chr30': '', 'chr31': '', 'chr32': '', 'chr33': '',
         'chr34': '', 'chr35': '', 'chr36': '', 'chr37': '', 'chr38': '', 'chrX': '', 'chrY': '', "chrM": '',
         'chrY_unplaced': '', 'chrUn': ''}
content = []
with open(args.vcf_list, "r") as file:
    f = file.read()
    print(f)
    for line in f.split('\n'):
        content.append(line.split('\n')[0])

common_part = os.path.commonprefix([content[0], content[-1]])

for path in content:
    i = 0
    suffix = ''
    while path[len(common_part) + i] != '_':
        suffix += path[len(common_part) + i]
        i += 1

    if suffix == 'Y' and (path[len(common_part) + i + 1]) == 'u':
        suffix = 'Y_unplaced'
    order['chr' + suffix] += path+"\n"

with open(args.vcf_list, "w") as file:
    for key in order.keys():
        for path in order[key]:
            file.write(path)

