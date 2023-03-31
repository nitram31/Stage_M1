import os

path = "SEPT4.vcf"
with open(path, "r") as file:
    body = "#CHROM\tPOS\tREF\tALT\tFILTER\t12675RE\t6958RE\tBC0480\tBC1029\tBC547\tBC548\tBC549\tB_10764\tB_13693\tB_7558\tB_8345\tALLELE\tEFFECT\tIMPACT\tGENE\tGENEID\tFEATURE\tFEATUREID\tBIOTYPE\tRANK\tHGVS_C\tHGVS_P\n"

    file_name = os.path.basename(path)
    file_name = os.path.splitext(file_name)[0]
    for line in file:
        text = line.split("\t")
        RE_1267 = text[5]
        RE_6958 = text[6]
        BC0480 = text[7]
        BC1029 = text[8]
        BC547 = text[9]
        BC548 = text[10]
        BC549 = text[11]
        B_10764 = text[12]
        B_13693 = text[13]
        B_7558 = text[14]
        B_8345 = text[15]
        if RE_1267 == "1/1":
            if RE_6958 == "1/1" or RE_6958 == "0/1" or RE_6958 == "1/0":
                if BC0480 == "0/1" or BC0480 == "0/0" or BC0480 == "1/0":
                    if BC1029 == "0/0" or BC1029 == "1/0" or BC1029 == "0/1":
                        if BC547 == "1/1" and BC548 == "1/1" and BC549 == "1/1":
                            if B_10764 == "1/1":
                                if B_13693 == "1/0" or B_13693 == "0/1":
                                    if B_7558 == "1/1":
                                        if B_8345 == "1/1":
                                            body += line
with open(file_name + "_candidates.vcf", "w") as file:
    file.write(body)

