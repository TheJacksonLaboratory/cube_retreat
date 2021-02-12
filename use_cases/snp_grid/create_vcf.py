"""
    Copyright 2020 The Jackson Laboratory.  All Rights Reserved.
"""
import glob
import os.path
from os import path
import logging


def run(argv=None):
    """
    This is spark job, and need to install the spark first
        https://spark.apache.org/downloads.html

    It take three csv files and generate a pusedo-vcf.
    to do: change to take argments, merge the header
    """
    logging.info("begin")

    # data directory
    input_dir = "csv/"
    output_dir = "vcf/"
    chr = "chr1"
    data_start_index = 6

    merged_files = glob.glob(input_dir + "merged_*")
    imputed_files = []
    for mg_file in merged_files:
        logging.info(mg_file)
        index = mg_file.rfind('_');
        file_name = mg_file[index:]
        im_file = input_dir + "imputed_prob_" + chr + file_name
        logging.info(im_file)
        if path.exists(im_file):
            logging.info(im_file + ": exist")
        else:
            logging.error(im_file + ": DO NOT exist")

        region = file_name[1:file_name.rfind('.')]
        logging.info(output_dir)

        vcf_file_name = output_dir + region + ".vcf"
        merge_two_file(im_file, mg_file, data_start_index, vcf_file_name)

def merge_two_file(imputed_file_name, merged_file_name, data_start_index, vcf_file_name):
    logging.info("merge_two_file")
    logging.info("imputed_file_name: " + imputed_file_name)
    logging.info("merged_file_name: " + merged_file_name)
    logging.info("vcf_file_name: " + vcf_file_name)
    tab = "\t"
   
    with open(vcf_file_name, 'w') as vcf_out:
        with open('vcf_prefix.txt') as prefix:
            vcf_out.write(prefix.read())

        with open(imputed_file_name, 'r') as imputed_file:
            with open(merged_file_name, 'r') as merged_file:
                imputed_first = imputed_file.readline()
                merged_first = merged_file.readline()
                imputed_firsts = remove_and_split(imputed_first)
                merged_firsts = remove_and_split(merged_first)

                num_of_columns = len(imputed_firsts)
                if not is_header_okay(imputed_firsts, merged_firsts):
                    logging.error("Two file header does not match")
                    return

                vcf_first = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
                for x in range(data_start_index, num_of_columns):
                    vcf_first += tab + imputed_firsts[x]
                vcf_out.writelines(vcf_first)

                for imputed_line, merged_line in zip(imputed_file, merged_file):
                    imputed_cells = remove_and_split(imputed_line)
                    merged_cells = remove_and_split(merged_line)
                    
                    if not is_data_okay(imputed_cells, merged_cells):
                        logging.error("Two row does not match")
                        logging.error(imputed_line)
                        logging.error(merged_line)
                        return

                    vcf_line = imputed_cells[0]         # CHROM
                    vcf_line += tab + imputed_cells[1]  # POS
                    vcf_line += tab + imputed_cells[2]  # ID
                    vcf_line += tab + imputed_cells[3]  # REF
                    vcf_line += tab                     # ALT
                    vcf_line += tab                     # QUAL
                    vcf_line += tab                     # FILTER
                    vcf_line += tab                     # INFO
                    vcf_line += tab + "GT:VA"           # FORMAT
                    for x in range(data_start_index, num_of_columns):
                        vcf_line += tab + get_value(imputed_cells, x) + ":" + get_value(merged_cells, x)
                    vcf_out.write(vcf_line + "\n")


def remove_and_split(str):
    out = str.replace("\"", "")
    return out.split(",")

def is_header_okay(list1, list2):
    if len(list1) != len(list2):
        return False
    for x in range(0, len(list1)):
        if list1[x] != list2[x]:
            return False
    return True

def is_data_okay(list1, list2):
    if len(list1) != len(list2):
        return False   
    if list1[0] != list2[0]:
            return False
    if list1[1] != list2[1]:
            return False
    return True

def get_value(list, index):
    value = list[index]
    if not value:
        return value
    return value.rstrip('\r\n')

if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    run()