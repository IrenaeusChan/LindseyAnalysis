import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
pd.options.mode.chained_assignment = None  # default='warn'

#Handle the various arguments
parser = argparse.ArgumentParser(description='filter merged fastq reads to those of appropiate length')
parser.add_argument("-s1", "--string1", type=str, required=True)
parser.add_argument("-s2", "--string2", type=str, required=True)
parser.add_argument("-l", "--library", type=str, required=True)
parser.add_argument("--sequence", type=str, required=True)
print('start tabulation')

args = parser.parse_args()
#Define variables for libraries
codon = 0
WT = args.sequence
# if args.library == 'Library_1':
#     WT = 'EGVNWLAQRFHCQNGCILGDEMGLGKTCQTIALFIYLAGRLNDEGPFLIL'
#     codon = 0
# elif args.library == '':
#     WT = ''
#     codon = 0
# elif args.library == '':
#     WT = ''
#     codon = 0
# else:
#     WT = ''
#     codon = 0

#Functions
#ratioizer takes the ratio between RCPM of the selection/no-selection datasets and takes then takes the log10 value.
def ratioizer(table1, table2):
    output_table = table1 / table2
    output_table_log = output_table.apply(lambda x: np.log10(x) if np.issubdtype(x.dtype, np.number) else x)
    return output_table_log

#subtractor will subtract the mutant log10 ratio by the WT log10 ratio for each position on the protein.
def Substractor(table, codon):
    for x in WT:
        codon += 1
        row = table.loc[table["Amino Acid"] == x].index[0]
        WT_value = table.loc[row, "codon{}".format(codon)]
        table["codon{}".format(codon)] = table["codon{}".format(codon)].apply(lambda x: x - WT_value)
    return(table)

Length1 =0
Length2 =0


#read in the two data sets (table 1 = selection condition; table 2 = no-selection condition)
table1 = pd.read_csv("{}.csv".format(args.string1), keep_default_na=False)

table1_num = table1.drop('Amino_Acid', axis=1)
table2 = pd.read_csv("{}.csv".format(args.string2), keep_default_na=False)
table2_num = table2.drop('Amino_Acid', axis=1)

#Call the ratioizer function for the two tables
table_log = ratioizer(table1_num, table2_num)

#Insert a column at position 0 in the table
AminoAcid =('A', 'R', 'N', 'D', 'C', 'E', 'Q', \
          'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T','W', 'Y','V', "*")
AminoAci_new =('*', 'A', 'C', 'D', 'E', 'F', \
          'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R','S', 'T','V', "W", "X","Y")

table_log.insert(0, 'Amino Acid', AminoAci_new)

#call the subtractor function on the log10 ratio table
table_log_dif = Substractor(table_log, codon)

#export final table
table_log_dif_removeX = table_log_dif.drop(axis=0, index=20)
table_log_dif_removeX.to_csv('diflog_{}_{}.csv'.format(args.string1,args.string2), index=False)
