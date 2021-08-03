# Python 3
'''
# Plot HLA motif from 9mer cores
##Will generate and plot position specific frequencies of amino acids in 9mer peptides

-This will take text file containing list of peptides in one column. Specify this file at start of script

-separates out the aminoacids into position specific columns (Pos- 1-9) into a dataframe

-calculates the number of amino acids at each position and its percentage.
-plots a heatmap

Example file DR_cores.txt has 9-mer minimal cores.
Output will be stores in folder 'output'

'''
import collections
import pandas as pd
import os
import seaborn as sns
import re
import matplotlib.pyplot as plt

#########################################################
# specify file name here or pass to the function at the end
# Alternatively can be modified to run on all .txt files in directory.
filename = 'DR_cores.txt'
#########################################################
# output files will be placed here
output_folder = './output/'
os.makedirs(output_folder, exist_ok=True)
#########################################################
# main function
def count_freq_motif(file):
    '''

    :param file: takes text files with one column containing 9-mer peptides
    :return: heatmap and AA frequency tables
    '''
    print('Now processing :', file)

    name = file[:-4]
    f2 = output_folder+ name + '_AA_Counts_and_perc.csv'
    f2 = open(f2, 'w', newline='')  # for windows have newline '' argument
    f3 = output_folder+ name + '_AA_table.csv' # each peptide split into AA per position
    f3 = open(f3, 'w', newline='')  # for windows have newline '' argument
    f4 = output_folder+ name + '_AA_Perc.csv'
    f4 = open(f4, 'w', newline='')  # for windows have newline '' argument
    # generate a pandas dataframe
    columns = ['Sequence', 'length', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    df = pd.DataFrame(columns=columns)
    # generate lists
    sequence = []
    length = []

    pos1 = []
    pos2 = []
    pos3 = []
    pos4 = []
    pos5 = []
    pos6 = []
    pos7 = []
    pos8 = []
    pos9 = []

    lines = [line.strip() for line in open(file)]
    for line in lines:
        # print(line)
        # print('position 1=' + line[0])
        # print('position 9=' + line[8])

        sequence.append(line)
        length.append(len(line))

        pos1.append(line[0])
        pos2.append(line[1])
        pos3.append(line[2])
        pos4.append(line[3])
        pos5.append(line[4])
        pos6.append(line[5])
        pos7.append(line[6])
        pos8.append(line[7])
        pos9.append(line[8])

        letter1 = collections.Counter(line)

    df['Sequence'] = sequence
    df['length'] = length

    df['1'] = pos1
    df['2'] = pos2
    df['3'] = pos3
    df['4'] = pos4
    df['5'] = pos5
    df['6'] = pos6
    df['7'] = pos7
    df['8'] = pos8
    df['9'] = pos9

    # make a new dataframe
    columns = ['Counts_1', 'Counts_2', 'Counts_3', 'Counts_4', 'Counts_5', 'Counts_6', 'Counts_7', 'Counts_8',
               'Counts_9']
    aa = ['G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T']
    df_freq = pd.DataFrame(index=aa, columns=columns)

    df_freq['Counts_1'] = df['1'].value_counts()
    df_freq['Counts_2'] = df['2'].value_counts()
    df_freq['Counts_3'] = df['3'].value_counts()
    df_freq['Counts_4'] = df['4'].value_counts()
    df_freq['Counts_5'] = df['5'].value_counts()
    df_freq['Counts_6'] = df['6'].value_counts()
    df_freq['Counts_7'] = df['7'].value_counts()
    df_freq['Counts_8'] = df['8'].value_counts()
    df_freq['Counts_9'] = df['9'].value_counts()

    df_freq['Perc_1'] = (df_freq['Counts_1'] / df_freq['Counts_1'].sum()) * 100
    df_freq['Perc_2'] = (df_freq['Counts_2'] / df_freq['Counts_2'].sum()) * 100
    df_freq['Perc_3'] = (df_freq['Counts_3'] / df_freq['Counts_3'].sum()) * 100
    df_freq['Perc_4'] = (df_freq['Counts_4'] / df_freq['Counts_4'].sum()) * 100
    df_freq['Perc_5'] = (df_freq['Counts_5'] / df_freq['Counts_5'].sum()) * 100
    df_freq['Perc_6'] = (df_freq['Counts_6'] / df_freq['Counts_6'].sum()) * 100
    df_freq['Perc_7'] = (df_freq['Counts_7'] / df_freq['Counts_7'].sum()) * 100
    df_freq['Perc_8'] = (df_freq['Counts_8'] / df_freq['Counts_8'].sum()) * 100
    df_freq['Perc_9'] = (df_freq['Counts_9'] / df_freq['Counts_9'].sum()) * 100
    
    #export files
    df_freq.to_csv(f2)
    df.to_csv(f3, index=False)
    # print(Freq.columns)

    # make another df with only percentages
    df_perc = df_freq.filter(regex=r'^Perc', axis=1)
    df_perc = df_perc.round(2)  # rounds everything to 2 decimals.
    df_perc.to_csv(f4)

    # change column names for better labels
    df_perc = df_perc.rename(columns=lambda x: re.sub('^Perc_', '', x))

    sns.set()

    # Draw a heatmap with the numeric values in each cell

    freq_heatmap = sns.heatmap(data=df_perc, annot=True, linewidths=.5, vmax=100)

    plt.yticks(rotation=0)  # to turn Y labels on heatmap to horizontal
    # plt.show() # shows plot
    plt.savefig(output_folder+ name + '_heatmap' +'.svg', format='svg', dpi=600)
    plt.gcf().clear()


count_freq_motif(filename)
