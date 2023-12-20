# usage: python3.10 guide_find_SpCas9.py exon.fasta
"""
This program get all the possible guide sequences for an input exon sequence.
The output is two txt files in the file named potential_guide,
one for the sense strand and the other for the antisense strand.
"""

import os
import sys
import shutil
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq



current_directory = os.getcwd()


# get potential guide sequences
def get_guide_sequence(exon_string):
    potential_guide_list = []
    for i in range(len(exon_string) - 22):  
        substring = exon_string[i:i+23]           
        if substring[-2:] == 'GG':
            tmp_list = []
            
            tmp_list.append(i+1) # start
            tmp_list.append(substring[0:20]) # sequence
            tmp_list.append(substring[20:23]) # PAM
            # tmp_list.append(i+20) # end

            potential_guide_list.append(tmp_list)

            tmp_list = []

    return potential_guide_list






"""""""""""""""""""""""""""""""""""""""""""""""""""
"*** *** *** ***"     doing the first half job
"""""""""""""""""""""""""""""""""""""""""""""""""""
# get the exon sequence
print('\n',end='')
exon_string = str(SeqIO.read(sys.argv[1], 'fasta').seq)

# get guides            
guide_list = get_guide_sequence(exon_string)         
            
"""            
deal with the results
make a directory for storing result txt files, if the folder has alreadly existed, first
delete, then make the directory.
"""
foldername = 'potential_guide'
dirs = os.listdir(current_directory)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)

# go to the off_target folder
os.chdir(foldername)

guide_file = open("guide_in_sense_strand.txt",'w')

# write results
print('\r\n',end='',file=guide_file)
print("start".ljust(15,' '),\
      "sequence".ljust(30,' '),\
      "end".ljust(15,' '),\
      "PAM".ljust(15,' '),\
      sep="",\
      file=guide_file
      )
print('\r\n',end='',file=guide_file)


for i in range(len(guide_list)):
        C1_print_line_str = str(guide_list[i][0])
        C2_print_line_str = str(guide_list[i][1])
        
        C3_print_line_str = str(guide_list[i][0]+19)
        C4_print_line_str = str(guide_list[i][2])

            
        print(C1_print_line_str.ljust(15,' '),\
              C2_print_line_str.ljust(30,' '),\
              C3_print_line_str.ljust(15,' '),\
              C4_print_line_str.ljust(15,' '),\
              sep="",\
              file=guide_file
              )
             
guide_file.close()



"""""""""""""""""""""""""""""""""""""""""""""""""""
"*** *** *** ***"     doing the other half job
"""""""""""""""""""""""""""""""""""""""""""""""""""

os.chdir(current_directory)
exon_string = str(SeqIO.read('exon.fasta', 'fasta').seq.reverse_complement())

# get guides            
guide_list = get_guide_sequence(exon_string)     

# reverse the list
guide_list_sorted = sorted(guide_list,key=itemgetter(0),reverse=True) # in reverse order


os.chdir(foldername)
guide_file = open("guide_in_antisense_strand.txt",'w')

# write results
print('\r\n',end='',file=guide_file)
print("start".ljust(15,' '),\
      "sequence".ljust(30,' '),\
      "end".ljust(15,' '),\
      "PAM".ljust(15,' '),\
      sep="",\
      file=guide_file
      )
print('\r\n',end='',file=guide_file)


for i in range(len(guide_list)):
        C1_print_line_str = str(len(exon_string)+1-guide_list_sorted[i][0]) # position transformation
        C2_print_line_str = str(guide_list_sorted[i][1])             
        C3_print_line_str = str(len(exon_string)+1-guide_list_sorted[i][0]-19) # position transformation
        C4_print_line_str = str(guide_list_sorted[i][2])

            
        print(C1_print_line_str.ljust(15,' '),\
              C2_print_line_str.ljust(30,' '),\
              C3_print_line_str.ljust(15,' '),\
              C4_print_line_str.ljust(15,' '),\
              sep="",\
              file=guide_file
              )
             
guide_file.close()










