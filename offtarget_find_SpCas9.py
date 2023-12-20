# usage: python3.10 offtarget_find_SpCas9.py

"""
This program searches the potential off-targets in the chromosome sequences.
The chromosome fasta sequence and the input guide sequence are all in
uppercase DNAs. The output files are two txt files in the file named off_target,
one for the sense strand and the other for the antisense strand.
"""

import sys
import os
import shutil
import time
import random
from multiprocessing import Pool,Process,set_start_method,Queue
from operator import itemgetter
from Bio.Seq import Seq
from Bio import SeqIO


current_directory = os.getcwd()


# get the guide sequence from user input
print('\n',end='')
target_str = input("Please input the guide sequence (20 bases in uppercase DNA):")

# example: target_str = 'ATTGTAGTAGCCAGCGTTCC' # the 3rd guide in guide_in_antisense_strand.txt
if len(target_str) != 20:
        print("The length of the guide sequence is not 20.")
        print("Please check the guide sequence and rerun.")
        sys.exit()
        
if (set(list(target_str)).union({'A', 'T', 'C', 'G'})== {'A', 'T', 'C', 'G'}) == False:
        print("There are bases other than 'A', 'T', 'C', 'G' in the guide sequence.")
        print("Please check the guide sequence and rerun.")
        sys.exit()
        
print('\n',end='')        
print("The target sequence is %s." % target_str)
print('\n',end='')

 
"""
The string is the chromosome DNA sequence.
The m is the segment number.
The overlap_length is the number of overlap bases between the neighboring segments.
For SpCas9, the number is 20 (guide) + 3(NGG) - 1.
The function returns a list recording the segments with the start index (by Python list) and the segment string.
"""
def split_string(string, m, overlap_length):
    # make sure that the segment length is over 100
    n = len(string)
    if n < 100 * m:
 
        print("The segmentation number is too great.")
        print("Please decrease the segmentation (from 1000 to 100 or less, for example) number and rerun.")
        sys.exit()

         
    segment_length = n // m  # get the segment length
    
    segments = []
    position_and_segment = []
    
    for i in range(m - 1):
        if i == 0:
            start_index = 0
            end_index = segment_length
            segment = string[start_index:end_index]
            
            position_and_segment.append(0)
            position_and_segment.append(segment)

            segments.append(position_and_segment)

            position_and_segment = []
        else:
            start_index = i * segment_length
            start_index_new = start_index - overlap_length  # add necessary bases
            end_index = start_index + segment_length
            segment = string[start_index_new:end_index]
            
            position_and_segment.append(start_index_new)
            position_and_segment.append(segment)

            segments.append(position_and_segment)

            position_and_segment = []
            
    # deal with the last segment
    start_index = (m-1) * segment_length
    start_index_new = start_index - overlap_length  # add necessary bases
    end_index = n
    segment = string[start_index_new:end_index]

    position_and_segment.append(start_index_new)
    position_and_segment.append(segment)

    segments.append(position_and_segment)

    position_and_segment = []

    
    return segments



# search eligible off-target and store in the result queue
def find_off_target_positions(task_id,result_queue):   
    for chromosome_segment in chromosome_segments[task_id::sub_process_count]:
        start_index = chromosome_segment[0]
        input_string = chromosome_segment[1]

        for i in range(len(input_string) - 22):  # overlap_length
            substring = input_string[i:i+23]
                
            if substring[-2:] == 'GG':
                count = 0
                for char1, char2 in zip(substring[0:20], target_str[0:20]):
                    if char1 == char2:
                        count += 1
                          
                if count >= 15:      
                        """
                        The 1st element in the tuple is the position (by natural seuqence, starting from 1).
                        The 2nd element in the tuple is the sequence.
                        The 3rd element in the tuple is the number of match bases.
                        This does not happen frequently, therefore, we can put each result directly to the result_queue.
                        If this happens very often, the tuples should be collected before being put to the result_queue.
                        """
                        result_queue.put((i+start_index+1,substring,count)) 




"""
compare the off-target sequence and the target sequence
return the new off-target sequence, with mismatched bases tranformed to lower case
"""
def get_lowercase_off_target(input_string):
        output_string = ''
        for char1, char2 in zip(input_string[0:20], target_str[0:20]):
                if char1 == char2:
                        output_string = output_string + char1
                else:
                        output_string = output_string + char1.lower()
                        
        output_string = output_string + input_string[20:23]
        
        return output_string        










# calculate running time
start = time.time()


"""
set parameters
divide each chromosome to 1000 segments. This number ensure multicore searching.
If the chromosome is too small in length, the split_string function will report error,
tell user to decrease the segmentation number, i.e., set m = 100, or 50, etc.,
and rerun the program.
"""
m = 1000
overlap_length = 22


# list for the genome data
os.chdir(current_directory)
file_names_list = os.listdir('genome')

# get only the .fasta files to list
chromosome_file_names_list = []   
for file_name in file_names_list:
    file_name_split = file_name.split('.')
    if file_name_split[1] == 'fasta':
        chromosome_file_names_list.append(file_name)

chromosome_file_names_list_len = len(chromosome_file_names_list)



"""""""""""""""""""""""""""""""""""""""""""""""""""
"*** *** *** ***"     doing the first half job
"""""""""""""""""""""""""""""""""""""""""""""""""""


###search the sense strands
os.chdir(current_directory)
os.chdir('genome')
chromosome_off_target_result_list = []
sub_process_count = 8
# start subprocesses
set_start_method('fork')

for i in range(chromosome_file_names_list_len):
    input_string = str(SeqIO.read(chromosome_file_names_list[i], 'fasta').seq)

    print("Searching %s......" % chromosome_file_names_list[i])

    # divide the chromosome into m segments
    chromosome_segments = split_string(input_string, m, overlap_length)


    # multiprocessing
    result_queue = Queue()
    
    sub_process_list = []   
    for task_id in range(sub_process_count):
        p = Process(target=find_off_target_positions, args=(task_id, result_queue))
        sub_process_list.append(p)
        p.start()

    for p in sub_process_list:
        p.join()
    sub_process_list = []

    # get all the results   
    results = []
    while not result_queue.empty():  # get result from result_queue
        result = result_queue.get()
        results.append(result)

    # store the results in a new list
    tmp_list = []

    tmp_list.append(int(chromosome_file_names_list[i].split('.')[0]))    
    tmp_list.append(results)
    chromosome_off_target_result_list.append(tmp_list)
    
    tmp_list = []

  

chromosome_off_target_result_list_sorted = sorted(chromosome_off_target_result_list,key=itemgetter(0))
"""
print('\n',end='')
print(chromosome_off_target_result_list_sorted)
print('\n',end='')
"""


  
# get the list for output
output_list = []
for i in range(len(chromosome_off_target_result_list_sorted)):
        tmp_list = []

        tmp_list.append(chromosome_off_target_result_list_sorted[i][0])
        tmp_list.append(sorted(chromosome_off_target_result_list_sorted[i][1],key=itemgetter(0)))

        output_list.append(tmp_list)
        
        tmp_list = []
   
print('\n',end='')
# print(output_list)     

       
# store the result to txt files
# get current path
"""
print("\n", end = '')
print("Current directory is %s." % current_directory)
print("\n", end='')
"""
os.chdir(current_directory)

# make a directory for storing result txt files, if the folder has alreadly existed, first
# delete, then make the directory. 
foldername = 'off_target'
dirs = os.listdir(current_directory)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)

# go to the off_target folder
os.chdir(foldername)

sense_file = open("sense_strands.txt",'w')

# write results
print('\r\n',end='',file=sense_file)
print("chromosome".ljust(15,' '),\
      "start".ljust(15,' '),\
      "sequence".ljust(30,' '),\
      "end".ljust(15,' '),\
      "match".ljust(10,' '),\
      sep="",\
      file=sense_file
      )
print('\r\n',end='',file=sense_file)

print(30*' ',target_str.ljust(30,' '),sep="",file=sense_file)

print('\r\n',end='',file=sense_file)

for i in range(len(output_list)):
        for j in range(len(output_list[i][1])):
            C1_print_line_str = str(output_list[i][0])
            C2_print_line_str = str(output_list[i][1][j][0])
            C3_print_line_str = get_lowercase_off_target(output_list[i][1][j][1])
            C4_print_line_str = str(output_list[i][1][j][0]+22)
            C5_print_line_str = str(output_list[i][1][j][2])
            
            print(C1_print_line_str.ljust(15,' '),\
                  C2_print_line_str.ljust(15,' '),\
                  C3_print_line_str.ljust(30,' '),\
                  C4_print_line_str.ljust(15,' '),\
                  C5_print_line_str.ljust(10,' '),\
                  sep="",\
                  file=sense_file
                  )
        print('\r\n',end='',file=sense_file)
        
sense_file.close()


 
"""""""""""""""""""""""""""""""""""""""""""""""""""
"*** *** *** ***"     doing the other half job
"""""""""""""""""""""""""""""""""""""""""""""""""""


###search the antisense strands
os.chdir(current_directory)
os.chdir('genome')
chromosome_off_target_result_list = []
for i in range(chromosome_file_names_list_len):
    input_string = str(SeqIO.read(chromosome_file_names_list[i], 'fasta').seq.reverse_complement())

    print("Searching the reverse complement of %s......" % chromosome_file_names_list[i])

    # divide the chromosome into m segments
    chromosome_segments = split_string(input_string, m, overlap_length)


    # multiprocessing
    result_queue = Queue()
    
    sub_process_list = []   
    for task_id in range(sub_process_count):
        p = Process(target=find_off_target_positions, args=(task_id, result_queue))
        sub_process_list.append(p)
        p.start()

    for p in sub_process_list:
        p.join()
    sub_process_list = []

    # get all the results   
    results = []
    while not result_queue.empty():  # get result from result_queue
        result = result_queue.get()
        results.append(result)

    # store the results in a new list
    tmp_list = []

    tmp_list.append(int(chromosome_file_names_list[i].split('.')[0]))    
    tmp_list.append(results)
    tmp_list.append(len(input_string)) # record the chromosome length
    chromosome_off_target_result_list.append(tmp_list)
    
    tmp_list = []

  

chromosome_off_target_result_list_sorted = sorted(chromosome_off_target_result_list,key=itemgetter(0))
"""
print('\n',end='')
print(chromosome_off_target_result_list_sorted)
print('\n',end='')
"""


# get the list for output
output_list = []
for i in range(len(chromosome_off_target_result_list_sorted)):
        tmp_list = []

        tmp_list.append(chromosome_off_target_result_list_sorted[i][0])
        tmp_list.append(sorted(chromosome_off_target_result_list_sorted[i][1],key=itemgetter(0),reverse=True)) # in reverse order
        tmp_list.append(chromosome_off_target_result_list_sorted[i][2])

        output_list.append(tmp_list)
        
        tmp_list = []
   
print('\n',end='')
# print(output_list)     

       
# store the result to txt files
os.chdir(current_directory)
os.chdir(foldername)

antisense_file = open("antisense_strands.txt",'w')

# write results
print('\r\n',end='',file=antisense_file)
print("chromosome".ljust(15,' '),\
      "start".ljust(15,' '),\
      "sequence".ljust(30,' '),\
      "end".ljust(15,' '),\
      "match".ljust(10,' '),\
      sep="",\
      file=antisense_file
      )
print('\r\n',end='',file=antisense_file)

print(30*' ',target_str.ljust(30,' '),sep="",file=antisense_file)

print('\r\n',end='',file=antisense_file)

for i in range(len(output_list)):
        for j in range(len(output_list[i][1])):
            C1_print_line_str = str(output_list[i][0])
            C2_print_line_str = str(output_list[i][2]+1-output_list[i][1][j][0]) # position transformation
            C3_print_line_str = get_lowercase_off_target(output_list[i][1][j][1])
            C4_print_line_str = str(output_list[i][2]+1-output_list[i][1][j][0]-22)
            C5_print_line_str = str(output_list[i][1][j][2])
            
            print(C1_print_line_str.ljust(15,' '),\
                  C2_print_line_str.ljust(15,' '),\
                  C3_print_line_str.ljust(30,' '),\
                  C4_print_line_str.ljust(15,' '),\
                  C5_print_line_str.ljust(10,' '),\
                  sep="",\
                  file=antisense_file
                  )
        print('\r\n',end='',file=antisense_file)
        
antisense_file.close()


end = time.time()
print('Running time: %s Seconds'%(end-start))
print('\n',end='')





