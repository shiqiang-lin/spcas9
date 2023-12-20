# usage: python3.10 get_random_genome.py

"""
This script generates randome genome containing 10 chromosomes.
All the chromosome sequences are stored in fasta format.
We use the random.seed function to guarantee that the generated
chromosome sequences are identical everytime. If you need to
change, you can delete the random.seed from the generate_random_string function.
"""

import os
import random
import shutil
import time
from Bio.Seq import Seq
from Bio import SeqIO


def generate_random_string(length,seed_num):
    random.seed(seed_num)
    letters = "ATCG"
    return ''.join(random.choice(letters) for _ in range(length))



current_directory = os.getcwd()

# calculate running time
start = time.time()

os.chdir(current_directory)
# make a directory for storing result txt files, if the folder has alreadly existed, first
# delete, then make the directory. 
foldername = 'genome'
dirs = os.listdir(current_directory)

if(foldername not in dirs):
    os.mkdir(foldername)
else:
    shutil.rmtree(foldername)
    os.mkdir(foldername)

# go to the off_target folder
os.chdir(foldername)


print('\n',end='')
print("Generating random chromosomes.")
print('\n',end='')

# Very long chromosome, you can change to another number if you like.
string_length = 100000000 

for i in range(10):
    input_string = generate_random_string(string_length,100+i)

    record = SeqIO.SeqRecord(Seq(input_string), id='example_sequence_'+str(i+1), description='')
    SeqIO.write(record, str(i+1)+'.fasta', 'fasta')
    print(str(i+1)+'.fasta')

os.chdir(current_directory)

end = time.time()
print('\n',end='')
print('Running time: %s Seconds'%(end-start))
print('\n',end='')



    
