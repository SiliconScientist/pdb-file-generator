# Title: pdb file generator
# Description: Turn ali file input into pdb output
# Author: Avery Hill
# Last modified: 08/09/2021

import os
from modeller import *
from modeller.automodel import *
import shutil

def remove_ipynb_checkpoints(path):
    flag = False
    
    for i in range(len(os.listdir(path))):
        if (os.listdir(path)[i] == '.ipynb_checkpoints'):
            flag = True
    if flag:
        os.rmdir(f'{path}/.ipynb_checkpoints')



directory_name = 'test_ali_files'
path = f'/Users/averyhill/MyDocuments/schoeffler_research_summer_2021/pdbs/pdb_generation_code/{directory_name}'
remove_ipynb_checkpoints(path)

file_array = os.listdir(path)
file_path = '../GyraBCC_library'
template_id = '1ab4'
num_files = len(file_array)
sequence_id = [None] * num_files
auto_model = [None] * num_files

# Run Dr. Schoeffler's collected scripts/files
# Run seq2struc_alignment.py code

env = Environ()
env.io.atom_files_directory = ['.', file_path]

len(file_array)

# File 0-100 have been generated. Next batch: 101-200
start = 0
end = len(file_array)

# Run Dr. Schoeffler Scripts on loop for each file in file_array
for i in range(start, end):
    fh = open(f'./{directory_name}/{file_array[i]}')
    # sequence_id example: The sequence_id of "to_align_AsubGaBCC.ali" is "AsubGaBCC"
    sequence_id[i] = file_array[i].split('_')[2].split('.')[0] + '_' + file_array[i].split('_')[3].split('.')[0]

    aln = Alignment(env, file = fh, align_codes = 'all')
    aln.salign(rr_file='$(LIB)/as1.sim.mat', 
               output='',
               max_gap_length=50,
               gap_function=True,          
               feature_weights=(1., 0., 0., 0., 0., 0.),
               gap_penalties_1d=(-100, 0),
               gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0.0, 0.0),
               similarity_flag=True)
    aln.write(file=f'{sequence_id[i]}_align2d.ali', alignment_format='PIR')

    env = Environ()

    env.io.atom_files_directory = ['.', file_path]

    # an array of automodels
    auto_model[i] = AutoModel(env, alnfile = f'{sequence_id[i]}_align2d.ali', knowns = template_id, sequence = sequence_id[i])
    auto_model[i].starting_model = 1
    auto_model[i].ending_model = 1
    auto_model[i].make()

# Remove extraneous files
dir_name = os.getcwd()
test = os.listdir(dir_name)

os.mkdir(f'{dir_name}/pdb_results')
os.mkdir(f'{dir_name}/align_results')

for item in test:
    
    # Move pdb files to folder
    if (item.endswith('.B99990001.pdb')):
        shutil.move("./" + item, "./pdb_results")
        
    # Move align2d files to folder
    if (item.endswith('.ali')):
        shutil.move('./' + item, './align_results')
        
    if (item.endswith(".V99990001") or item.endswith(".D00000001") or item.endswith(".rsr") or 
           item.endswith(".sch") or item.endswith(".ini") or item.endswith(".pap")):
            os.remove(os.path.join(dir_name, item))