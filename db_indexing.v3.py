# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 08:54:08 2020

@author: Valerya
"""

import re
import numpy as np
import h5py

def word_to_num(letters):
    word_num = {}
    for i,l in enumerate(letters):
        word_num[l] = i
    return word_num

def word_to_idx(word_num,word,word_len):
    word_index = 0
    for i,l in enumerate(word,1):
        word_index += (len(word_num))**(word_len-i)*word_num[l]
    return word_index
    
def build_index(seq_name,seq,idx,indexed_db,seq_idx,seq_names,word_num):
    print(seq_name)
    seq_idx.append(idx)
    seq_names.append(seq_name)
    for i in range(len(seq)-word_len+1):
        word = seq[:word_len]
        seq = seq[1:] #the first AA will not be read again, so we discard it, and 'seq' is gradually shortening
        word_index = word_to_idx(word_num,word,word_len)
        #A 4 AA sequence can be regarded as a 25 band digit. Next we will turn 25 band digit to 10 band index.
        indexed_db[word_index].append(idx)
        #print(idx)
        idx += 1
    return idx
        
def build_lib(db_path_list,output_path,letters,word_len):
    
    indexed_db = [[] for _ in range(len(letters)**word_len)]

    word_num = word_to_num(letters)
    
    seq_idx = []
    seq_names = []
    species_names = []
    species_idx = []
    idx = 0
    
    with open("database/"+output_path+"_seq.txt",'w') as f:
        for db_path in db_path_list:
            species_names.append(db_path.split("/")[0])
            species_idx.append(idx)
            print(db_path.split("/")[0])
            with open("database/"+db_path+".fasta") as db:
                line = db.readline()
                seq_name = line.split(" ")[-3]
                if re.match("^GN=",seq_name):
                    seq_name = seq_name[3:]
                else:
                    seq_name = "unnamed"
                seq = ""
                for line in db:
                    if re.match(r">",line):
                        f.write(seq)
                        idx = build_index(seq_name,seq,idx,indexed_db,seq_idx,seq_names,word_num) + 3
                        seq = ""
                        seq_name = line.split(" ")[-3]
                        if re.match("^GN=",seq_name):
                            seq_name = seq_name[3:]
                        else:
                            seq_name = "unnamed"
                    else:
                        seq += line[:-1]
                idx = build_index(seq_name,seq,idx,indexed_db,seq_idx,seq_names,word_num) + 3
                f.write(seq)
    
    print("\nSaving library...")
    seq_idx.append(idx)
        
    p = 0
    with h5py.File("database/"+output_path+".h5", "w") as f:
        for i in letters:
            for j in letters:
                for k in letters:
                    for l in letters:
                        if len(indexed_db[p]) > 0:
                            f.create_dataset(i+j+k+l, data = indexed_db[p])
                            p += 1
                        else:
                            indexed_db.pop(p)
                        

    #indexed_db = np.array(indexed_db)
    seq_idx = np.array(seq_idx)
    seq_names = np.array(seq_names)
    species_names = np.array(species_names)
    species_idx = np.array(species_idx)
    #np.save("database/"+db_path+".npy",indexed_db)
    np.save("database/"+output_path+"_seq_id.npy",seq_idx)
    np.save("database/"+output_path+"_seq_name.npy",seq_names)
    np.save("database/"+output_path+"_species_idx.npy",species_idx)
    np.save("database/"+output_path+"_species_names.npy",species_names)
    print("genes:",len(seq_idx))
    print("k-mers:",len(indexed_db))
    
    return indexed_db,seq_idx
    
if __name__ == "__main__":
    
    matrix_name = 'BLOSUM62'
    word_len = 4
    
    input_file = input("Input file: ")
    if input_file == "":
        input_file = "database/database_list.txt"
        
    output_path = input("Output path: ")
    if output_path == "":
        output_path = "database/combined"
    #Homo_sapiens/uniprot-proteome_UP000005640
    #C.elegans/uniprot-proteome_UP000001940
    #Mus_musculus/uniprot-proteome_UP000000589
    #Zebrafish/uniprot-proteome_UP000000437
    #Gorilla/uniprot-proteome_UP000001519
    #Drosophila/uniprot-proteome_UP000000803
    
    with open("matrix/"+matrix_name+".txt") as f:
        letters = f.readline()[:-1]
    
    with open(input_file, 'r') as f:
        db_path_list = f.read()
        db_path_list = db_path_list.split('\n')
    
    
    indexed_db,seq_idx = build_lib(db_path_list,output_path,letters,word_len)
    
