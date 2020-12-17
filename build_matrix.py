# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 09:23:00 2020

@author: Valerya
"""
import re

filename = input("matrix name:")
file = "matrix/"+filename+"_ncbi.txt"

with open(file) as f:
    s = f.read()
    s = s.split("\n")
    s = list(map(lambda x:re.split('\\s+',x),s))
    letters = ''.join(s[0][:-1])+'U'
    

matrix = {}
for i in range(1,len(s)-1):
    for j in range(1,len(s)-1):
        matrix[s[0][i]+s[j][0]] = int(s[j][i])
for i in range(1,len(s)-1):
    matrix['U'+s[0][i]] = matrix[s[0][i]+'U'] = matrix[s[0][i]+'C']
matrix['UU'] = matrix['CC']
print(matrix)

with open("matrix/"+filename+".txt",'w') as f:
    f.write(letters+"\n")
    for key in matrix.keys():
        line = key+","+str(matrix[key])+"\n"
        f.write(line)
