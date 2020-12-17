# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 16:07:03 2020

@author: Valerya
"""

import pandas as pd
from math import e

class Alignment():
    def __init__(self, species, name, score, s1, s2, matrix, query_length, db_length, k=0.1, l=0.5):
        self.species = species
        self.name = name
        self.s1 = s1
        self.s2 = s2
        self.align = []
        self.score = score
        
        self.n = len(self.s1)
        self.identity = 0
        self.gap = 0
        self.query_cover = 0
        for i in range(self.n):
            if self.s1[i]==self.s2[i]:
                self.align.append('|')
                self.identity += 1
                self.query_cover += 1
            elif self.s1[i] == '-':
                self.gap += 1
                self.align.append(' ')
            elif self.s2[i] == '-':
                self.gap += 1
                self.query_cover += 1
                self.align.append(' ')
            elif matrix[self.s1[i]+self.s2[i]] > 0:
                self.align.append(':')
                self.query_cover += 1
            else:
                self.align.append('.')
                self.query_cover += 1
        
        self.query_cover = self.query_cover*100/query_length
        self.E = self.e_value(db_length, k, l)
    
    def visualize(self):
        m = (18 if self.n > 100 else 24)
        l = 0
        print('_'*(3*m+15)+'\n')
        print(self.species, self.name, "Score:%d"%self.score, "Query cover:%.d%%"%self.query_cover, "Identity:(%.d/%.d)%.2f"%(self.identity,self.n,self.identity/self.n), "Gap:(%.d/%.d)%.2f"%(self.gap,self.n,self.gap/self.n), "Expect:{:.2e}".format(self.E))
        print('_'*(3*m+15)+'\n')
        df = pd.DataFrame([self.s1,self.align,self.s2],index=['Query','','Sbjct'])
        pd.set_option('display.max_columns', None)
        while l < self.n:
            print(df.iloc[:,l:l+m],'\n')
            l += m
    
    def save(self, filepath):
        m = (18 if self.n > 100 else 24)
        with open(filepath,"at") as f :
            l = 0
            print('_'*(3*m+15)+'\n', file=f)
            print(self.species, self.name, "Score:%d"%self.score, "Query cover:%.d%%"%self.query_cover, "Identity:(%.d/%.d)%.2f"%(self.identity,self.n,self.identity/self.n), "Gap:(%.d/%.d)%.2f"%(self.gap,self.n,self.gap/self.n), "Expect:{:.2e}".format(self.E), file=f)
            print('_'*(3*m+15)+'\n', file=f)
            df = pd.DataFrame([self.s1,self.align,self.s2],index=['Query','','Sbjct'])
            pd.set_option('display.max_columns', None)
            while l < self.n:
                print(df.iloc[:,l:l+m],'\n', file=f)
                l += m
    
    def e_value(self, db_length, k, l):
        E = (k*self.n*db_length)*e**(-l*self.score)
        return E
    
    def verify_score(self,matrix,gap_penalty,extend_penalty):
        score = 0
        gap1 = False
        gap2 = False
        for i in range(self.n):
            if self.s1[i] == '-':
                if gap1 == False:
                    gap1 = True
                    gap2 = False
                    score += gap_penalty
                else:
                    score += extend_penalty
            elif self.s2[i] == '-':
                if gap2 == False:
                    gap2 = True
                    gap1 = False
                    score += gap_penalty
                else:
                    score += extend_penalty
            else:
                gap1 = False
                gap2 = False
                score += matrix[self.s1[i]+self.s2[i]]
        if score == self.score:
            print("The score is correct.")
        else:
            print(score)
        return score == self.score
        
if __name__ == "__main__":
    matrix = {}
    with open("matrix/BLOSUM62.txt") as f:
        letters = f.readline()[:-1]
        for line in f:
            line = line.split(',')
            matrix[line[0]] = int(line[1][:-1])
    al = Alignment("","",237,["","",""],["","",""],matrix,1E9, 1, 1)
    print(al.E)
    
    