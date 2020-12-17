# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 17:38:33 2020

@author: Valerya
"""

import numpy as np
import time
from bisect import bisect_right
import re
from Alignment import Alignment
import h5py

class BLAST():
    def __init__(self, query, database_name, matrix_name, gap_penalty, extend_penalty, LIM, X, A, H, N, evalue_limit = 1e-10, word_len = 4):
        self.word_len = word_len
        self.query = query.upper()
        self.gap_penalty = gap_penalty
        self.extend_penalty = extend_penalty
        self.database_name = database_name
        
        self.matrix_name = matrix_name
        self.matrix = {}
        with open("matrix/"+matrix_name+".txt") as f:
            letters = f.readline()[:-1]
            for line in f:
                line = line.split(',')
                self.matrix[line[0]] = int(line[1][:-1])
        self.word_num = {}
        for i,l in enumerate(letters):
            self.word_num[l] = i # each AA matched to a number (0~24)
        
        self.indexed_db = h5py.File("database/"+database_name+".h5", "r")
        
        self.LIM = LIM # max displayed alignment
        self.X = X # during extending, lowest score accepted below best score
        self.A = A # max distance between two hits when they can be considered as one hit
        self.H = H # during denoising, min number of seeds on a hit start point that can be collected for furthur analysis
        self.N = N # during hit calling, min continuous seeds that can be considered a hit
        self.evalue_limit = evalue_limit
        
        # store the position of each sequence
        self.seq_pos = np.load("database/"+database_name+"_seq_id.npy")
        # store the gene name of each sequence
        self.seq_names = np.load("database/"+database_name+"_seq_name.npy")
        # store the start position of each organism
        self.species_idx = np.load("database/"+database_name+"_species_idx.npy")
        # store the name of each organism
        self.species_names = np.load("database/"+database_name+"_species_names.npy")

    
    def word_to_idx(self, word):
        word_index = 0
        for i,l in enumerate(word,1):
            word_index += (len(self.word_num))**(self.word_len-i)*self.word_num[l]
        return word_index
    
    def split_query(self):
        seeds = []
        query = self.query
        for i in range(len(query)-self.word_len+1):
            seed = query[:self.word_len]
            query = query[1:]
            seeds.append(seed)
        return seeds
    
    def extract_seq(self, pos, length):
        with open("database/"+self.database_name+"_seq.txt") as f:
            f.seek(pos,0)
            seq = f.read(length)
        return seq
    
    def get_seq_end(self,pos):
        idx = bisect_right(self.seq_pos, pos)
        seq_end = self.seq_pos[idx] if idx < len(self.seq_pos) else float('inf')
        return seq_end
    
    def get_seq_start(self,pos):
        idx = bisect_right(self.seq_pos, pos) - 1
        return self.seq_pos[idx]
    
    def get_seed_position(self,seeds):
        positions = []
        for seed in seeds:
            try:
                pos = self.indexed_db[seed][:]
            except KeyError:
                pos = []
            positions.append(pos)
        return positions

    def find_hit(self,positions):
        occurs = {}
        i = 0
        for position in positions:
            for pos in position:
                try:
                    occurs[pos-i].append(i)
                except KeyError:
                    occurs[pos-i] = [i]
            i += 1
            
        #print(occurs)
        occurs_count = []
        occurs_clear = {}
        for key,value in occurs.items():
            if len(value) > self.H:
                occurs_count.append((key,len(value)))
                occurs_clear[key] = value
        occurs.clear()
        occurs_count.sort()
        #print(occurs_count)
        #print(occurs_clear)
                    
        def next_gene(tmp,hit_lists,last_hit_finished,hit_start,seq_end):
            if hits[j-1] - hit_start > self.N:
                #print(hits[j-1] - hit_start+self.word_len)
                tmp.append((hit_start,hit_start+occurs_count[i][0],hits[j-1]-hit_start+self.word_len))
            tmp, last_hit_finished, hit_lists = add_hit(tmp, last_hit_finished, hit_lists)
            tmp = []
            last_hit_finished = True
            hit_start = hits[j]
            seq_end = self.get_seq_end(hit_start+occurs_count[i][0])
            
            return tmp,hit_lists,last_hit_finished,hit_start,seq_end
        
        def add_hit(tmp, last_hit_finished, hit_lists):
            if len(tmp) == 0:
                return tmp, last_hit_finished, hit_lists
            
            #print(last_hit_finished)
            if last_hit_finished:
                hit_lists.append(tmp)
                last_hit_finished = False
            else:
                last_seq_end = hit_lists[-1][-1][2] + hit_lists[-1][-1][1]
                gap = tmp[0][1] - last_seq_end
                #print(gap)
                if abs(gap-self.word_len) < self.A and seq_end == self.get_seq_end(last_seq_end):
                    while len(tmp) > 0:
                        last_hit_end = hit_lists[-1][-1][0] + hit_lists[-1][-1][2]
                        if last_hit_end >= tmp[0][0]:
                            if last_hit_end >= tmp[0][0] + tmp[0][2]:
                                #print("drop", tmp, hit_lists[-1][-1])
                                tmp = tmp[1:]
                            else:
                                overlap = last_hit_end-tmp[0][0]
                                #print("overlap",overlap,hit_lists[-1][-1],tmp[0])
                                tmp[0] = (tmp[0][0]+overlap, tmp[0][1]+overlap, tmp[0][2]-overlap)
                                hit_lists[-1] += tmp
                                break
                        else:
                            hit_lists[-1] += tmp
                            #print("overlap", tmp, hit_lists[-1])
                            break
                            
                else:
                    hit_lists.append(tmp)
            return tmp, last_hit_finished, hit_lists
            
        last_hit_finished = True
        hit_lists = []
        for i in range(len(occurs_count)):
            hits = occurs_clear[occurs_count[i][0]]
            hit_start = hits[0]
            seq_end = self.get_seq_end(hit_start+occurs_count[i][0])
            tmp = []
            
            for j in range(1,len(hits)):
                if hits[j]-hits[j-1] > 1:
                    if hits[j] + occurs_count[i][0] + self.word_len > seq_end:
                        tmp,hit_lists,last_hit_finished,hit_start,seq_end = next_gene(tmp,hit_lists,last_hit_finished,hit_start,seq_end)
                    elif hits[j] < hits[j-1]+self.A:
                        if hits[j-1] - hit_start > self.N:
                            #print(hits[j-1] - hit_start+self.word_len)
                            tmp.append((hit_start,hit_start+occurs_count[i][0],hits[j-1]-hit_start+self.word_len))
                    else:
                        tmp,hit_lists,last_hit_finished,hit_start,seq_end = next_gene(tmp,hit_lists,last_hit_finished,hit_start,seq_end)
                    hit_start = hits[j]
            if hits[j] < hits[j-1]+self.A and hits[j-1] - hit_start > self.N:
                #print(hits[j-1] - hit_start)
                tmp.append((hit_start, hit_start+occurs_count[i][0],hits[j-1]-hit_start+self.word_len))
            
            tmp, last_hit_finished, hit_lists = add_hit(tmp, last_hit_finished, hit_lists)
            
        if len(hit_lists)==0:
            return -1
        else:
            return hit_lists
    
    # needleman-wunsch algorithm
    def needle(self, Q, S):
    
        m = len(Q)
        n = len(S)
    
        if m == 0:
            score = self.gap_penalty + (n-1)*self.extend_penalty
            return score, ['-' for _ in range(n)], list(S)
        elif n == 0:
            score = self.gap_penalty + (m-1)*self.extend_penalty
            return score, list(Q), ['-' for _ in range(m)]
    
        match = self.matrix[Q[0]+S[0]]
        
        #DP scoring matrix
        x = np.zeros((m,n)) #each value in the matrix computed from the left value, last step = gap(s) in query
        y = np.zeros((m,n)) #each value computed from the upper value, last step = gap(s) in reference sequence
        mt = np.zeros((m,n)) #each value computed from the upper left value, last step = match or substitution
        # vertical:length of s1(query)  horizontal:length of s2(Refseq)
        
        #complete start point
        x[0,0] = 2*self.gap_penalty #horizontal gap + vertical gap
        y[0,0] = 2*self.gap_penalty
        mt[0,0] = match
        
        #complete first column
        for i in range(1,m):
            match = self.matrix[Q[i]+S[0]]
            pointer = (i,0)
            p_pointer = (i-1,0)
            open_gap = mt[p_pointer] + self.gap_penalty
            extend_gap = y[p_pointer] + self.extend_penalty
        
            if open_gap >= extend_gap:
                y[pointer] = open_gap
            else:
                y[pointer] = extend_gap
        
            mt[pointer] = match + (self.gap_penalty + (i - 1)*self.extend_penalty)
            x[pointer] = i*self.extend_penalty + 2*self.gap_penalty
        
        #complete first row
        for j in range(1,n):
            match = self.matrix[Q[0]+S[j]]
            pointer = (0,j)
            p_pointer = (0,j-1)
         
            open_gap = mt[p_pointer] + self.gap_penalty
            extend_gap = x[p_pointer] + self.extend_penalty
            
            if open_gap >= extend_gap:
                x[pointer] = open_gap
            else:
                x[pointer] = extend_gap

            mt[pointer] = match + (self.gap_penalty + (j - 1)*self.extend_penalty)
            y[pointer] = j*self.extend_penalty + 2*self.gap_penalty
        
        #loop through the matrix
        j = 1
        while j < n:
            i = 1
            S_base = S[j]
        
            p_pointer = (0,j-1)
            pointer = (0,j)
            j += 1
            while i < m:
                match = self.matrix[Q[i]+S_base]
                i += 1
                
                pointer = (pointer[0]+1,pointer[1])
                
                # for matrix mt, p_pointer = upper left
                mt[pointer] = max(mt[p_pointer],x[p_pointer],y[p_pointer]) + match
                
                # for matrix y, p_pointer = upper
                p_pointer = (p_pointer[0], p_pointer[1]+1)
                
                open_gap = max(x[p_pointer],mt[p_pointer]) + self.gap_penalty
                extend_gap = y[p_pointer] + self.extend_penalty
            
                if open_gap > extend_gap:
                    y[pointer] = open_gap
                else:
                    y[pointer] = extend_gap
                
                # for matrix x, p_pointer = left
                p_pointer = (p_pointer[0]+1, p_pointer[1]-1)
                
                open_gap = max(mt[p_pointer],y[p_pointer]) + self.gap_penalty
                extend_gap = x[p_pointer] + self.extend_penalty
            
                if open_gap > extend_gap:
                    x[pointer] = open_gap
                else:
                    x[pointer] = extend_gap
            
        #print(mt)
        #print(x)
        #print(y)
        #print(Q,S)
        return self.traceback(Q, S, mt, x, y, m, n)
    
    def aboveX(self, n, max_score):
        if n < max_score+self.X:
            n = -99
        return n
    
    def fill_matrix(self, t, query, refseq, mt, x, y, ypointer, xpointer, max_score, max_pos):
        # upper left
        mt[ypointer].append(self.aboveX(max(mt[ypointer-1][xpointer-1], x[ypointer-1][xpointer-1], y[ypointer-1][xpointer-1]) + self.matrix[query[ypointer]+refseq[xpointer]], max_score))
        # left
        x[ypointer].append(self.aboveX(max(max(mt[ypointer][xpointer-1],y[ypointer][xpointer-1]) + self.gap_penalty, x[ypointer][xpointer-1] + self.extend_penalty), max_score))
        # upper
        y[ypointer].append(self.aboveX(max(max(mt[ypointer-1][xpointer],x[ypointer-1][xpointer]) + self.gap_penalty, y[ypointer-1][xpointer] + self.extend_penalty), max_score))
        current = max(mt[ypointer][xpointer], x[ypointer][xpointer], y[ypointer][xpointer])
        #print(current)
        if current > self.X:
            t = True
            if current > max_score:
                # update max score and position
                max_score = current
                max_pos = (ypointer,xpointer)
        return t, mt, x, y, max_score, max_pos
    
    def extend_tail(self,tail,pos,X=-15):
        
        l = len(tail)
        if l == 0:
            return 0, [], []
    
        idx = bisect_right(self.seq_pos, pos-1)
        end = float('inf')
        if idx < len(self.seq_pos):
            end = self.seq_pos[idx]
            
        segment = 20
        if end-pos < segment:
            segment = end-pos
        seq = self.extract_seq(pos,segment)
        L = len(seq)
        pos += L
        #print(seq,L,pos)
        if L == 0:
            return 0, [], []
    
        match = self.matrix[tail[0]+seq[0]]
        if match >= 0:
            max_score = match
            max_pos = (0,0)
        else:
            max_score = 0
            max_pos = None
        
        y = [[2*self.gap_penalty]]
        x = [[2*self.gap_penalty]]
        mt = [[self.aboveX(match, max_score)]]
    
        pointer = 1
        t = True
    
        while t:
            t = False
            if pointer == L:
                if end-pos < segment:
                    segment = end-pos
                add_seq = self.extract_seq(pos,segment)
                if len(add_seq)==0 and pointer >= l:
                    break
                L += len(add_seq)
                pos += len(add_seq)
                seq += add_seq
                #print(seq,L,pos)
                
            # Complete first point in vertical edge
            if pointer < L:
                mt[0].append(self.aboveX(self.matrix[tail[0]+seq[pointer]]+self.gap_penalty+(pointer-1)*self.extend_penalty, max_score))
                x[0].append(self.aboveX(max(x[0][pointer-1]+self.extend_penalty, mt[0][pointer-1]+self.gap_penalty), max_score))
                y[0].append(self.aboveX(y[0][pointer-1]+self.extend_penalty, max_score))
                current = max(mt[0][pointer], x[0][pointer], y[0][pointer])
                #print(current)
                if current > X:
                    t = True
                    if current > max_score:
                        max_score = current
                        max_pos = (0,pointer)
            # Complete first point in horizontal edge
            if pointer < l:
                mt.append([self.aboveX(self.matrix[tail[pointer]+seq[0]]+self.gap_penalty+(pointer-1)*self.extend_penalty, max_score)])
                x.append([self.aboveX(x[pointer-1][0]+self.extend_penalty, max_score)])
                y.append([self.aboveX(max(y[pointer-1][0]+self.extend_penalty, mt[pointer-1][0]+self.gap_penalty), max_score)])
                current = max(mt[pointer][0], x[pointer][0], y[pointer][0])
                if current > X:
                    t = True
                    if current > max_score:
                        max_score = current
                        max_pos = (pointer,0)
            
            for i in range(1,pointer):
                #complete vertical edge
                if i < l and pointer < L:
                    t, mt, x, y, max_score, max_pos = self.fill_matrix(t, tail, seq, mt, x, y, i, pointer, max_score, max_pos)
                #complete horizontal edge
                if pointer < l and i < L:
                    t, mt, x, y, max_score, max_pos = self.fill_matrix(t, tail, seq, mt, x, y, pointer, i, max_score, max_pos)
            
            #Complete intersection point of edges
            if pointer < l and pointer < L:
                t, mt, x, y, max_score, max_pos = self.fill_matrix(t, tail, seq, mt, x, y, pointer, pointer, max_score, max_pos)
                
            pointer += 1
    
        #print(pointer,l,L)
        #print(max_pos)
        #print(max_score)
        
        if max_pos == None:
            return max_score, [], []
        
        i = max_pos[0] + 1
        j = max_pos[1] + 1
    
        mt = np.array(mt)
        x = np.array(x)
        y = np.array(y)
        
        #print("tail", seq, tail, i, j)
        return self.traceback(tail, seq, mt, x, y, i, j)
    
    def extend_head(self, head, pos): # same as extend_tail, but sequences are reversed
        
        l = len(head)
        if l == 0:
            return 0, [], []
            
        head = head[::-1] #reversed
    
        idx = bisect_right(self.seq_pos, pos) - 1
        end = self.seq_pos[idx]
    
        segment = 20 # length of the sequence segment retrieved fromdatabase each time
        if pos-end < segment: # near the end of the gene (distance < segment length)
            segment = pos-end
        pos -= segment
        seq = self.extract_seq(pos,segment)[::-1] #reversed
        L = len(seq)
        #print(seq,L,pos)
        if L == 0:
            return 0, [], []
        
        match = self.matrix[head[0]+seq[0]]
        if match >= 0:
            max_score = match
            max_pos = (0,0)
        else:
            max_score = 0
            max_pos = None
        
        y = [[2*self.gap_penalty]]
        x = [[2*self.gap_penalty]]
        mt = [[self.aboveX(match, max_score)]] #size is unfixed
        
        pointer = 1
        t = True
    
        while t:
            t = False
            # When refseq comes to an end, retrive more sequence from database
            if pointer == L:
                if pos-end < segment: # near the end of the gene (distance < segment length)
                    segment = pos-end
                if segment==0 and pointer >= l: # at the end of the gene sequence & query sequence, end loop
                    break
                pos -= segment 
                add_seq = self.extract_seq(pos,segment)[::-1]
                L += segment 
                seq += add_seq
                #print(seq,L,pos)
            if pointer < L:
                mt[0].append(self.aboveX(self.matrix[head[0]+seq[pointer]]+self.gap_penalty+(pointer-1)*self.extend_penalty, max_score))
                x[0].append(self.aboveX(max(x[0][pointer-1]+self.extend_penalty, mt[0][pointer-1]+self.gap_penalty), max_score))
                y[0].append(self.aboveX(y[0][pointer-1]+self.extend_penalty, max_score))
                current = max(mt[0][pointer], x[0][pointer], y[0][pointer])
                #print(current)
                if current > self.X:
                    t = True
                    if current > max_score:
                        max_score = current
                        max_pos = (0,pointer)
            if pointer < l:
                mt.append([self.aboveX(self.matrix[head[pointer]+seq[0]]+self.gap_penalty+(pointer-1)*self.extend_penalty, max_score)])
                x.append([self.aboveX(x[pointer-1][0]+self.extend_penalty, max_score)])
                y.append([self.aboveX(max(y[pointer-1][0]+self.extend_penalty, mt[pointer-1][0]+self.gap_penalty), max_score)])
                current = max(mt[pointer][0], x[pointer][0], y[pointer][0])
                if current > self.X:
                    t = True
                    if current > max_score:
                        max_score = current
                        max_pos = (pointer,0)
            for i in range(1,pointer):
                if i < l and pointer < L:
                    t, mt, x, y, max_score, max_pos = self.fill_matrix(t, head, seq, mt, x, y, i, pointer, max_score, max_pos)
                if pointer < l and i < L:
                    t, mt, x, y, max_score, max_pos = self.fill_matrix(t, head, seq, mt, x, y, pointer, i, max_score, max_pos)
            if pointer < l and pointer < L:
                t, mt, x, y, max_score, max_pos = self.fill_matrix(t, head, seq, mt, x, y, pointer, pointer, max_score, max_pos)
            pointer += 1
    
        #print(pointer,l,L)
        #print(max_pos)
        #print(max_score)
        
        if max_pos == None:
            return max_score, [], []
        
        i = max_pos[0] + 1
        j = max_pos[1] + 1
    
        mt = np.array(mt)
        x = np.array(x)
        y = np.array(y)
        #print("head", seq, head, i, j)
        return self.traceback(head, seq, mt, x, y, i, j, True)
    
    def traceback(self, s1, s2, mt, x, y, m, n, reverse = False):
    
        p_pointer = 0 # direction of previous step

        i = m - 1 
        j = n - 1
    
        trace = np.empty((m,n), dtype=int) # vertical:length of s1(query)  horizontal:length of s2(Refseq)
        
        while i>=0 and j>=0:
            pointer = (i,j)
            mp = mt[pointer]
            # previous step is 1 and the right value - this value in x matrix is extend penalty
            if p_pointer == 1 and self.extend_penalty == x[pointer[0],pointer[1]+1]-x[pointer]:
                trace[pointer] = 1 #Left, gap(s) in query sequence
                j -= 1
            # previous step is 2 and the lower value - this value in y matrix is extend penalty
            elif p_pointer == 2 and self.extend_penalty == y[pointer[0]+1,pointer[1]]-y[pointer]:
                trace[pointer] = 2 #Up, gap(s) in reference sequence
                i -= 1
            elif mp >= x[pointer] and mp >= y[pointer]: #mp is max
                '''
                if p_pointer == 1 and mp == x[pointer]:
                    print("1")
                    trace[pointer] = 1
                    j -= 1
                elif p_pointer == 2 and mp == y[pointer]:
                    print("2")
                    trace[pointer] = 2
                    i -= 1
                else:
                    print("0")
                '''
                trace[pointer] = 0; # Diagonal, matched
                i -= 1
                j -= 1

            elif x[pointer]>=y[pointer] and j > -1: #x[pointer] is max
                trace[pointer] = 1
                j -= 1
            elif i > -1: #y[pointer] is max
                trace[pointer] = 2
                i -= 1
            else:
                print("Warning: Error encountered in Needleman-Wunsch traceback algorithm")
                return -1
            p_pointer = trace[pointer]
    
        S1 = []
        S2 = []
    
        i = m - 1
        j = n - 1
        
        while j >= 0 and i >= 0:
            #print(i,j)
            pointer = (i,j)
            if trace[pointer] == 0: #Diagonal, matched
                S1.append(s1[i])
                S2.append(s2[j])
                #print(align_s,align_q)
                i -= 1
                j -= 1
                continue
            elif trace[pointer] == 1: #Left, gap(s) in query
                S1.append('-')
                S2.append(s2[j])
                #print(align_s,align_q)
                j -= 1
                continue
            elif trace[pointer] == 2: #Up, gap(s) in reference sequence
                S1.append(s1[i])
                S2.append('-')
                #print(align_s,align_q)
                i -= 1
                continue
            else:
                print("Warning: Error encountered in Needleman-Wunsch traceback algorithm");
                return -1
        
        # complete gaps in the end of the sequence
        while j >= 0:
            S2.append(s2[j])
            S1.append('-')
            j -= 1
        while i >= 0:
            S1.append(s1[i])
            S2.append('-')
            i -= 1
        if not reverse:   
            S1.reverse()
            S2.reverse()
    
        pointer = (m-1,n-1);
        if mt[pointer] > x[pointer] and mt[pointer] > y[pointer]:
            score = mt[pointer]
        elif x[pointer] > y[pointer]:
            score = x[pointer]
        else:
            score = y[pointer]
        
        #print("".join(S1))
        #print("".join(S2))
        return score, S2, S1
    
    def extend_hit(self,hits):
        alignments = []
        idx = 0
        db_length = self.seq_pos[-1]
        #print(db_length)
        for hit in hits:
            name = self.seq_names[bisect_right(self.seq_pos,hit[0][1])-1]
            species = self.species_names[bisect_right(self.species_idx,hit[0][1])-1]
            score, align_s, align_q = self.extend_head(self.query[:hit[0][0]], hit[0][1])
            #print(idx, score, align_s, align_q)
            for i in range(len(hit)-1):
                tmp_s = list(self.extract_seq(hit[i][1], hit[i][2]))
                tmp_q = list(self.query[hit[i][0]:hit[i][0]+hit[i][2]])
                tmp_score = 0
                for j in range(hit[i][2]):
                    tmp_score += self.matrix[tmp_s[j]+tmp_q[j]]
                #print(idx, tmp_score, tmp_s, tmp_q)
                align_s += tmp_s
                align_q += tmp_q
                score += tmp_score
                
                tmp_s = self.extract_seq(hit[i][1]+hit[i][2], hit[i+1][1]-hit[i][1]-hit[i][2])
                tmp_q = self.query[hit[i][0]+hit[i][2]:hit[i+1][0]]
                tmp_score, tmp_q, tmp_s = self.needle(tmp_s,tmp_q)
                #print(idx, tmp_score, tmp_s, tmp_q)
                align_s += tmp_s
                align_q += tmp_q
                score += tmp_score
            tmp_s = list(self.extract_seq(hit[-1][1], hit[-1][2]))
            tmp_q = list(self.query[hit[-1][0]:hit[-1][0]+hit[-1][2]])
            tmp_score = 0
            for j in range(hit[-1][2]):
                tmp_score += self.matrix[tmp_s[j]+tmp_q[j]]
            align_s += tmp_s
            align_q += tmp_q
            score += tmp_score
            tmp_score, tmp_s, tmp_q = self.extend_tail(self.query[hit[-1][0]+hit[-1][2]:], hit[-1][1]+hit[-1][2])
            align_s += tmp_s
            align_q += tmp_q
            score += tmp_score
            
            idx += 1
            alignment = Alignment(species, name, score, align_q, align_s, self.matrix, len(self.query), db_length)
            if alignment.E < self.evalue_limit:
                alignments.append(alignment)
                
        alignments.sort(key=lambda alignments: alignments.score, reverse=True)
        alignments = alignments[:self.LIM]
        return alignments
    
    def blast(self,filename):
        filename = "result/"+filename+".txt"
        seeds = self.split_query()
        positions = self.get_seed_position(seeds)
        #print(positions)
        hits = BLAST.find_hit(self,positions)
        if hits == -1:
            print("Not found")
            return -1
        alignments = self.extend_hit(hits)
        with open(filename, 'w') as f:
            f.write("BLAST output ")
            f.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+"\n")
            f.write("Query = "+self.query+"\n")
            f.write("Similarity score matrix: "+self.matrix_name+"\n")
            f.write("Use settings: gap penalty="+str(self.gap_penalty)+" gap extend penalty="+str(self.extend_penalty)+" X="+str(self.X)+" A="+str(self.A)+" H="+str(self.H)+" N="+str(self.N)+"\n")
            f.write("Display max "+str(self.LIM)+" matches.\n")
        for alignment in alignments:
            alignment.visualize()
            alignment.save(filename)
            ###alignment.verify_score(self.matrix,self.gap_penalty,self.extend_penalty)
        ###print(hits)
        self.indexed_db.close()
        print("The saved result filepath is "+filename)
        
        
if __name__ == '__main__':
    query = input("Query: ")
    
    filename = input("Output result file name: ")
    if filename == "":
        filename = time.strftime("%Y.%m.%d_%H.%M.%S", time.localtime())+"_blast"
    
    settings = input("Othersettings: ")
    settings = settings.split("-")[1:]
    
    database = "database/combined"
    matrix = "BLOSUM62"
    
    gap_penalty = -10
    extend_penalty = -0.5
    LIM = 20
    X = -15
    A=70
    H=4
    N=2
    
    try:
        for setting in settings:
            if re.match("^[Gg] ", setting):
                gap_penalty = -float(re.sub(" ","",setting[1:]))
            elif re.match("^[Ee] ", setting):
                extend_penalty = -float(re.sub(" ","",setting[1:]))
            elif re.match("^[Ll] ", setting):
                LIM = int(re.sub(" ","",setting[1:]))
            elif re.match("^[Xx] ", setting):
                X = -float(re.sub(" ","",setting[1:]))
            elif re.match("^[Aa] ", setting):
                A = int(re.sub(" ","",setting[1:]))
            elif re.match("^[Hh] ", setting):
                H = int(re.sub(" ","",setting[1:]))
            elif re.match("^[Nn] ", setting):
                N = int(re.sub(" ","",setting[1:]))
            elif re.match("^[Mm] ", setting):
                matrix = re.sub(" ","",setting[1:]).upper()
            elif re.match("^[Dd] ", setting):
                database = re.sub(" ","",setting[1:]).upper()
            else:
                print("invalid setting")
                break
    except ValueError:
        print("invalid setting")
    
    ###print(gap_penalty, extend_penalty, LIM, X, A, H, N)
    
    blastwave = BLAST(query,database,matrix,gap_penalty,extend_penalty,LIM,X,A,H,N)
    
    blastwave.blast(filename)

   
   