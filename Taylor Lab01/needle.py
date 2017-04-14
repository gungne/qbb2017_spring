#!/usr/bin/env python

"""
Perform Needleman-Wunsch global alignment of two nucleotide sequences.
usage: needle.py sequence_1 sequence_2
"""

import sys
import numpy as np
global sigma
# HoxD70 matrix of Chiaromonte, Yap, Miller 2002,
#              A     C     G     T
sigma = [ [   91, -114,  -31, -123 ],
          [ -114,  100, -125,  -31 ],
          [  -31, -125,  100, -114 ],
          [ -123,  -31, -114,   91 ] ]

gap = 300
def trace_back(tb_value,s_letter,s_value,t_letter,t_value):
    if tb_value == 0:
        s_letter = s_letter
        t_letter = t_letter
        s_value = s_value-1
        t_value = t_value-1
    elif tb_value == 1:
        s_letter = '-'
        t_letter = t_letter
        s_value = s_value
        t_value = t_value-1
    elif tb_value == 2:
        s_letter = s_letter
        t_letter = '-'
        s_value = s_value-1
        t_value = t_value
    else:
        print('warning')
            
        
    return s_letter,s_value,t_letter,t_value

def nt_value(nt):
    if nt == 'A':
        value = 0
    elif nt == 'C':
        value = 1
    elif nt == 'G':
        value = 2
    elif nt == 'T':
        value = 3
    else:
        print('warning nt error')
    return int(value)

def fill_matrix(seg_matrix,s_nt,t_nt):
    #segment the matrix into 2x2
    global sigma
    diag_dp = seg_matrix[0,0] + sigma[nt_value(s_nt)][nt_value(t_nt)]
    hori_dp = seg_matrix[1,0] - gap
    vent_dp = seg_matrix[0,1] - gap

    dp_value = max(diag_dp,hori_dp,vent_dp)
    if dp_value == diag_dp:
        tb_value = 0
    elif dp_value == vent_dp :
        tb_value = 1
    elif dp_value== hori_dp:
        tb_value = 2
    else:
        print('error')

    return dp_value,tb_value

def compute_matrix( s, t ):
    """
    Fill in similarity score matrix, keeping traceback.
    """
    n = len( s )
    m = len( t )
    dp_matrix = np.zeros( (m+1,n+1), float )
    tb_matrix = np.zeros( (m+1,n+1), int )

    # ...
    for n_count in range(1,n+1):
        dp_matrix[0,n_count] = n_count * -1*gap      
        tb_matrix[0,n_count] = 2
    for m_count in range(1,m+1):
        dp_matrix[m_count,0] = m_count * -1*gap
        tb_matrix[m_count,0] = 1   
    for n_count in range(1,n+1):
        for m_count in range(1,m+1):
            # print(dp_matrix[np.ix_([m_count-1,m_count],[n_count-1,n_count])])
            dp_matrix[m_count,n_count],tb_matrix[m_count,n_count] = fill_matrix(dp_matrix[np.ix_([m_count-1,m_count],[n_count-1,n_count])],s[n_count-1],t[m_count-1]);



    return dp_matrix, tb_matrix

def print_alignment( dp_matrix, tb_matrix, s, t ):
    s_align = ''
    t_align = ''
    t_value = len(t)
    s_value = len(s)
    while s_value>0 or t_value>0:
        tb_value = tb_matrix[t_value,s_value]
        t_letter = t[t_value-1]
        s_letter = s[s_value-1]
        s_letter,s_value,t_letter,t_value = trace_back(tb_value,s_letter,s_value,t_letter,t_value)
        s_align+=s_letter
        t_align+=t_letter
        

    # print(start_sitem,start_siten)


    print(s_align[::-1])
    print(t_align[::-1])
    print("Score: ", dp_matrix[-1,-1])

def main():
    s = sys.argv[1]
    t = sys.argv[2]
    dp_matrix, tb_matrix = compute_matrix( s, t )
    print(dp_matrix)
    print(tb_matrix)
    print_alignment( dp_matrix, tb_matrix, s, t )

if __name__ == "__main__":
    main()