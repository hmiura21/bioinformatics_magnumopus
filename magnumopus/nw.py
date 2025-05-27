#!/usr/bin/env python3


def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    #make matrix:
    v = len(seq_a)
    h = len(seq_b)
    score_matrix=[]
    direction_matrix=[]
    for i in range(v+1):
        score_matrix.append((h+1) * [0])
        direction_matrix.append((h+1) * [0])
    

    #fill first row and column with gap scores
    for i in range(1,v+1): 
        score_matrix[i][0]=score_matrix[i-1][0]+gap
    for j in range(1,h+1):
        score_matrix[0][j]=score_matrix[0][j-1]+gap

    #calculate diagnol match/mismatch scores
    for i in range(1,v+1): 
        for j in range(1,h+1): 
            if seq_b[j-1]==seq_a[i-1]:
                match_score=match
            else:
                match_score=mismatch
            #calculate scores from all direction and choose highest score 
            score_fromD=score_matrix[i-1][j-1]+match_score
            score_fromL=score_matrix[i][j-1]+gap
            score_fromU=score_matrix[i-1][j]+gap
            chosen_score=max(score_fromD,score_fromL,score_fromU)
            score_matrix[i][j]=chosen_score
            if chosen_score==score_fromD:
                direction_matrix[i][j]="D"
            elif chosen_score==score_fromL:
                direction_matrix[i][j]="L"
            elif chosen_score==score_fromU:
                direction_matrix[i][j]="U"



    #initialize direction matrix starting with last row last column and make empty strings for final seqs
    updated_seq_a=""
    updated_seq_b=""
    i=len(direction_matrix)-1
    j=len(direction_matrix[0])-1

    score=score_matrix[i][j]

    #add seq to string while moving back to top left 
    while (i!=0 and j!=0): 
        location=direction_matrix[i][j] 
        if location=="D":
            updated_seq_b+=seq_b[j-1]
            updated_seq_a+=seq_a[i-1]
            i-=1
            j-=1
        elif location=="L":
            updated_seq_a+="-"
            updated_seq_b+=seq_b[j-1]
            j-=1
        elif location=="U":
            updated_seq_b+="-"
            updated_seq_a+=seq_a[i-1]
            i-=1

    reversed_seq_a = updated_seq_a[::-1]
    reversed_seq_b = updated_seq_b[::-1]

    final=((reversed_seq_a,reversed_seq_b),score)
    return final

   

