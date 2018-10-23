from operator import add
from functools import reduce
import re,textwrap,csv

def translate_dna(sequence):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
    'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    proteinsequence = ''

    for n in range(0,len(sequence),3):
        if sequence[n:n+3] in codontable:
            proteinsequence += codontable[sequence[n:n+3]]
        else:   
            proteinsequence += ''
    return proteinsequence

def atgc(sequence):
    freq = [0,0,0,0]
    for i in sequence:
        if i == 'A':
            freq[0]+=1
        elif i == 'T':
            freq[1]+=1
        elif i == 'G':
            freq[2]+=1
        elif i == 'C':
            freq[3]+=1
        else:
            continue
    return freq

def tm(freq):
    return ( 2*(freq[0]+freq[1]) + 4*(freq[2]+freq[3]))

seq_dict = {'A':'T','T':'A','G':'C','C':'G','*':'*','Z':'Z'}

def primer(sequence,lmno):
    l = []
    flag = 0
    cunt = 0
    for j in range(18,30):
        flag = 0
        for i in range(0, len(sequence)):
            sum = 0
            if (j+i > len(sequence)):
                break

            seq = (sequence)[i:j+i]
            user = atgc(seq)
            t = tm(user)
            if (t > 60 and t < 70):
                sum +=1

            ratio = (user[2]+user[3]) / (user[0] + user[1] + user[2] + user[3])
            if (ratio >= .40 and ratio <= .60):
                adity= 1
            else:
                adity= 0
            sum += adity
            if not (re.compile(r'([CG]{4,30})').match(seq)):
                sum+=1
            if not (re.compile(r'([A]{5,30})').match(seq) or re.compile(r'([T]{5,30})').match(seq) or re.compile(r'([G]{5,30})').match(seq) or re.compile(r'([C]{5,30})').match(seq)):
                sum+=1
            if (seq[-1] == 'A'):
                sum+=1
            elif (seq[-1] == 'C'):
                sum+=1
            if ((user[1]+user[0]) / (user[0] + user[1] + user[2] + user[3]) >= .50):
                sum+=1
            if not (seq[0:9] == ("".join([seq_dict[base] for base in reversed(seq[-10:-1])]))[::-1]):
                sum+= 1
            if (sum == 7):
                l.append(seq)
                cunt +=1

            if (cunt == lmno+2):
                flag = 1
                break
        
        if (flag):
            break
    flag = 0
    
    k = l[lmno]
    for j in range(18,30):
        flag = 0
        for i in range(0, len(sequence)):
            sum = 0
            seq = ("".join([seq_dict[base] for base in reversed(sequence)]))[i:j+i]
            if (j+i > len(sequence)):
                break
            user = atgc(seq)
            t = tm(user)
            if (t > 60 and t < 70):
                sum +=1

            ratio = (user[2]+user[3]) / (user[0] + user[1] + user[2] + user[3])
            if (ratio >= .40 and ratio <= .60):
                adity= 1
            else:
                adity= 0
            sum += adity
            if not (re.compile(r'([CG]{4,30})').match(seq)):
                sum+=1
            if not (re.compile(r'([A]{5,30})').match(seq) or re.compile(r'([T]{5,30})').match(seq) or re.compile(r'([G]{5,30})').match(seq) or re.compile(r'([C]{5,30})').match(seq)):
                sum+=1
            if (seq[-1] == 'A'):
                sum+=1
            elif (seq[-1] == 'C'):
                sum+=1
            if ((user[1]+user[0]) / (user[0] + user[1] + user[2] + user[3]) >= .50):
                sum+=1
            if not (seq[0:9] == ("".join([seq_dict[base] for base in reversed(seq[-10:-1])]))[::-1]):
                sum+= 1
            jpg = tm(atgc(k))
            if (sum >=5):
                if (abs(jpg-t) < 5):
                    string = "Forward Primer - " + str(k) + "\t" + "Tm(in deg C) - " + str(jpg) + "\n" + "Reverse Primer - " + str(seq) + "\t" + "Tm(in deg C) - " + str(tm(user)) + "\n"
                    open4_file.write(string)
                    flag = 1
                    break
        if (flag):
            break
    return 

in_file1 = open("sequence.fasta","r")
s1 = in_file1.read()
in_file1.close()
in_file2 = open("Introns.csv")
s2 = list(csv.reader(in_file2))
in_file2.close()
s2.pop(0)
s2 = reduce(add, s2, [])

l2=['*']

info = []

l1 = s1.split(">")
l1.pop(0)
sequence = []
for i in l1:
    info.append((i.split("\n"))[0])
    l3 = i.split("\n")
    l3.pop(0)
    sequence.append(''.join(l3))

isequence = []
temp=''
for i in sequence:
    l3 = []
    temp = ''
    if not (int(s2[0])-1 < 0):
        temp = i[0:int(s2[0])-1]
        for j in range(0, len(s2),2):
            for l in range(int(s2[j]),int(s2[j+1])+1):
                temp += '*'
            if (j+2 < len(s2)):
                temp += i[int(s2[j+1]):int(s2[j+2])-1]
            if (j+2 >= len(s2)):
                break
        if (int(s2[-1]) < len(i)):
            temp += i[int(s2[j-1]):-1]
            temp += i[-1]
        break
isequence.append(temp)
refrence_frame1 = []
refrence_frame4 = []

for i in range(0,3):
    refrence_frame1.append(sequence[0][i:]) 
    refrence_frame4.append(isequence[0][i:])

for i in range(0,3):
    refrence_frame1.append(("".join([seq_dict[base] for base in reversed(sequence[0])]))[i:])
    refrence_frame4.append(("".join([seq_dict[base] for base in reversed(isequence[0])]))[i:])

counter = 0

open4_file = open("4.txt","w")
for i in refrence_frame1:
    primer(str(i),counter)
    counter+=1

refrence_frame2 = []
for i in refrence_frame4:
    stseq = ''
    count = 0
    for j in i:
        count+=1
        stseq += j 
        if (count%3 == 0):
            stseq += 'Z'
    refrence_frame2.append(stseq)
pmer = (("".join([seq_dict[base] for base in reversed(refrence_frame2[0])]))[::-1])
pmer = pmer.replace("*","")
pmer = pmer.replace("Z","")
refrence_frame7 = []
for i in refrence_frame1:
    stseq = ''
    count = 0
    for j in i:
        count+=1
        stseq += j 
        if (count%3 == 0):
            stseq += 'Z'
    refrence_frame7.append(stseq)
counter = 0
open1_file = open("1.txt","w")
for i in refrence_frame7:
    counter += 1
    slr = "> Reading Frame " + str(counter) + "\n"
    open1_file.write(slr)
    open1_file.write(textwrap.fill(i.replace("Z",""),70))
    slr= "\n\n"
    open1_file.write(slr)


counter = 0
length = 0
cunt  = 0
fsquirt = ''

pattern = re.compile(r'(ATG(?:[ATGC]{3})*?([ATGCZ*])*?(?:TAA|TAG|TGA))')
for i in refrence_frame2:
    a = i.replace("Z"," ")
    a = a.replace("*","X")
    for key1 in range(0,len(i),4):
        if (i[key1:key1+3] == 'ATG'):
            break;
    for key2 in range(len(i)-1,-1,-1):
        if (i[key2] != '*' and i[key2] != 'Z'):
            break;
    i = i[key1:key2+1]
    a = i.replace("Z","")
    a = a.replace("*","")
    stseq = ''
    count = 0
    for mem in a:
        count+=1
        stseq += mem 
        if (count%3 == 0):
            stseq += 'Z'
    a = stseq
    a = list(pattern.findall(a))
    counter += 1
    lona = 0
    for j in a:
        j = list(j)
        j.pop(-1)
        j[0] = j[0].replace("Z","")
        j[0] = j[0].replace("*","X")
        if (length < len(j[0])):
            cunt  = counter
            fsquirt = j[0]
            length = len(j[0])

open2_file = open("2.txt","w")
slr = "> Reading Frame " + str(cunt) + "\n"
open2_file.write(slr)
open2_file.write(textwrap.fill(refrence_frame7[cunt-1].replace("Z",""),70))
open2_file.write("\n")

open3_file = open("3.txt","w")
slr= "> Protein Sequence\n"
open3_file.write(slr)
open3_file.write(textwrap.fill(translate_dna(fsquirt),70))
open3_file.write("\n")

open5_file = open("5.txt","w")
slr= "> E.coli cloning sequence \n"
open5_file.write(slr)
open5_file.write(textwrap.fill(fsquirt,70))
slr="\n\n> Pichia pastoris cloning sequence \n"
open5_file.write(slr)
open5_file.write(textwrap.fill(refrence_frame7[cunt-1].replace("Z",""),70))
slr="\n\n> HEK293 cloning sequence \n"
open5_file.write(slr)
open5_file.write(textwrap.fill(refrence_frame7[cunt-1].replace("Z",""),70))
slr="\n\n"
open5_file.write(slr)
	