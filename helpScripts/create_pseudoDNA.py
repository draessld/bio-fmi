import sys
from random import choices, choice, randint
from functools import reduce

number_of_sequences, avg_sequence_length, change_pbt = sys.argv[-3:]

number_of_sequences = int(number_of_sequences)
avg_sequence_length = int(float(avg_sequence_length)*1000000) // number_of_sequences
change_pbt = float("0."+change_pbt)
# print(change_pbt)

alphabet = "ACGTN"
psts = [0.24,0.24,0.24,0.24,0.04]

chances = [change_pbt * i for i in [0.8,0.2]]
sequences = []
sequence = ""

sigma = round(avg_sequence_length / 10)

def gen_avg(expected_avg,n,a,b):
    while True:
        l = [randint(a,b) for i in range(n)]
        avg = reduce(lambda x,y: x+y,l)/len(l)

        if avg == expected_avg:
            return l


sequence_lengths = gen_avg(avg_sequence_length,number_of_sequences,avg_sequence_length - sigma,avg_sequence_length + sigma) 

#   create reference
for _ in range(sequence_lengths[0]):
    sequence += choices(list(alphabet),psts)[0]
for _ in range(5):
    random_loc = randint(0,sequence_lengths[0])
    random_length = randint(0,sigma)
    sequence = sequence[:random_loc] + '-'*random_length + sequence[random_loc:]

sequences.append(sequence)
reference_length = len(sequence)

for sequence_id in range(1,number_of_sequences):
    length = 0
    sequence = ""
    for i in range(reference_length):
        change_bool = choices([0,1],[1-change_pbt,change_pbt])[0]
        if change_bool:
            #   create change
            action = choices([0,1],chances)[0]
            if action == 0:
                #   SNP
                sequence += choices(list(alphabet),psts)[0]
                length +=1
            else:
                #   INDEL
                if sequences[0][i] != '-':
                    #   INSERT
                    sequence += choices(list(alphabet),psts)[0]
                    length +=1
                else:
                    #   DELETE
                    sequence += '-'
        else:
            sequence += sequences[0][i]
            if sequences[0][i] != '-':
                length+=1
        if length >= sequence_lengths[sequence_id]:
            break
    if length < sequence_lengths[sequence_id]:
        for i in range(sequence_lengths[sequence_id]-length):
            sequence += choices(list(alphabet),psts)[0]
            length +=1
        for i in range(len(sequences)):
            sequences[i] += '-'*(len(sequence)-len(sequences[i]))

    sequence += '-' * (len(sequences[0])-len(sequence))

    sequences.append(sequence)

for i,l in zip(sequences,sequence_lengths):
    print(i)