import sys

input_file = sys.argv[1]

with open(input_file,'r') as input:
    with open(input_file.split('.')[0]+'.txt','w') as out:
        for i in input.readlines():
            out.write(i.replace('-',''))