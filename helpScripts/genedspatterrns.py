import sys
import random

def next(curSeg,pattern,length,segments):
    # print(pattern,length)

    if length<=0 or curSeg >= len(segments):
        return ''

    if curSeg % 2 :
        #   block
        print("block",curSeg)
        curSeq = random.randint(0,len(segments[curSeg])-1)
        print(curSeq)
        if segments[curSeg][curSeq] == '':
            return next(curSeg+1,pattern,length,segments)
        elif length > len(segments[curSeg][curSeq]):
            return next(curSeg+1,pattern + segments[curSeg][curSeq],length - len(segments[curSeg][curSeq]),segments)
        else:
            return pattern+segments[curSeg][curSeq][:length]
    else:
        #   seed
        # print("seed",curSeg)
        # print(length > len(segments[curSeg]),segments[curSeg],length)
        if length > len(segments[curSeg]):
            # print("next")
            return next(curSeg+1,pattern + segments[curSeg],length - len(segments[curSeg]),segments)
        else:
            # print("stop")
            return pattern+segments[curSeg][:length]
    # return pattern
        

def run(length,number,segments):
    rand_start=0
    rand_loc=0     
    patterns = list()
    # for i in segments:
        # print(i)
    ranSegments = [random.randint(0,len(segments)-1) for _ in range(number)]

    pattern = ''
    for curSeg in ranSegments:
        # print("curSeg",curSeg)
        if curSeg  % 2 :
            # print("block")
            curSeq = random.randint(0,len(segments[curSeg])-1)
            curPos = random.randint(0,len(segments[curSeg][curSeq])-1) if len(segments[curSeg][curSeq]) != 0 else 0
            # print(curSeq,curPos)
            pattern = next(curSeg+1,segments[curSeg][curSeq][curPos:],length-len(segments[curSeg][curSeq][curPos:]),segments )
        else:
            # print("seed")
            curPos = random.randint(0,len(segments[curSeg])-1)
            # print(curPos)
            if length > len(segments[curSeg][curPos:]):
                pattern = next(curSeg+1, segments[curSeg][curPos:],length-len(segments[curSeg][curPos:]),segments)
            else:
                pattern = segments[curSeg][:length]

        if pattern == '':
            ranSegments.append(random.randint(0,len(segments)-1)) 
            continue
        if len(pattern) < length:
            ranSegments.append(random.randint(0,len(segments)-1)) 
        if len(pattern) > length:
            # print(pattern)
            exit
        else:
            patterns.append(pattern)

    return patterns

def __main__(args):

    input = args[1]
    length = int(args[2])
    number = int(args[3])
    output = args[4]

    # print(input,length,number,output)

    with open(input,'r') as input_file:
        segments = input_file.readlines()[0].replace('}','{').split('{')

    if segments[-1] == '':
        segments = segments[:-1]

    for i in range(len(segments)):
        segments[i] = segments[i].split(',') if ',' in segments[i] else segments[i]

    # print(segments)

    patterns = run(length,number,segments)

    # print(patterns)

    with open(output,'w') as output_file:
        output_file.write(f'# number={number} length={length} file={output}'+'\n')
        output_file.write(''.join(patterns))



#   run main
__main__(sys.argv)
