# -*- coding: utf-8 -*-

import argparse
import subprocess
import os
import requests as r
import xml.etree.ElementTree as ET
from tqdm import tqdm
import sys
import glob

#   GLOBALS
# allowed_fasta_extenstions = ['.fa','.fasta','.fna']
def transform(*args, **kwargs):
    file = args[0]
    ls = args[1]
    method = args[2]
    executable_path = f"{os.path.dirname(__file__)}/scripts/build/src/transform"

    if not ls:
        print(f"Transforming {file} into regular EDS")
        print(file)
    else:
        # print(f"Transforming {file} into l-EDS with l:{ls} using {method} method")
        for l in tqdm(ls):
            out_file = file[:-3] + str(l)
            out_file += '.l.leds' if method == 'linear' else '.c.leds'
            print(file, out_file,l)
            p = subprocess.Popen([executable_path,f"-i{file}", f"-o{out_file}", f"-l{l}" ,f"--method={method}"], stdout=subprocess.PIPE)
            # out, err = p.communicate()
            # print(err)

def stats(*args, **kwargs):
    files = args[0]
    executable_path = f"{os.path.dirname(__file__)}/scripts/build/src/stats"
    for file in files:
        p = subprocess.Popen([executable_path,file], stdout=subprocess.PIPE)
        out, err = p.communicate()
        yield out.decode(),file

def build(*args, **kwargs):
    files = args[0]
    executable_path = f"{os.path.dirname(__file__)}/scripts/build/src/build"
    for file in files:
        index_dir = file+'.index'
        print("Building",file,index_dir)
        p = subprocess.Popen([executable_path,file,index_dir], stdout=subprocess.PIPE)
        out, err = p.communicate()
        yield out.decode(),index_dir
        # yield "out",index_dir

def locate(*args, **kwargs):
    print(f"#       Locating...")
    if kwargs["type"] == 'SEQUENCE':
        print("... as sequence")
    else:
        print("... as file")

def parse_range(value):
    try:
        if value == None:
            return list()
        nums = [int(i) for i in value.split('-')]
        if len(nums) == 1:
            return nums
        else:
            return range(*nums)
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Invalid input: {value}. Must be a single number or a range in the format 'start-end'.")

def parse_inputs(inputs, type='f'):
    files = []
    sequences = []

    for input_path in inputs:
        input_path = Path(input_path)

        # Check if the input is a directory
        if input_path.is_dir():
            if type == 'f':
                # Collect all files under the directory
                files += [str(file) for file in input_path.rglob('*') if file.is_file()]
            elif type == 'd':
                # Collect all directories under the directory
                files += [str(dir) for dir in input_path.rglob('*') if dir.is_dir()]
        else:
            # If it's not a directory, consider it a sequence
            sequences.append(str(input_path))

    return sequences, files

def main():
    #   parse arguments
    parser = argparse.ArgumentParser(description="Create and index Elastic-Degenerat String")
    # parser.add_argument("command", help="The command to execute (e.g., build, download, get, find, minimize, kernelize).")
    subparsers = parser.add_subparsers(dest="command", help="Command to execute {stats/transform/build/locate}")

    parser_print_stats = subparsers.add_parser("stats", help="Get statistics about given EDS")
    parser_print_stats.add_argument("inputs", nargs='+', type=str, help="input in format of set of sequences, set of files or folder that contains given original files")
    parser_print_stats.add_argument("-o","--output", type=argparse.FileType('w'), default=sys.stdout, help="output file")
    parser_print_stats.add_argument("--csv", required=False, default=False, action='store_true', help="if minimized file already exists, it will be replaced by the new one")

    parser_transform = subparsers.add_parser("transform", help="Create EDS or lEDS from MSA or VCF")
    parser_transform.add_argument("-l", required=False, type=parse_range, help="context length, can be given range in format a-b-c as from a to b by c steps")
    parser_transform.add_argument("inputs", nargs='+', type=str, help="input in format of set of sequences, set of files or folder that contains given original files")
    parser_transform.add_argument("-o","--output", type=argparse.FileType('w'), default=sys.stdout, help="output file")
    parser_transform.add_argument("--method", required=None ,type=str, help="")

    # Build command
    parser_build = subparsers.add_parser("build", help="Build index")
    parser_build.add_argument("inputs", nargs='+', type=str, help="input in format of set of sequences, set of files or folder that contains given original files")
    parser_build.add_argument("--rebuild", required=False, default=False, action='store_true', help="if minimized file already exists, it will be replaced by the new one")

    # Locate command
    parser_locate = subparsers.add_parser("locate", help="find MEMs of given patterns with respect to the tree data")
    parser_locate.add_argument("index", nargs='+', type=str, help="path to the index directory")
    parser_locate.add_argument("-p", "--pattern", nargs='+', type=str, help="Patterns")
    parser_locate.add_argument("-P", "--pattern_files", nargs='+', type=str, help="A pattern file (one pattern per line)")
    # parser_locate.add_argument("--rebuild", required=False, default=False, action='store_true', help="if output file already exists, it will be replaced by the new one")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == "locate":
        raise Exception("Not Implemented!")
    elif args.command == "build":
        raise Exception("Not Implemented!")
    elif args.command == "stats":
        # print(f"#Printin statistics about given eds files on {args.inputs}")
        files = []
        for input in args.inputs:
            files += glob.glob(input)

        out_stream = args.output
        columns = "file,n_common,l_common,n,N,m,min_l,max_l,avg_l,n_empty_strings"
        if args.csv:
            out_stream.write(columns+'\n')
        for res,file in stats(files):
            if args.csv:
                out_stream.write(file+',')
                out_stream.write(','.join([i.split(':')[1] for i in res.split('\n')[:-1]])+'\n')
            else:
                out_stream.write(res+'\n')
    elif args.command == "transform":
        allowed_extensions = ["msa","vcf","eds"]
        allowed_methods = ["linear","cartesian"]
        method = args.method
        if not method and args.l:
            print("No method specified - will be used linear as default")
            method = "linear"

        if method not in allowed_methods:
            raise Exception(f"Unknown method: {args.method}, must be one of {allowed_methods}")

        # executable_path = "./scripts/build/src/transform"
        _,files = parse_inputs(args.inputs,type='f')
        for file in files: 
            ext = file.split('.')[-1] 
            if ext not in allowed_extensions:
                raise Exception(f"Unknown file extension: {ext}, supports {allowed_extensions}")

            transform(file,args.l,args.method)           
    else:
        raise Exception(f"Unknown command: {args.command}, check your arguments")
if __name__ == "__main__":
    main()
