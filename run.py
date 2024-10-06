# -*- coding: utf-8 -*-

import argparse
import subprocess
import os
import requests as r
import xml.etree.ElementTree as ET
from tqdm import tqdm
import sys
import json

#   GLOBALS
# allowed_fasta_extenstions = ['.fa','.fasta','.fna']
def transform(*args, **kwargs):
    file = args[0]
    ls = args[1]
    method = args[2]
    executable_path = "./scripts/build/src/transform"

    if not ls:
        print(f"Transforming {file} into regular EDS")
        print(file)
    else:
        print(f"Transforming {file} into l-EDS with l:{ls} using {method} method")
        for l in tqdm(ls):
            print(file,l)
            # p = subprocess.Popen([executable_path,"-i {file} -l"], stdout=subprocess.PIPE)
            # out, err = p.communicate()

def stats(*args, **kwargs):
    files = args[0]
    executable_path = "./scripts/build/src/stats"
    for file in files:
        p = subprocess.Popen([executable_path,file], stdout=subprocess.PIPE)
        out, err = p.communicate()
        yield out.decode(),file

def build(*args, **kwargs):
    files = args[0]
    executable_path = "./scripts/build/src/build"
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

def parse_inputs(inputs,type='f'):
    files = list()
    sequences = list()
    for input in inputs:
        try:
            out = subprocess.run([f'find {input} -type {type}'], capture_output = True, text = True, shell=True)
            out = out.stdout
            if out == '':
                sequences.append(input)
            else:
                files += out.split('\n')[:-1]
        except Exception as e:
            raise Exception(f"Input parsing unsucessful due to {e}")
    return sequences,files

def main():
    #   parse arguments
    parser = argparse.ArgumentParser(description="Create and index Elastic-Degenerat String")
    # parser.add_argument("command", help="The command to execute (e.g., build, download, get, find, minimize, kernelize).")
    subparsers = parser.add_subparsers(dest="command", help="Command to execute {stats/transform/build/locate}")

    parser_print_stats = subparsers.add_parser("stats", help="Get statistics about given EDS")
    parser_print_stats.add_argument("inputs", nargs='+', type=str, help="input in format of set of sequences, set of files or folder that contains given original files")
    parser_print_stats.add_argument("--csv", required=False, default=False, action='store_true', help="if minimized file already exists, it will be replaced by the new one")

    parser_transform = subparsers.add_parser("transform", help="Create EDS or lEDS from MSA or VCF")
    parser_transform.add_argument("-l", required=False, type=parse_range, help="context length, can be given range in format a-b-c as from a to b by c steps")
    parser_transform.add_argument("inputs", nargs='+', type=str, help="input in format of set of sequences, set of files or folder that contains given original files")
    parser_transform.add_argument("--method", required=None ,type=str, help="")

    # Build command
    parser_build = subparsers.add_parser("build", help="Build index")
    parser_build.add_argument("inputs", nargs='+', type=str, help="input in format of set of sequences, set of files or folder that contains given original files")
    parser_build.add_argument("--rebuild", required=False, default=False, action='store_true', help="if minimized file already exists, it will be replaced by the new one")

    # Locate command
    parser_locate = subparsers.add_parser("locate", help="find MEMs of given patterns with respect to the tree data")
    parser_locate.add_argument("index", nargs='+', type=str, help="path to the index directory")
    parser_locate.add_argument("-p", "--pattern", nargs='+', type=str, help="Single pattern or a pattern file (one pattern per line)")
    # parser_locate.add_argument("--rebuild", required=False, default=False, action='store_true', help="if output file already exists, it will be replaced by the new one")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == "locate":
        print(f"#Locate in indices {args.index} patterns {args.pattern}")
        _, indices_files = parse_inputs(args.index,type='d')
        sequences,patterns_files = parse_inputs(args.pattern,type='f')
        if sequences != list():
            print(f"#    Total {len(indices_files)} indices has been found")
        else:
            raise Exception("No index file, check your arguments")
        if patterns_files != list():
            print(f"#    Total {len(patterns_files)} pattern files has been found")
        if sequences != list():
            print(f"#    Total {len(sequences)} patterns has been found")

        for index in indices_files:
            for pattern in tqdm(patterns_files):
                out = locate(index,pattern,type='FILE')
                print(out)

            for pattern in tqdm(sequences):
                out = locate(index,pattern,type='SEQUENCE')
                print(out)
    elif args.command == "build":
        print(f"#Building indices on {args.inputs}")
        _,files = parse_inputs(args.inputs,type='f')
        for out,index_dir in tqdm(build(files)):
            print(out)
            print(out)
    elif args.command == "stats":
        # print(f"#Printin statistics about given eds files on {args.inputs}")
        _,files = parse_inputs(args.inputs,type='f')

        columns = "file, n_common, l_common, n, N, m, min_l, max_l, avg_l, n_empty_strings"
        if args.csv:
            print(columns)
        for out,file in stats(files):
            if args.csv:
                print(file,end=',')
                print(','.join([i.split(':')[1] for i in out.split('\n')[:-1]]))
            else:
                print(out)
    elif args.command == "transform":
        allowed_extensions = ["msa","vcf","eds"]
        allowed_methods = ["linear","cartesian"]
        method = args.method
        print(method)

        if not method and args.l:
            print("No method specified - will be used linear as default")
            method = "linear"

        if method not in allowed_methods:
            raise Exception(f"Unknown method: {args.method}, must be one of {allowed_methods}")

        executable_path = "./scripts/build/src/transform"
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
