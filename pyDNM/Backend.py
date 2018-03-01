#!/usr/bin/env python3
import os, sys
def err_fh(f):
    """checks if a file exists"""
    if not os.path.isfile(f):
        print('FATAL ERROR: {} does not exist'.format(f))
        sys.stderr.write('FATAL ERROR: {} does not exist\n'.format(f))
        sys.exit(1)
def tokenize(l):
    """splits a line by spaces or tabs, returns a list object"""
    r=l.rstrip().split('\t')
    if len(r)==1: r=l.rstrip().split(' ')
    if len(r)==1:
        sys.stderr.write('WARNING: {} not tokenized correctly. Check if {} is Tab or Space delimited\n'.format(l,l))
        return 0
    else: return r
def try_index(l,ind):
    v=None
    try: v=l[ind]
    except IndexError: v=None
    return v