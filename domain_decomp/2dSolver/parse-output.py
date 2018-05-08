#!/usr/bin/env python3
import argparse
import re

def parse_timings(f_iter,retval):
    while True:
        line = next(f_iter).strip('\n')
        try:
            if line != '':
                name = line
                while True:
                    try:
                        tokens = iter(next(f_iter).split())
                        curr_token = next(tokens)
                        if curr_token == 'average':
                            next(tokens)
                            retval[name] = float(next(tokens))
                            break
                    except StopIteration:
                        pass
        except StopIteration:
            pass



def read_file(file_name):
    retval = {}
    retval['Cores'] = re.findall('[0-9]+',file_name)

    with open(file_name) as f:
            content = f.readlines()

    f_iter = iter(content)
    while True:
        try:
            # get tokens for current line
            tokens = iter(next(f_iter).split())
            try:
                curr_token = next(tokens)

                if curr_token == 'Iterations:':
                    retval['Iterations'] = int(next(tokens))

                if curr_token == 'Error:':
                    retval['Error'] = float(next(tokens))

                if curr_token == 'Residual:':
                    retval['Residual'] = float(next(tokens))

                if curr_token == 'Total':
                    curr_token = next(tokens)
                    if curr_token == 'cells:':
                        retval['Degrees of freedom'] = int(next(tokens))

                if curr_token == 'TIMING':
                    next(f_iter)
                    parse_timings(f_iter,retval)

            except StopIteration:
                pass
        except StopIteration:
            break

    return retval

def refactor_cores(values):
    first_values=[]
    for num in values[0]['Cores']:
        first_values.append(int(num))

    last_values=[]
    for num in values[-1]['Cores']:
        last_values.append(int(num))

    divs = [b/a for a, b in zip(first_values,last_values)]
    core_index = divs.index(max(divs))

    for value in values:
        value['Cores']=int(value['Cores'][core_index])



# parse arguments
parser = argparse.ArgumentParser(description='Parse timing outputs.')
parser.add_argument('files' ,metavar='file', type=str, nargs='+', help='first output file')
args = parser.parse_args()

values = []
for curr_file in args.files:
     values.append(read_file(curr_file))


refactor_cores(values)

col_names = []
for key in values[0].keys():
    col_names.append(key)

print(', '.join(col_names))

for value in values:
    print(', '.join([str(s) for s in value.values()]))
