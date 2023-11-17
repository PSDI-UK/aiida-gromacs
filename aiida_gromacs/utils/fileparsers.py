#!/usr/bin/env python
"""
Functions for parsing various gromacs input/output files and extracting
metadata into a dictionary.
"""
import re
import json

def nested_dict_extractor(i, j, lines, input_dict, split_with, leading_space, 
            leading_space_check_list):
    """
    For each of the input parameers defined in the gromacs log file, output
    each indented section into dictionary format.
    """
    for k, line3 in enumerate(lines[i+j+1:]):
        if line3 == '\n':
            break
        leading_space_check = re.search('\S', line3).start()
        if ':' in line3 and leading_space_check in leading_space_check_list:
            break
        if split_with in line3:
            l = line3.strip().split(split_with)
            if leading_space_check > leading_space:
                input_dict[l[0].strip()] = l[1].strip()
    return input_dict

def gromacs_logfile_parser(handle):
    input_params = {}
    averages = {}
    with open(handle, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if re.match(r"(?i)Executable:", line):
                for j, line2 in enumerate(lines[i:]):
                    if line2 == "\n":
                        break
                    leading_space = re.search('\S', line2).start()
                    if ":" in line2 and leading_space == 0:
                        top = line2.strip().split(':')[0].strip()
                        val = line2.strip().split(':')[1].strip()
                        input_params[top] = val
            if "Command line:" in line:
                if i + 1 < len(lines):
                    command = lines[i+1].strip()
                    input_params["Command line"] = rf"{command}"
            if re.match(r"(?i)GROMACS version:", line):
                for j, line2 in enumerate(lines[i:]):
                    if line2 == "\n":
                        break
                    leading_space = re.search('\S', line2).start()
                    if ":" in line2 and leading_space == 0:
                        top = line2.strip().split(':')[0].strip()
                        val = line2.strip().split(':')[1].strip()
                        input_params[top] = val
            if re.match(r"(?i)Input\sParameters", line):
                for j, line2 in enumerate(lines[i:]):
                    if line2 == "\n":
                        break
                    leading_space = re.search('\S', line2).start()
                    if ":" in line2 and leading_space == 0:
                        top = line2.strip().split(':')[0]
                        input_params[top] = {}
                        nested_dict_extractor(i, j, lines, input_params[top], 
                                "=", leading_space, [0,3])
                    if ":" in line2 and leading_space == 3:
                        top2 = line2.strip().split(':')[0]
                        input_params[top][top2] = {}
                        nested_dict_extractor(i, j, lines, input_params[top][top2], 
                                "=", leading_space, [3,5])
                    if ":" in line2 and leading_space == 5:
                        top3 = line2.strip().split(':')[0]
                        input_params[top][top2][top3] = {}
                        nested_dict_extractor(i, j, lines, 
                                input_params[top][top2][top3], "=", leading_space, 
                                [5,7])
                    if ":" in line2 and leading_space == 7:
                        top4 = line2.strip().split(':')[0]
                        input_params[top][top2][top3][top4] = {}
                        nested_dict_extractor(i, j, lines, 
                                input_params[top][top2][top3][top4], 
                                "=", leading_space, [7])
            if "A V E R A G E S" in line:
                for j, line2 in enumerate(lines[i:]):
                    if "Statistics" in line2:
                        l = line2.strip().split()
                        averages["total-steps"] = l[2]
                        averages["total-frames"] = l[5]
                    if "M E G A - F L O P S" in line2:
                        break
                    numbers = re.findall(r'\d+\.\d+e[\+|-]\d+', line2, re.DOTALL)
                    if len(numbers) != 0:
                        possible_header = lines[i+j-1].strip() #remove \n
                        if re.match(r"[a-z,A-Z]", possible_header):    
                            header = re.split('\s{2}+', possible_header)
                            header = list(filter(None, header)) # remove '' entries          
                            if len(numbers) == len(header):
                                for hn in range(len(header)):
                                    averages[header[hn].strip()] = numbers[hn]
            if "Time:" in line:
                averages["Time"] = {}
                l = line.strip().split()[1:]
                head = list(filter(None, lines[i-1].strip().split('  ')))
                for hn in range(len(header)):
                    averages["Time"][head[hn]] = l[hn]
            if "Performance:" in line:
                averages["Performance"] = {}
                l = line.strip().split()[1:]
                head = list(filter(None, lines[i-1].strip().split('  ')))
                for hn in range(len(header)):
                    averages["Performance"][head[hn]] = l[hn]
                        
    averages = {"Summary": averages}
    # merge dicts
    all_dict = input_params | averages
    return all_dict
