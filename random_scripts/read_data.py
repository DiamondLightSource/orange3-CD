#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 13:32:36 2026

@author: ubx84221
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from natsort import natsorted



def test_empty_line(line):
    # if the line is only commas
    if len(''.join(line.strip().split(','))) == 0:
        return None
    else:
        return ', '.join(line.strip().split(','))

def parse_remarks(lines):
    """
    parse the remarks section at the start of a csv file
    """
    comments = {}
    for line in lines:
        tokens = line.split(":")
        
        if len(tokens) == 1:
            empty_line = test_empty_line(line)
            if empty_line is not None:
                comments["Header"] = empty_line            
        else:
            try:
                assert len(tokens) == 2
                entry_key = tokens[0].split('#')[1].strip()
                entry_item = tokens[1].strip().split(",")[0]
                comments[entry_key] = entry_item
            except AssertionError:
                comments[line.strip('#')] = None

    return comments

def parse_data(lines):
    """
    Parse the data section of an input csv file. 
    
    Return a dict of {property: pd.DataFrame}
    """
    data_dict = {}
    sections = [idx for idx, i in enumerate(lines) if "Wavelength" in i]
    for section in range(len(sections)-1):
        section_title = lines[sections[section]-1].split(",")[0].strip()
        section_data = lines[sections[section]+2:sections[section+1]-2]
        data = np.array([[float(i) for i in line.strip().split(",") if len(i)>0] for line in section_data])
        
        data_av = data[:,1:].mean(axis=1)
        
        data_out = np.hstack((data[:,1:], data_av[:,None]))

        df = pd.DataFrame(data = data_out, index=data[:,0],
                          columns=[str(i) for i in range(data[:,1:].shape[1])]+["Average"])
        df.index.name = "Wavelength"
        
        data_dict[section_title] = df
        
    return data_dict
    
def default_parse(lines):
    return None
            

def file_parser(file):
    """
    Parse a single file, return the data as appropriate
    """
    with open(file) as f:
        lines = f.readlines()
    
    """
    find the lines which are the beginnings of the sections
    """
    section_starts = [idx for idx, i in enumerate(lines) if ":" in i.split(",")[0] and "#" not in i.split(",")[0]]
    
    parsers = {"Remarks": parse_remarks,
               "Data": parse_data,
               }
    
    data = {}
    for idx in range(len(section_starts[:-1])):
        section_title = lines[section_starts[idx]]    
        needed = lines[section_starts[idx]+1:section_starts[idx+1]]
        section = parsers.get(section_title.split(":")[0], default_parse)(needed)
    
        data[section_title.split(":")[0]] = section

    return data


if __name__ == "__main":
    fs = natsorted(glob.glob("/dls/science/users/ubx84221/b23/test_data/*.csv"))


    solution_A = '/dls/science/users/ubx84221/b23/test_data/HKD_WT_B30SR_1_mg_mL_00000.csv'
    solution_B = '/dls/science/users/ubx84221/b23/test_data/dabrafenib_2_equiv00000.csv'
    buffer = '/dls/science/users/ubx84221/b23/test_data/buffer00000.csv'
    
    reference_files = [solution_A, solution_B, buffer]
    experiment_files = [i for i in fs if i not in reference_files]


    all_data = {i.split('/')[-1]: file_parser(i) for i in experiment_files}