"""
DICOM input/output operations for DefaceRT.

This module handles reading and writing DICOM files, including CT scans,
RT SStructure Sets, and RT Dose.
"""

import os
import pydicom as dcm
from datetime import datetime


def get_RS(patient_path, CT_file):
    '''
    Get the RT Structure Set (RS) for a given CT scan.

    :param patient_path: path to patient directories.
    :param CT_file: name of the CT file (assumes it contains an RT Struct file starting with 'RS').
    
    :returns: tuple of the opened RS file and the filename
    '''
    image_path = os.path.join(patient_path,CT_file)
    # CT = dcm.dcmread(image_path)
    RS_file = [f for f in os.listdir(image_path) if f[0:2] == 'RS'][0]
    RS = dcm.dcmread(os.path.join(image_path,RS_file))

    #TO DO:load dose and RP fiels
    return RS,RS_file


def get_first_CT(patient_path):
    
    CT_list = [d for d in os.listdir(patient_path) if d[9:11] == 'CT' and len(d) == 23]
    
    CT_list.sort(key=lambda x: datetime.strptime(x[12:], "%d_%b_%Y"))
    # print(CT_list[0])
    return CT_list[0]

def get_CT_list(patient_path,CT_keyword='',CT_dir_name_min_len=0,CT_dir_name_max_len=100,ignore_terms=[]):
    '''
    Get list of CT scan directories matching criteria.
    
    :param patient_path: Path to patient directory
    :param CT_keyword: Keyword that must be present in CT name
    :param CT_dir_name_min_len: Minimum directory name length
    :param CT_dir_name_max_len: Maximum directory name length
    :param ignore_terms: Terms to exclude from CT names (ignores directories containing these keywords)
        
    :Returns: List of CT directory names
    '''
    CT_list = sorted([d for d in os.listdir(patient_path) if CT_keyword in d and len(d) >= CT_dir_name_min_len and len(d) <= CT_dir_name_max_len and all(substring.lower() not in d.lower() for substring in ignore_terms) ])
    
    # CT_list.sort(key=lambda x: datetime.strptime(x[12:], "%d_%b_%Y"))
    return CT_list


def find_dose_file(CT_path):
    '''
    Find and load dose file from CT directory.
    
    :param CT_path: Path to CT directory

    :returns: loaded Dose file
    '''
    dose_files = [f for f in os.listdir(CT_path) if 'RD' in f]
    num_dose_files = len(dose_files)
    
    if num_dose_files == 0:
        raise FileNotFoundError("ERROR: NO DOSE FILES FOUND")
    
    RD = dcm.dcmread(CT_path+'/'+dose_files[0]) 
    
    if num_dose_files == 1:
        return RD
    
    elif num_dose_files > 1:
        found = True
        smallest_spacing = RD[0x0028, 0x0030]
        most_frames = RD.NumberOfFrames
        
        for dose_file in dose_files[1:]:
            rd = dcm.dcmread(CT_path+'/'+dose_file)
            spacing = rd[0x0028, 0x0030]
            frames = rd.NumberOfFrames
            
            if spacing[0] < smallest_spacing[0] and spacing[1] < smallest_spacing[1]:
                smallest_spacing = spacing
                RD = rd 
            elif frames > most_frames:
                most_frames = frames
                RD = rd
    return RD

