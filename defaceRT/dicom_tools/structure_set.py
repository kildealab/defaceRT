"""
RT Structure Set operations for DefaceRT.

This module handles RT structure sets, including finding ROI names,
extracting contours, and modifying structure data for anonymization.
"""

import os
import numpy as np

from dicom_tools.tools import get_image_slice

def find_ROI_names(RS, keyword=''):
    '''
    find_ROI_names  finds all contour names in RT Structure Set File containing keyword, 
                    while ignoring those containing 'nos' and 'z_'.
    
    :param RS: the RS file opened by pydicom
    :param keyword: The keyword to search the ROIs for. If blank returns all ROIs.
    
    :returns: list of ROI names containing keyword.
    '''
    ROI_names = []

    for seq in RS.StructureSetROISequence:
        roi_name = seq.ROIName
        if keyword.lower() in roi_name.lower() and 'nos' not in roi_name.lower() and 'z_' not in roi_name.lower():
            ROI_names.append(seq.ROIName)
    return ROI_names


def get_all_ROI_contours(list_ROIs,RS):
    '''
    Get dictionary of contours for each of the ROIs in list_ROIs.

    :param list_ROIs: List of ROI names to extract contours for.
    :param RS: The RT Structure Set file opened by pydicom.

    :returns: A dictionary of contours for each ROI and a list of z-coordinates.
    '''
    dict_contours = {}
    z_lists = []
#     print(list_ROIs)

    for roi_name in list_ROIs:
        # print(roi_name)
        all_contours_ROI = get_contour_from_ROI_name(roi_name, RS)
        if len(all_contours_ROI)==0 or len(all_contours_ROI[0])==0: # if contour is empty
            # return {},[]
            dict_contours[roi_name] = []
        else:
            all_contours_ROI = sorted(all_contours_ROI, key=lambda x: x[2])
            dict_contours[roi_name] = all_contours_ROI
            z_lists.append(sorted([float(roi[2]) for roi in all_contours_ROI]))
    return dict_contours, z_lists


def get_contour_from_ROI_name(ROI_name, RS):
    '''
    Gets the contour for the requested ROI_name in the RS file.

    :param ROI_name: The name of the ROI to extract contours for.
    :param RS: The RT Structure Set file opened by pydicom.

    :returns: A list of contour coordinates for the requested ROI.
    '''
    index = -1
    for i, seq in enumerate(RS.StructureSetROISequence):
        # print(seq.ROIName,ROI_name)
        if seq.ROIName == ROI_name:
            index = i
            break
    if index == -1:
        print("Warning: No contour named", ROI_name)
        return []

    contour_coords = [] 
    
    if 'ContourSequence' in RS.ROIContourSequence[index]:
        
        for ROI_contour_seq in RS.ROIContourSequence[index].ContourSequence:
            contour_coords.append(ROI_contour_seq.ContourData) 
    else:
        print("Warning:",ROI_name,"has no contour sequence.")
            
    return contour_coords


'''
THESE ARE ALL FROM SLICE_SELECTION
'''

def get_avg_ROI_z_and_slice(z_lists):
    '''
    Gets the average (ie middle) z-value of the set of contour slices.

    :param z_lists: A list of z-coordinates for each contour slice.

    :returns: The average z-coordinate and the slice index.
    '''
    z_avg = 0
    # print("Z:",z_lists)

    for z_list in z_lists:
        z_avg += (np.sum(z_list)/len(z_list))

    z_avg = z_avg/len(z_lists)
    roi_slice = np.argmin(abs(z_lists[0] - z_avg)) #TO DO: what if not in this list? but it should ebf ine
    z_smg = z_lists[0][roi_slice]

    return roi_slice, z_smg

def get_ROI_colour_dict(RS):
    '''
    Get a dictionary mapping ROI names to their display colours.

    :param RS: The RT Structure Set file opened by pydicom.

    :returns: A dictionary mapping ROI names to their display colours.
    '''
    index_colour_dict = {}
    for ROI_contour_seq in RS.ROIContourSequence:
        index_colour_dict[ROI_contour_seq.ReferencedROINumber] = ROI_contour_seq.ROIDisplayColor
   
    name_colour_dict = {}
    for seq in RS.StructureSetROISequence:
        name_colour_dict[seq.ROIName] = index_colour_dict[seq.ROINumber]

    return name_colour_dict


def get_eye_contours(RS,start_x,start_y,start_z,z_spacing,pixel_spacing):
    '''
    Get the contours for the eye regions from the RT Structure Set.

    :param RS: The RT Structure Set file opened by pydicom.
    :param start_x, start_y, start_z: Image origin coordinates.
    :param z_spacing: Slice thickness.
    :param pixel_spacing: The spacing between pixels in the x and y directions.

    :returns: A tuple containing:
            - y-cutoff: y-coord cutoff for image masking
            - z-lists: List of z-coordinates for each eye contour slice
            - y-cutoff ROI: y-coord cutoff in ROI space
            - average z: Middle eye z-coordinate
    '''
    eye_ROI_names = sorted(find_ROI_names(RS,keyword='eye') + find_ROI_names(RS,'globe') + find_ROI_names(RS,'orbit'), key=len)
    print(eye_ROI_names)
    if len(eye_ROI_names) == 0:
        print("ERROR: NO EYES.")
        print(find_ROI_names(RS))
        raise Exception("No eye contours found.")
   

    dict_contours, z_lists = get_all_ROI_contours(eye_ROI_names, RS)
    roi_slice, z_smg = get_avg_ROI_z_and_slice(z_lists)
    all_z_slices = []
    for z_list in z_lists:
        all_z_slices = sorted(all_z_slices + [z for z in z_list if z not in all_z_slices])

        img_slice = get_image_slice(start_z, z_smg, [0,0,z_spacing])
    print("Centre ROI slice:", img_slice)

    roi_x, roi_y = get_ROI_pixel_array(dict_contours[eye_ROI_names[0]][roi_slice],start_x,start_y,pixel_spacing)

    if len(eye_ROI_names) > 1: # In case only one eye contoured
        roi_x2, roi_y2 = get_ROI_pixel_array(dict_contours[eye_ROI_names[1]][roi_slice],start_x,start_y,pixel_spacing)

    y_eye = []
    roi_eye = dict_contours[eye_ROI_names[0]][roi_slice]
    for eye in eye_ROI_names:
        try:
            roi_eye = dict_contours[eye][roi_slice] 
            for i in range(0,len(roi_eye),3):
                y_eye.append(roi_eye[i+1])
        except:
            print("No contours for",eye, "at ROI slice",roi_slice)
    max_eye = max(y_eye)
    min_eye = min(y_eye)
    y_cutoff_roi = ((max_eye-min_eye)/2+min_eye)

    if len(eye_ROI_names) > 1:
            
        min_x = (min(min(roi_x),min(roi_x2)))
        max_x = (max(max(roi_x),max(roi_x2)))
        min_y = (min(min(roi_y),min(roi_y2)))
        max_y = (max(max(roi_y),max(roi_y2)))
    else:
        min_x = min(roi_x)
        max_x = max(roi_x)
        min_y = min(roi_y)
        max_y = max(roi_y)

    y_cutoff = ((max_y-min_y)/2+min_y)

    return y_cutoff, z_lists, y_cutoff_roi, z_smg

def get_ROI_pixel_array(roi_array,start_x,start_y,pixel_spacing):
    '''
    Get the pixel positions (rather than the x,y coords) of the contour array so it can be plotted.
    '''
#     roi_array = dict_contours[subgland_ROI_names[0]][roi_slice]
    x = []
    y = []

    for i in range(0,len(roi_array),3):
        x.append((roi_array[i] - start_x)/pixel_spacing[0])
        y.append((roi_array[i+1] - start_y)/pixel_spacing[1]) 
        
    return x, y


def get_ROI_slice(z_target,z_list):
    
    slice_num = 0
    slice_nums = []
    for z in z_list:
        if z == z_target:
            slice_nums.append(slice_num)
#             return slice_num
        slice_num +=1
    return slice_nums

