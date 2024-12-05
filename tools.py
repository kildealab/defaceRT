import os
from datetime import datetime
import pydicom as dcm
import numpy as np
from scipy import interpolate

def xyz_to_image_coords(X,Y,Z,spacing,origin):
    X_new, Y_new, Z_new = [], [], []
    for x,y,z in zip(X,Y,Z):
        X_new.append((x-origin[0])/spacing[0])
        Y_new.append((y-origin[1])/spacing[1])
        Z_new.append((z-origin[2])/spacing[2])
    
    return X_new, Y_new, Z_new


def get_uni_spline(xi,yi):
    
        # Fit spline with s=0 (passing through all points)
    tck = interpolate.UnivariateSpline(xi, yi, s=0)
        
#         # Evaluate spline for 1000 evenly spaced points
#         xj, yj = interpolate.splev(np.linspace(0, 1, 1000), tck)
    return tck
        
        
def reverse_pixels_hu(anon,scans):
#     image = np.stack([s.pixel_array for s in scans])
    # Convert to int16 (from sometimes int16), 
    # should be possible as values should always be low enough (<32k)
    image = anon
    image = image.astype(np.int16)
    
    # Convert to Hounsfield units (HU)
    intercept = scans[0].RescaleIntercept
    slope = scans[0].RescaleSlope
    print("Intersept:", intercept)
    print("slope:", slope)
    
    if slope != 1:
        image = image.astype(np.float64) / slope
        image = image.astype(np.int16)
        
    image -= np.int16(intercept)
    
    return np.array(image, dtype=np.int16)

def get_pixels_hu(scans):
    image = np.stack([s.pixel_array for s in scans])
    # Convert to int16 (from sometimes int16), 
    # should be possible as values should always be low enough (<32k)
    image = image.astype(np.int16)
    
    # Convert to Hounsfield units (HU)
    intercept = scans[0].RescaleIntercept
    slope = scans[0].RescaleSlope
    print("Intersept:", intercept)
    print("slope:", slope)
    
    if slope != 1:
        image = slope * image.astype(np.float64)
        image = image.astype(np.int16)
        
    image += np.int16(intercept)
    
    return np.array(image, dtype=np.int16)

def get_RS(patient_path, CT_file):
    image_path = os.path.join(patient_path,CT_file)
    # CT = dcm.read_file(image_path)
    RS_file = [f for f in os.listdir(image_path) if f[0:2] == 'RS'][0]
    RS = dcm.read_file(os.path.join(image_path,RS_file))

    #TO DO:load dose and RP fiels
    return RS,RS_file


def get_first_CT(patient_path):
    
    CT_list = [d for d in os.listdir(patient_path) if d[9:11] == 'CT' and len(d) == 23]
    
    CT_list.sort(key=lambda x: datetime.strptime(x[12:], "%d_%b_%Y"))
    print(CT_list[0])
    return CT_list[0]


def get_eye_contours(RS,start_x,start_y,start_z,z_spacing,pixel_spacing):
    eye_ROI_names = find_ROI_names(RS,keyword='eye') + find_ROI_names(RS,'globe') + find_ROI_names(RS,'orbit')
    print(eye_ROI_names)

    dict_contours, z_lists = get_all_ROI_contours(eye_ROI_names, RS)
    roi_slice, z_smg = get_avg_ROI_z_and_slice(z_lists)
    all_z_slices = []
    for z_list in z_lists:
        all_z_slices = sorted(all_z_slices + [z for z in z_list if z not in all_z_slices])

        img_slice = get_image_slice(start_z, z_smg, [0,0,z_spacing])
    print("Centre ROI slice:", img_slice)

    roi_x, roi_y = get_ROI_pixel_array(dict_contours[eye_ROI_names[0]][roi_slice],start_x,start_y,pixel_spacing)
    roi_x2, roi_y2 = get_ROI_pixel_array(dict_contours[eye_ROI_names[1]][roi_slice],start_x,start_y,pixel_spacing)

    y_eye = []
    roi_eye = dict_contours[eye_ROI_names[1]][roi_slice]
    for eye in eye_ROI_names:
        roi_eye = dict_contours[eye][roi_slice] 
        for i in range(0,len(roi_eye),3):
            y_eye.append(roi_eye[i+1])
    max_eye = max(y_eye)
    min_eye = min(y_eye)
    y_cutoff_roi = ((max_eye-min_eye)/2+min_eye)

    min_x = (min(min(roi_x),min(roi_x2)))
    max_x = (max(max(roi_x),max(roi_x2)))
    min_y = (min(min(roi_y),min(roi_y2)))
    max_y = (max(max(roi_y),max(roi_y2)))

    y_cutoff = ((max_y-min_y)/2+min_y)

    return y_cutoff, z_lists, y_cutoff_roi, z_smg


'''
THESE ARE ALL FROM SLICE_SELECTION
'''


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
    '''
    dict_contours = {}
    z_lists = []
#     print(list_ROIs)

    for roi_name in list_ROIs:
        # print(roi_name)
        all_contours_ROI = get_contour_from_ROI_name(roi_name, RS)
        if len(all_contours_ROI)==0 or len(all_contours_ROI[0])==0: # if contour is empty
            return {},[]
        all_contours_ROI = sorted(all_contours_ROI, key=lambda x: x[2])
        dict_contours[roi_name] = all_contours_ROI
        z_lists.append(sorted([float(roi[2]) for roi in all_contours_ROI]))
    return dict_contours, z_lists


def get_contour_from_ROI_name(ROI_name, RS):
    '''
    Gets the contour for the requested ROI_name in the RS file.
    '''
    for i, seq in enumerate(RS.StructureSetROISequence):
        if seq.ROIName == ROI_name:
            index = i
            break

    contour_coords = [] 
    
    if 'ContourSequence' in RS.ROIContourSequence[index]:
        
        for ROI_contour_seq in RS.ROIContourSequence[index].ContourSequence:
            contour_coords.append(ROI_contour_seq.ContourData) 
    else:
        print("Warning:",ROI_name,"has no contour sequence.")
            
    return contour_coords


def get_avg_ROI_z_and_slice(z_lists):
    '''
    Gets the average (ie middle) z-value of the set of contour slices.
    
    TO DO: divide by zero error encountered when submandibular gland contour doesn't exist.
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
    index_colour_dict = {}
    for ROI_contour_seq in RS.ROIContourSequence:
        index_colour_dict[ROI_contour_seq.ReferencedROINumber] = ROI_contour_seq.ROIDisplayColor
   
    name_colour_dict = {}
    for seq in RS.StructureSetROISequence:
        name_colour_dict[seq.ROIName] = index_colour_dict[seq.ROINumber]

    return name_colour_dict

def image_to_xyz_coords_single(X,spacing,origin):
    if type(X) is not list:
        X = [X]
    X_new = []
    for x in X:
        X_new.append(x*spacing+origin)
        
    return X_new
    
def get_image_slice(start_z, z_smg, spacing):
    img_slice = int((abs(start_z - z_smg)/spacing[2]))
    return img_slice

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



'''
THESE ARE ALL FROM REGISTRATION
'''


'''
THESE ARE ALL FROM anon notebook
'''


def get_ROI_slice(z_target,z_list):
    
    slice_num = 0
    slice_nums = []
    for z in z_list:
        if z == z_target:
            slice_nums.append(slice_num)
#             return slice_num
        slice_num +=1
    return slice_nums



#new dose stuff

def find_dose_file(CT_path):
    dose_files = [f for f in os.listdir(CT_path) if 'RD' in f]
    num_dose_files = len(dose_files)
    
    if num_dose_files == 0:
        raise FileNotFoundError("ERROR: NO DOSE FILES FOUND")
    
    RD = dcm.read_file(CT_path+'/'+dose_files[0]) 
    
    if num_dose_files == 1:
        return RD
    
    elif num_dose_files > 1:
        
        smallest_spacing = RD[0x0028, 0x0030]
        
        for dose_file in dose_files[1:]:
            rd = dcm.read_file(CT_path+'/'+dose_file)
            spacing = rd[0x0028, 0x0030]
            
            if spacing[0] <= smallest_spacing[0] and spacing[1] <= smallest_spacing[1]:
                smallest_spacing = spacing
                RD = rd 
    return RD




def resample_dose_map_2D(dose_map, new_spacing, old_spacing):

    scaling_factors = [old/new for old, new in zip(old_spacing,new_spacing)]
    original_size = np.array([len(dose_map), len(dose_map[0])]) # change to 3D
    
    new_size = np.round(original_size*scaling_factors).astype(int)
    resampled_dose_map = np.zeros(new_size)

    # to do make 3d
    for i in range(new_size[0]):
        for j in range(new_size[1]):

            original_index = np.array([i,j])/scaling_factors
            resampled_dose_map[i][j] = dose_map[int(original_index[0]), int(original_index[1])]            
        
    return resampled_dose_map


def resample_dose_map_3D(dose_map, new_spacing, old_spacing):

    scaling_factors = [old/new for old, new in zip(old_spacing,new_spacing)]
    original_size = np.array([len(dose_map[0]), len(dose_map[0][0]),len(dose_map)]) # change to 3D
    
    print("OG SIZE",original_size)
    new_size_xyz = np.round(original_size*scaling_factors).astype(int)
    new_size = [new_size_xyz[2],new_size_xyz[0],new_size_xyz[1]]
    print("NEW SIZE",new_size)
    print("scaling factors",scaling_factors)
    print(original_size*scaling_factors)
    resampled_dose_map = np.zeros(new_size)
#     return resampled_dose_map
    print(resampled_dose_map.shape)

    # to do make 3d
    for k in range(new_size[0]):
        for i in range(new_size[1]):
            for j in range(new_size[2]):
           
#                 print(k)
                original_index = (np.array([i,j,k])/scaling_factors).astype(int)
#                 print(original_index)
#                 print(int(original_index[2]), int(original_index[0]),int(original_index[1]))
                resampled_dose_map[k][i][j] = dose_map[ int(original_index[2]), int(original_index[0]),int(original_index[1])]            

    return resampled_dose_map


def resize_dose_map(dose_map,new_size, spacing, new_origin, old_origin,default=0):

    resized_dose_map = np.zeros(new_size)

    x_img = int((new_origin[0]-old_origin[0])/spacing[0])
    y_img = int((new_origin[1]-old_origin[1])/spacing[1])
    # print(x_img,y_img)
    # print((new_origin[0]-old_origin[0])/spacing[0],(new_origin[1]-old_origin[1])/spacing[1])
    # print("X",x_img, len(dose_map[0]))
    # print("Y",y_img, len(dose_map))
    # print(x_img+len(dose_map[0]))
    
    y_end = y_img+len(dose_map)
    x_end = x_img+len(dose_map[0])
    


    print('xy ends',x_end,y_end)
    if y_end > 512:
        dose_map = dose_map[:(new_size[1]-y_end)]
    if x_end > 512:
        dose_map = dose_map[:,:(new_size[0]-x_end)]
          
    resized_dose_map[y_img:y_end,x_img:x_end] = dose_map
    return resized_dose_map

def resize_dose_map_3D(dose_map,new_size, spacing, new_origin, old_origin,default=0):
    
    new_size_zxy = new_size[2],new_size[0],new_size[1]
    print(new_size_zxy)
    resized_dose_map = np.zeros(new_size_zxy)
    z_image = int((new_origin[2] - old_origin[2])/spacing[2])
    print("z_img",z_image, (new_origin[2] - old_origin[2])/spacing[2])
    
    len_dose_map = len(dose_map) 
    print("len dose map",len_dose_map)
    # Crop dose map if starting index is negative
    if z_image < 0:       
        z_image = 0
        len_dose_map = len_dose_map + z_image
    print("len dose map after < 0",len_dose_map)     
    
    for i,resized_index in enumerate(range(z_image,z_image+len_dose_map)):  # added zimage + lendose map. idk anymore
        print(i,resized_index)
        print("new or (DOSE)", new_origin[2])
        print("old or (IMG)", old_origin[2])
        print("spacing", spacing)
        resized_dose_map[resized_index] = resize_dose_map(dose_map[i],[new_size[0],new_size[1]],spacing,new_origin,old_origin,default=0)
    
    return resized_dose_map



def reverse_resize_dose_map_3D(dose_map,new_size, spacing, new_origin, old_origin,default=0):
    new_size_zxy = new_size[0],new_size[1],new_size[2]
    # print(new_size_zxy)
    resized_dose_map = np.zeros(new_size_zxy)
    z_image = int((new_origin[2] - old_origin[2])/spacing[2])
    # print("z_img",z_image, (new_origin[2] - old_origin[2])/spacing[2])
    
    len_dose_map = len(dose_map)
    # print("len dose map",len_dose_map)
    # Crop dose map if starting index is negative
    if z_image < 0:       
        z_image = 0
        len_dose_map = len_dose_map + z_image
    # print("len dose map after < 0",len_dose_map)    
    
    z_start = int((new_origin[2]-old_origin[2])/spacing[2])
    
    for i,resized_index in enumerate(range(z_image,z_image+new_size[0])):  # added zimage + lendose map. idk anymore
        # if i > lennew_size[0]
        # print(i,resized_index)
        # print("new or (DOSE)", new_origin[2])
        # print("old or (IMG)", old_origin[2])
        # print("spacing", spacing)
        resized_dose_map[resized_index] = reverse_resize_dose_map(dose_map[i],[new_size[1],new_size[2]],spacing,new_origin,old_origin)
    
    return resized_dose_map

def reverse_resize_dose_map(dose_map,new_size, spacing,  old_origin,new_origin):
    resized_dose_map = np.zeros(new_size)

    x_start = int((new_origin[0]-old_origin[0])/spacing[0])
    y_start = int((new_origin[1]-old_origin[1])/spacing[1])
  
    
    y_end = y_start+new_size[0]#len(dose_map)
    x_end = x_start+new_size[1]#len(dose_map[0])
        
    resized_dose_map = dose_map[y_start:y_end,x_start:x_end]
    


    return resized_dose_map