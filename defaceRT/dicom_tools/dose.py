"""
RT Dose operations for DefaceRT.

This module handles RT dose data, including resizing and rescaling the mask to the dose map.
"""

import numpy as np


debug=False

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


def resample_dose_map_3D(dose_map, new_spacing, old_spacing,new_size):

    scaling_factors = [old/new for old, new in zip(old_spacing,new_spacing)]
    original_size = np.array([len(dose_map[0]), len(dose_map[0][0]),len(dose_map)]) # change to 3D
    if debug:
        print("OG SIZE",original_size)
        # new_size_xyz = np.round(original_size*scaling_factors).astype(int)
        # new_size = [new_size_xyz[2],new_size_xyz[0],new_size_xyz[1]]
        print("NEW SIZE",new_size)
        print("scaling factors",scaling_factors)
        print(original_size*scaling_factors)
    resampled_dose_map = np.zeros(new_size)
#     return resampled_dose_map
    if debug:
        print(resampled_dose_map.shape)
    # resampled_dose_map[:new_size[2], :new_size[0], :new_size[1]] = dose_map[:new_size[2], :new_size[0], :new_size[1]]

#     # to do make 3d # to do -- this is exctremely slow.
#     for k in range(new_size[0]):
#         for i in range(new_size[1]):
#             for j in range(new_size[2]):

#                 original_index = (np.array([i,j,k])/scaling_factors).astype(int)
# #                 print(original_index)
# #                 print(int(original_index[2]), int(original_index[0]),int(original_index[1]))
#                 resampled_dose_map[k][i][j] = dose_map[ int(original_index[2]), int(original_index[0]),int(original_index[1])]            
    grid_k, grid_i, grid_j = np.meshgrid(np.arange(new_size[0]), np.arange(new_size[1]), np.arange(new_size[2]), indexing='ij') 
    original_indices = (np.array([grid_i, grid_j, grid_k]).transpose(1, 2, 3, 0) / scaling_factors).astype(int)
    # print(original_indices)
      # Debugging: Print the maximum values of original indices
    if debug:
        print("Max original indices:", np.max(original_indices, axis=(0, 1, 2)))

    # Ensure indices are within bounds
    original_indices = np.clip(original_indices, 0, [len(dose_map[0]) - 1, len(dose_map[0][0]) - 1, len(dose_map) - 1])
    if debug:
        print("clip to:", [len(dose_map) - 1, len(dose_map[0][0]) - 1, len(dose_map[0]) - 1])
    # Debugging: Print the modified original indices
    # print("Clipped original indices:", original_indices)
    resampled_dose_map = dose_map[original_indices[..., 2], original_indices[..., 0], original_indices[..., 1]]
    return resampled_dose_map


def resize_dose_map(dose_map,new_size, spacing, new_origin, old_origin,default=0):

    resized_dose_map = np.zeros(new_size)

    x_img = int((new_origin[0]-old_origin[0])/spacing[0])
    y_img = int((new_origin[1]-old_origin[1])/spacing[1])
    if debug:
        print(x_img,y_img)
        print((new_origin[0]-old_origin[0])/spacing[0],(new_origin[1]-old_origin[1])/spacing[1])
        print("X",x_img, len(dose_map[0]))
        print("Y",y_img, len(dose_map))
        print(x_img+len(dose_map[0]))
    
    y_end = y_img+len(dose_map)
    x_end = x_img+len(dose_map[0])
    

    if debug:
        print('xy ends',x_end,y_end)
    if y_end > 512:
        dose_map = dose_map[:(new_size[1]-y_end)]
    if x_end > 512:
        dose_map = dose_map[:,:(new_size[0]-x_end)]
          
    resized_dose_map[y_img:y_end,x_img:x_end] = dose_map
    return resized_dose_map

def resize_dose_map_3D(dose_map,new_size, spacing, new_origin, old_origin,default=0):
    
    new_size_zxy = new_size[2],new_size[0],new_size[1]
    # print(new_size_zxy)
    resized_dose_map = np.zeros(new_size_zxy)
    z_image = int((new_origin[2] - old_origin[2])/spacing[2])
    if debug:
        print("z_img",z_image, (new_origin[2] - old_origin[2])/spacing[2])
    
    len_dose_map = len(dose_map) 
    if debug:
        print("len dose map",len_dose_map)
    # Crop dose map if starting index is negative
    if z_image < 0:       
        z_image = 0
        len_dose_map = len_dose_map + z_image
    # print("len dose map after < 0",len_dose_map)     
    
    for i,resized_index in enumerate(range(z_image,z_image+len_dose_map)):  # added zimage + lendose map. idk anymore
        if debug:
            print(i,resized_index)
            print("new or (DOSE)", new_origin[2])
            print("old or (IMG)", old_origin[2])
            print("spacing", spacing)
        resized_dose_map[resized_index] = resize_dose_map(dose_map[i],[new_size[0],new_size[1]],spacing,new_origin,old_origin,default=0)
    
    return resized_dose_map



def reverse_resize_dose_map_3D(dose_map,new_size, spacing, new_origin, old_origin,default=0):
    new_size_zxy = new_size[0],new_size[1],new_size[2]
    if debug:
        print(new_size_zxy)
    resized_dose_map = np.zeros(new_size_zxy)
    z_image = int((new_origin[2] - old_origin[2])/spacing[2])
    if debug:
        print("z_img",z_image, (new_origin[2] - old_origin[2])/spacing[2])
    # 
    len_dose_map = len(dose_map)
    if debug:
        print("len dose map",len_dose_map)
    # Crop dose map if starting index is negative
    if z_image < 0:       
        z_image = 0
        len_dose_map = len_dose_map + z_image
    if debug:
        print("len dose map after < 0",len_dose_map)    
    
    # z_start = int((new_origin[2]-old_origin[2])/spacing[2])

    
    for i,resized_index in enumerate(range(z_image,z_image+new_size[0]-1)):  # added zimage + lendose map. idk anymore
        # print(i,resized_index)
        # if i > lennew_size[0]
        if debug:
        
            print(i,resized_index)
            print("new or (DOSE)", new_origin[2])
            print("old or (IMG)", old_origin[2])
            print("spacing", spacing)
            print("len dose map resized:",len(resized_dose_map))
            print("len dose map or?:", len(dose_map))
            print()
        #if dose map is smaller than resized dose map, leave with zeros?
        if i < len(dose_map):
            try:
                resized_dose_map[resized_index] = reverse_resize_dose_map(dose_map[i],[new_size[1],new_size[2]],spacing,new_origin,old_origin)
            except Exception as e:

                print("Leaving zero in dose map at index", resized_index,"because :",e)
    return resized_dose_map

def reverse_resize_dose_map(dose_map,new_size, spacing,  old_origin,new_origin):
    resized_dose_map = np.zeros(new_size)

    x_start = int((new_origin[0]-old_origin[0])/spacing[0])
    y_start = int((new_origin[1]-old_origin[1])/spacing[1])
  
    
    y_end = y_start+new_size[0]#len(dose_map)
    x_end = x_start+new_size[1]#len(dose_map[0])
        
    resized_dose_map = dose_map[y_start:y_end,x_start:x_end]
    


    return resized_dose_map