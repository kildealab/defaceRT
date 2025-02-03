
import os, time, sys
# from datetime import datetime
import pydicom as dcm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cv2

# sys.path.append('../')
# from Slice_Selection.slice_selection import *
# from tools import find_ROI_names, get_first_CT, get_RS, get_pixels_hu, get_eye_contours, get_all_ROI_contours, get_contour_from_ROI_name, get_avg_ROI_z_and_slice
from tools import *
from defaceRT.plotting import plot_3_views, plot_all_contours






# to do, fix z mess for when its revesrsed
def generate_anon_image(image,z_lists,spacing, start_z,y_cutoff, cut_above = True,isodose=[], contours_to_keep = {},origin=[0,0,0],reverse_z=False):
    anon_image = image.copy()
    if len(z_lists) > 1:
        all_z_slices = sorted(z_lists[0] + [z for z in z_lists[1] if z not in z_lists[0]])
    else:
        all_z_slices = sorted(z_lists[0])
    img_slices = [get_image_slice(start_z, z, spacing) for z in all_z_slices]
    # print("iamge slices",img_slices)
    
    
    z_min = min(img_slices)
    z_max = max(img_slices)
    if z_max >= image.shape[0]:
        z_max = image.shape[0]-1
    
    if cut_above:
        if reverse_z:
            z_min = 0
        else:
            z_max = len(image)-1

    # print("image_size:",image.shape)
    mask = np.full_like(image,1)
    mask[z_min:z_max+1,0:int(y_cutoff)+1,:] = 0
    # print("mask size:",mask.shape)


    for contour in contours_to_keep:
        for c_slice in contours_to_keep[contour]:
            slice_num = get_image_slice(start_z, c_slice[2], spacing)
            # print("slice num:",slice_num)
            # print(z_min,z_max)
            if slice_num >= z_min and slice_num <= z_max:
                xi, yi, zi = c_slice[::3], c_slice[1::3], c_slice[2::3]
                X,Y,Z = xyz_to_image_coords(xi,yi,zi,spacing,origin)

                cont_xy = np.array([tuple([int(x),int(y)]) for x,y in zip(X,Y)])
                cont_xy = cont_xy.reshape((-1,1,2))
                cv2.drawContours(mask[slice_num],[cont_xy],-1,1, -1)
                #TO DO IUSE THIS TOCROP OTEHR CINTIURS

  # TO DO: ISODOSE BELOW      
    # print("anon*MASK************")    
    # print((anon_image*mask)[0])
    # print(mask)
    # print(mask*-1+1)
    # print(anon_image)
    return anon_image*mask++(mask*-1+1)*(-1000), mask

def reshape_contours_in_region():

    for contour in contours_to_keep:
        for c_slice in contours_to_keep[contour]:
            slice_num = get_image_slice(start_z, c_slice[2], spacing)
            if slice_num >= z_min and slice_num <= z_max:
                xi, yi, zi = c_slice[::3], c_slice[1::3], c_slice[2::3]
                X,Y,Z = xyz_to_image_coords(xi,yi,zi,spacing,origin)

                cont_xy = np.array([tuple([int(x),int(y)]) for x,y in zip(X,Y)])
                cont_xy = cont_xy.reshape((-1,1,2))
                cv2.drawContours(mask[slice_num],[cont_xy],-1,1, -1)

def get_contour_curves(slices, all_z_slices, y_cutoff_roi):
    dict_curves = {}
    
    for c_slice in slices:
#         print(c_slice[2] )
        z = c_slice[2]
        if z in all_z_slices or z>min(all_z_slices):
#             print(z)
            
            

            xo = []
            yo = []
            x_prev = -1000

            xi, yi, zi = c_slice[::3], c_slice[1::3], c_slice[2::3]
            for k,y_i in enumerate(yi):
            #         print(xi[k])
                if y_i < y_cutoff_roi and x_prev != xi[k]:
                    xo.append(xi[k])
                    yo.append(yi[k])
                    x_prev = xi[k]
                    
            
            if len(xo)> 3 and len(yo) > 3:
                
                

                xy = [[xs,ys] for xs, ys in sorted(zip(xo, yo))]
                # print("xy",xy)
                xs, ys = zip(*xy)
                # print("xs",xs)
#                 x_prev = -1000
                xs = list(xs)
                ys = list(ys)
                
                ind_to_del = []
                
#                 i = 1
                for i,x in enumerate(xs):
                    # print(i,x)
                    # print(x, xs[i-1])
                    if x == xs[i-1]:
                        # print("YRD")
                        if  ys[i] > ys[i-1]:
                            ind_to_del.append(i-1)

                        else:
                            ind_to_del.append(i)
#                     i+=1
                
                offset =0
                for ind in ind_to_del:
                    del xs[ind-offset]
                    del ys[ind-offset]
                    offset+=1
                    
#                         print(len(xs))
#                         if  ys[i] > ys[i-1]:
#                             del ys[i-1]
#                             del xs[i-1]

#                         else:
#                             del ys[i]
#                             del xs[i]
#                         i-=1
                   
    

                        
                # print(xs)
                # print(ys)

                s =  get_uni_spline(xs, ys)
    #             print(s)
                dict_curves[z] = {}
                dict_curves[z]['spline'] = s
                dict_curves[z]['x_max'] = max(xs)
                dict_curves[z]['x_min'] = min(xs)
    return dict_curves



def generate_anon_body(dict_contours_body,z_lists,y_cutoff, isodose=[], body_name = 'BODY',contours_to_keep=[],get_ind = False,reverse_z=False):
    if len(dict_contours_body[body_name]) == 0:
        return None

    if len(z_lists) > 1:
        all_z_slices = sorted(z_lists[0] + [z for z in z_lists[1] if z not in z_lists[0]])
    else:
        all_z_slices = sorted(z_lists[0])
    # print(body_name)
    edited = False
#     x_iso = X_new
#     y_iso = Y_new
#     z_iso = Z_new

#TODO:add isodoes
    contour_curves = {}
  
    for contour in contours_to_keep:
        contour_curves[contour] = get_contour_curves(contours_to_keep[contour], all_z_slices, y_cutoff)
        # print(contour_curves)
 
        
    
    x_iso = []
    y_iso = []
    z_iso = []
    
    
    # test = dict_contours[body_ROI_names[0]][roi_slice]
    # print(len(test))
    XN = []
    YN = []
    ZN = []
    fullN = []
    full_stack_N = []
    y_prev = 0
    x_prev = 0
    z_prev = -1000
    z_cont = -1000
    for c in dict_contours_body[body_name]:
    #     print(c)
        test = c
    #     full = full + list(c)
        slice = []

#         if c[2] < -650:
#             continue

        for i in range(0,len(c),3):

            x = test[i]
            y = test[i+1]
            z = test[i+2]
#             print(x)


            add_N = True
            if z in all_z_slices or (not reverse_z and z>min(all_z_slices)) or (reverse_z and z<min(all_z_slices)):
                if y < y_cutoff:
                    # if z==137:
                    #     print(y)
                    y = y_cutoff
#                     print(z,z_prev)
#                     if True:#z != z_prev:
#                         iso_z_ind =  [ind for ind, e in enumerate(z_iso) if e == z]
#                         print(z_iso)
#                     for index in iso_z_ind:
#                         if round(x_iso[index]) == round(x) and y_iso[index] < y:
#                             y = y_iso[index]
                    edited = True
                    for contour in contours_to_keep:
                        if z in contour_curves[contour]:  
                            contour_slice = contour_curves[contour][z]
                            if x > contour_slice['x_min'] and  x < contour_slice['x_max']:
                                y = contour_slice['spline'](x)
                                


    #         if y < min_eye + 30:
    #             print(y, max_y)
#             print(add_N)
            if not (y== y_prev and x==x_prev):# and add_N:
                XN.append(x)
                YN.append(y)
                ZN.append(z)
                fullN.append([x,y,z])
                slice = slice + [x,y,z]
            
            y_prev = y
            x_prev = x
            z_prev = z

        # if z==137:
        #     print(slice)
        full_stack_N.append(slice)
    if edited:
        return full_stack_N
    else:
        return None

def generate_anon_dose(RD, mask,CT_spacing,CT_origin,CT_size,z_image_slice):
    dose_array = RD.pixel_array
    original_dose_size = dose_array.shape

    # no need to scale to Gy since just applying mask, but to scale: dose_array * RD.DoseGridScaling
    dose_spacing = [RD.PixelSpacing[0],RD.PixelSpacing[1],RD.GridFrameOffsetVector[1]-RD.GridFrameOffsetVector[0]]



    scaling_factors = [old/new for old, new in zip(dose_spacing,CT_spacing)]
    original_size = np.array([len(dose_array[0]), len(dose_array[0][0]),len(dose_array)]) # change to 3D
    
    print("OG SIZE",original_size)
    new_size_xyz = np.round(original_size*scaling_factors).astype(int)
    new_size = [new_size_xyz[2],new_size_xyz[0],new_size_xyz[1]]

    # img_slice_resized = get_image_slice(RD.ImagePositionPatient[2], z_image_slice, CT_spacing)
    img_slice_resized = get_image_slice(RD.ImagePositionPatient[2], z_image_slice, dose_spacing)
    print("Slice resiszed",img_slice_resized)
    print("SPACING ct dose:", CT_spacing, dose_spacing)



    # rewrite to instead just resize/resample mask to dose map -- might lose info if changing dose map to image and back
    resized_mask = reverse_resize_dose_map_3D(mask,new_size,CT_spacing, CT_origin,RD.ImagePositionPatient)
    print("*******")
    print(original_dose_size,CT_spacing, CT_origin,RD.ImagePositionPatient)
    print(resized_mask.shape)
    print("*******")
    resampled_mask = resample_dose_map_3D(resized_mask, dose_spacing,CT_spacing,original_dose_size)
    defaced_dose = dose_array * resampled_mask

    plt.subplot(2,2,4)
    if len(defaced_dose) > img_slice_resized:
        plt.tick_params(axis='y', which='both', labelleft=False, labelright=True)
        plt.gca().yaxis.set_ticks_position('right')
        plt.imshow(defaced_dose[img_slice_resized],cmap=plt.cm.plasma)

        # plt.colorbar(label="")
        

        # move label
        # plt.ylabel('Your label here', labelpad=-725, fontsize=18)

    else:
        for arr in reversed(defaced_dose):
            # print(arr)
            if arr.any():

                plt.imshow(arr,cmap=plt.cm.plasma)
                plt.title("Defaced slice not in dose map.")
                break


    return defaced_dose

    # resampled_dose = resample_dose_map_3D(dose_array,CT_spacing, dose_spacing) 
    # resized_dose = resize_dose_map_3D(resampled_dose,CT_size,CT_spacing, RD.ImagePositionPatient,CT_origin)
    # anon_dose = resized_dose * mask

    # original_dose_size = dose_array.shape
    # desized_dose = reverse_resize_dose_map_3D(anon_dose,original_dose_size,CT_spacing, CT_origin,RD.ImagePositionPatient)
    # desampled_dose = resample_dose_map_3D(desized_dose, dose_spacing,CT_spacing)

    # plt.subplot(2,2,4)
    # plt.imshow(desampled_dose[img_slice],cmap=plt.cm.plasma)
    # return desampled_dose


def init_data(patient_path, CT_file, RS):

    slices = [dcm.read_file(patient_path+'/'+CT_file + '/'+ s) for s in os.listdir(patient_path+'/'+CT_file) if s[0:2] =='CT']
    CT_dcm = slices[0]
    # print(CT_dcm)
    # Order slices
    slices.sort(key = lambda x: (x.InstanceNumber))
    image = get_pixels_hu(slices)
    reverse_z = (slices[0].ImagePositionPatient[2]-slices[1].ImagePositionPatient[2] > 0) # Z increasung

    # Image Sacping Parameters
    origin = slices[0].ImagePositionPatient
    start_z = origin[2]
    start_x = origin[0]
    start_y = origin[1]
    z_spacing = slices[0].SliceThickness 
    pixel_spacing = slices[0].PixelSpacing
    spacing = [pixel_spacing[0],pixel_spacing[1],z_spacing]

    # Spacing in Z,X,Y format
    # dose_spacing = [RD.GridFrameOffsetVector[1]-RD.GridFrameOffsetVector[0],RD.PixelSpacing[0],RD.PixelSpacing[1]]
    CT_spacing = [CT_dcm.SliceThickness,pixel_spacing[0],pixel_spacing[1]]

    return slices,image,reverse_z, start_x, start_y, start_z, pixel_spacing, z_spacing, spacing, origin




def anon_dicom(anon,slices,slope=1.0,intercept=-1000):
    # print(len(anon), len(slices))
    slices_copy = slices.copy()
    # anon_rev = reverse_pixels_hu(anon, slices)
    # slope = 1.0
    # intercept = -1000
    # print(slices_copy[0].pixel_array)
    # print("anon")
    # print(anon[0])
    # print("anon-inter")
    # print(anon[0] - intercept)
    # print(slices_copy[0].RescaleIntercept)
    # print(anon_rev[0])
    for i,slice in enumerate(slices_copy):
        slices_copy[i].PixelData = ((anon[i]-intercept)).tobytes() #NOTE: divide by SLOPE IS FUCKED

    return slices_copy


def save_dicom(slices,save_path,patient, CT_file):
    for slice in slices:
        if not os.path.exists(save_path+patient):
            os.makedirs(save_path+patient)
        if not os.path.exists(save_path+patient+'/'+CT_file):
            os.makedirs(save_path+patient+'/'+CT_file)
        slice.save_as(save_path+patient+'/'+CT_file+'/'+'CT.'+slice.SOPInstanceUID+'.dcm')
    
def save_RT_struct(RS, RS_file,contour_stack,save_path,patient,CT_file,body_name='BODY'):
    #TODO: make it choose roi name, putting body for now
    if contour_stack == []:
        print("NO BODY RS TO SAVE")
    else:
        for i, seq in enumerate(RS.StructureSetROISequence):
            if seq.ROIName == body_name:
                    index = i
                    break
        
        dict_stack = {}
        z_s = []
        z_refs = []
        
        # print("CONTOUR STACK HERE")
        # print(contour_stack)
        # print("END OF CONTOUR STACK")

        for row in contour_stack:
            z_s.append(row[2])
            dict_stack[row[2]] = row
            # if row[2] == -545:
                # print("******************************************")
                # print("z-s after cropping")
                # print(row)
                # print("******************************************")

        z_s_done =[]
        for i, ROI_contour_seq in enumerate(RS.ROIContourSequence[index].ContourSequence):
            z_ref = ROI_contour_seq.ContourData[2] 
            z_refs.append(z_ref)
            # if z_ref == -545:
                # print("******************************************")
                # print("z-refs")
                # print(ROI_contour_seq.ContourData)
                # print(ROI_contour_seq)
                # print("******************************************")
            # if z_ref not in list(dict_stack.keys()):
            #     print("******",z_ref,"not in dict")

            replacement_contour = []
            for row in contour_stack:
                if float(z_ref) == float(row[2]):
                    print("YES")
                    # if int(row[2]) == int(-497):
                    #     print(row)
                    #     print("^^^^^497")
                    replacement_contour = row
                    contour_stack.remove(row)
                    z_s_done.append(z_ref)
                    continue
            if len(replacement_contour)!=0:
                # print("******",z_ref,"not in dict")
                # print()

            # else:
                RS.ROIContourSequence[index].ContourSequence[i].ContourData = replacement_contour
            # RS.ROIContourSequence[index].ContourSequence[i].ContourData = dict_stack[z_ref]

    # remove eyes, lenses, corneas, more?
    indices_remove = []
    for i, seq in enumerate(RS.StructureSetROISequence):
        if any(substring in seq.ROIName.lower() for substring in ('eye','globe','orbit','lens','cornea')):
            # print(seq.ROIName, i)
            indices_remove.append(i)

    for index in indices_remove:
        # print(index)
        for i, ROI_contour_seq in enumerate(RS.ROIContourSequence[index].ContourSequence):
            RS.ROIContourSequence[index].ContourSequence[i].ContourData = []

                
    # print(sorted(z_s))
    # print(sorted(z_refs))
    # print(sorted(z_s_done))
    # z_prev = sorted(z_s)[0]
    # for z in sorted(z_s[1:]):
    #     if z - z_prev !=3 :
    #         print("z ERROR",z, z_prev)
    #     z_prev = z

    # z_prev = sorted(z_refs)[0]
    # for z in sorted(z_refs[1:]):
    #     if z - z_prev !=3 :
    #         print("z_ref ERROR",z, z_prev)
    #     z_prev = z

    # print("len dict", len(dict_stack.keys()))
    # print("len zs", len(z_s))
    # print("len z_refs", len(z_refs))
    # print("len z_s_done", len(z_s_done))
    #todotoday - uncomment below
    # RS.save_as(os.path.join(save_path,patient,CT_file,RS_file))
        # contour_coords.append(ROI_contour_seq.ContourData)
    return RS

def save_RT_struct_all(RS, RS_file,contour_stack_dict,save_path,patient,CT_file):
    #TODO: make it choose roi name, putting body for now
    for name in contour_stack_dict.keys():
        body_name = name
        contour_stack = contour_stack_dict[name]
        if contour_stack == []:
            print("NO NEW STRUCTURE TO SAVE")
        else:
            for i, seq in enumerate(RS.StructureSetROISequence):
                if seq.ROIName == body_name:
                        index = i
                        break
            
            dict_stack = {}
            z_s = []
            z_refs = []

            for row in contour_stack:
                z_s.append(row[2])
                dict_stack[row[2]] = row
                # if row[2] == -545:
                    # print("******************************************")
                    # print("z-s after cropping")
                    # print(row)
                    # print("******************************************")

            z_s_done =[]
            for i, ROI_contour_seq in enumerate(RS.ROIContourSequence[index].ContourSequence):
                # print(i, len(ROI_contour_seq.ContourData))
                if len(ROI_contour_seq.ContourData) >= 3:
                    z_ref = ROI_contour_seq.ContourData[2] 
                    # print(name,"cont seq",i,':z=',z_ref)
                    z_refs.append(z_ref)
                    # if z_ref == -545:
                        # print("******************************************")
                        # print("z-refs")
                        # print(ROI_contour_seq.ContourData)
                        # print(ROI_contour_seq)
                        # print("******************************************")
                    # if z_ref not in list(dict_stack.keys()):
                    #     print("******",z_ref,"not in dict")

                    replacement_contour = []
                    for row in contour_stack:
                        if float(z_ref) == float(row[2]):
                            # print("YES")
                            # if int(row[2]) == int(-497):
                            #     print(row)
                            #     print("^^^^^497")
                            replacement_contour = row
                            contour_stack.remove(row)
                            z_s_done.append(z_ref)
                            break
                            # continue
                    if len(replacement_contour)!=0:
                        # print("******",z_ref,"not in dict")
                    #     print()

                    # else:
                        RS.ROIContourSequence[index].ContourSequence[i].ContourData = replacement_contour
                    # RS.ROIContourSequence[index].ContourSequence[i].ContourData = dict_stack[z_ref]

        # remove eyes, lenses, corneas, more?
        indices_remove = []
        dict_names = {}
        for i, seq in enumerate(RS.StructureSetROISequence):
            if any(substring in seq.ROIName.lower() for substring in ('eye','globe','orbit','lens','cornea')):
                # print(seq.ROIName, i)
                dict_names[i] = seq.ROIName
                indices_remove.append(i)


        for index in indices_remove:
            # print(index)
            # error handling for empty roi contour sequence
          
            try:
                
                for i, ROI_contour_seq in enumerate(RS.ROIContourSequence[index].ContourSequence):
                    RS.ROIContourSequence[index].ContourSequence[i].ContourData = []
            except:
                print("Warning: No contour sequence for structure", dict_names[index])


                
    # print(sorted(z_s))
    # print(sorted(z_refs))
    # print(sorted(z_s_done))
    # z_prev = sorted(z_s)[0]
    # for z in sorted(z_s[1:]):
    #     if z - z_prev !=3 :
    #         print("z ERROR",z, z_prev)
    #     z_prev = z

    # z_prev = sorted(z_refs)[0]
    # for z in sorted(z_refs[1:]):
    #     if z - z_prev !=3 :
    #         print("z_ref ERROR",z, z_prev)
    #     z_prev = z

    # print("len dict", len(dict_stack.keys()))
    # print("len zs", len(z_s))
    # print("len z_refs", len(z_refs))
    # print("len z_s_done", len(z_s_done))
    #todotoday - uncomment below
    # RS.save_as(os.path.join(save_path,patient,CT_file,RS_file))
        # contour_coords.append(ROI_contour_seq.ContourData)
    return RS

def run_anonymization(PATH,patient,save_path,keywords_keep = [],CT_name='',produce_pdf=True):
    patient_path = os.path.join(PATH,patient)
    if len(CT_name)==0:
        CT_file = get_first_CT(patient_path)
    else:
        CT_file = CT_name
    RS,RS_file = get_RS(patient_path, CT_file)
    
    # print(RS.ROIContourSequence[0].ContourSequence)
    #TODO: fixx dreadful return below to global var
    slices,image,reverse_z, start_x, start_y, start_z, pixel_spacing, z_spacing, spacing,origin = init_data(patient_path, CT_file, RS)
    # print(image[0])
    origin = [start_x,start_y,start_z]
    y_cutoff,z_lists,y_cutoff_roi,z_smg =  get_eye_contours(RS,start_x,start_y,start_z,z_spacing,pixel_spacing)
    CT_spacing = [pixel_spacing[0],pixel_spacing[1],z_spacing]
    CT_size = [len(image[0]),len(image[0][0]),len(image)] #to do double check x y are correct positon

    img_slice = get_image_slice(start_z, np.mean(z_lists[0]), spacing)
    if len(keywords_keep) == 0:
        dict_contours_keep = {}
    else:
        list_names_keep = []
        for keyword in keywords_keep:
            list_names_keep = list_names_keep + find_ROI_names(RS,keyword=keyword)

        list_remove = []
        # note to make customized -- doesn't include structures starting with z and cases where it is OAR-PTV, which should refer to things outside the PTV
        for name in list_names_keep:
            if '-PTV' in name or 'nonptv' in name.lower() or name[0].lower()=='z':
                list_remove.append(name)

        # Note - separated as can't remove from list while looping through list
        for name in list_remove:
            list_names_keep.remove(name)
        print("Keeping:",list_names_keep)
        dict_contours_keep,_ = get_all_ROI_contours(list_names_keep,RS)

    # print(find_ROI_names(RS)

    anon, mask = generate_anon_image(image,z_lists,spacing,start_z,y_cutoff,reverse_z=reverse_z,contours_to_keep=dict_contours_keep,origin=origin)
    # print("anon")
    # print(anon[0])




    new_dicom = anon_dicom(anon, slices)
    all_contour_names = find_ROI_names(RS)
    contour_names = [c for c in all_contour_names if c not in dict_contours_keep.keys()]
    dict_contours_all,_ = get_all_ROI_contours(contour_names,RS)

    dict_new_contours = {}


    for contour_name in contour_names:
        full_stack_N = generate_anon_body(dict_contours_all,body_name = contour_name,z_lists=z_lists,y_cutoff=y_cutoff_roi,reverse_z=reverse_z,contours_to_keep=dict_contours_keep)
        # print(contour_name, full_stack_N)
        if full_stack_N != None:
            # print("to change",contour_name)
            dict_new_contours[contour_name] = full_stack_N

    RS_new = save_RT_struct_all(RS, RS_file,dict_new_contours,save_path,patient,CT_file)
    #todotoday - uncomment below
    # save_dicom(new_dicom,save_path,patient,CT_file)
    '''
    body_names = find_ROI_names(RS,'brainstem')
    print(body_names)
    if len(body_names) == 0:
        body_names = find_ROI_names(RS,'external')
    
    if len(body_names)==0:
        print("NO BODY CONTOUR FOUND")
        full_stack_N = []
    else:
        dict_contours_body,_ = get_all_ROI_contours([body_names[0]],RS)
        # print(len(dict_contours_body['BODY']))
        # print("******************")
        dict_new_contours = {}

        full_stack_N = generate_anon_body(dict_contours_body,body_name = body_names[0],z_lists=z_lists,y_cutoff=y_cutoff_roi,reverse_z=reverse_z,contours_to_keep=dict_contours_keep)
            # for slice in full_stack_N:
        # if slice[2] == -497:#z_smg:
        #     print(slice)
        #     break
    # x_new_body = slice[::3] 
    # y_new_body = slice[1::3] 
    # z_new_body =slice[2::3] 
    # x_new, y_new, z_new = xyz_to_image_coords(x_new_body,y_new_body,z_new_body,spacing,origin)

    # print(len(full_stack_N))
    #todotoday - uncomment below
    RS_new = save_RT_struct(RS, RS_file,full_stack_N,save_path,patient,CT_file,body_names[0])
    '''

    # print(RS)
    # dict_contours_body,_ = get_all_ROI_contours(['BODY'],RS)


    # full_stack_N = generate_anon_body(dict_contours_body,z_lists=z_lists,y_cutoff=y_cutoff_roi,reverse_z=reverse_z,contours_to_keep=dict_contours_keep)
    # for slice in full_stack_N:
    #     if slice[2] == z_smg:
    # #         print(slice)
    #         break
    # x_new_body = slice[::3] 
    # y_new_body = slice[1::3] 
    # z_new_body =slice[2::3] 
    # x_new, y_new, z_new = xyz_to_image_coords(x_new_body,y_new_body,z_new_body,spacing,origin)

    # generate_anon_image(image,z_lists,spacing, y_cutoff,reverse_z=False,contours_to_keep=dict_contours_keep,origin=origin)

    if produce_pdf:
        # if i % batch_size == 0: # Open PDF for appending 
        #     pdf_mode = 'w' if i == 0 else 'a' 
        #     pdf = PdfPages(pdf_filename, mode=pdf_mode) # Plot the 3 views 
        # print("RS")
        plot_3_views(slices, img_slice, anon, patient=patient + ' ' + CT_file) 
        # print("ploted NEW")
        plt.subplot(2, 2, 3)
        plot_all_contours(RS_new,anon,get_image_slice(start_z, np.mean(z_lists[0]),spacing),origin,spacing,legend=True)
            # Save the current figure to the PDF
        # produce_pdfs()
        # pdf.savefig() 
        # plt.close()
        
        # # produce_pdfs()
        # if (i + 1) % batch_size == 0: 
        #     pdf.savefig(plt.figure()) # Add a blank page p
        #     df.close() # Close the current batch


    try:
        RD = find_dose_file(patient_path+'/'+CT_file)
        anon_dose = generate_anon_dose(RD, mask,CT_spacing,origin,CT_size,np.mean(z_lists[0]))#img_slice)
    except Exception as e:
        plt.subplot(2,2,4)
        # plt.title(patient+" ERROR")
        plt.text(0.5,0.5,e,dict(ha='center',va='center',fontsize=14,color='blue'))


if __name__ == "__main__":
    start = time.time()
    
    from Anon_config import *

    PATH = config["path"]
    keywords_keep = config["contours_to_keep"]
    save_path = config["save_path"]
    CT_name = config["CT_name"]
    produce_pdf = config['print_PDFs']

    for name in keywords_keep:
        if '-PTV' in name or name[0].lower()=='z':
            keywords_keep.remove(name)
    # print("Keeping these ")

    for i,patient in enumerate(config["patient_list"]):
        print("****************************************************")
        print(patient)
        plt.figure(i)
        patient = str(patient)
        if os.path.exists(PATH+patient):
            # try:
                run_anonymization(PATH,patient, save_path,keywords_keep,CT_name,produce_pdf)
            # except Exception as e:
            #     print("ERROR WITH PATIENT",patient,":",e)
            #     if produce_pdf:
            #         plt.subplots(figsize=(10,2))
            #         plt.title(patient+" ERROR")
            #         plt.text(0.5,0.5,e,dict(ha='center',va='center',fontsize=14,color='blue'))
                   
        else:   
            print("Patient directory "+ PATH+patient + " does not exist.")
    
    pdf = PdfPages("output.pdf")
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        pdf.savefig(fig,bbox_inches='tight')
    pdf.close()

    # if sys.argv[1:] == "all"
    # for patient in sys.argv[1:]:
    #     if patient == "all":
    #         list_patients_to_sort = sorted([f for f in os.listdir(PATH) if 'b' not in f and 'old' not in f],key=int)
        
        # Check if command line arguments correspond to existing patient directories
    #     elif os.path.exists(PATH+patient):
    #         list_patients_to_sort.append(patient)
    #     else:   
    #         print("Patient directory "+ PATH+patient + " does not exist.")
    

    # organize_multiple_patients(list_patients_to_sort, PATH)    

    # TO DO: fix issue where mutiple CT with diff names

    end = time.time()
    print("***TOTAL TIME***")
    print(end-start,"seconds")

