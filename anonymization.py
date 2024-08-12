
import os, time, sys
from datetime import datetime
import pydicom as dcm
import numpy as np
from scipy import interpolate

sys.path.append('../')
from Slice_Selection.slice_selection import *

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

# to do, fix z mess for when its revesrsed
def generate_anon_image(image,z_lists,spacing, start_z,y_cutoff, cut_above = True,isodose=[], contours_to_keep = {},origin=[0,0,0],reverse_z=False):
    anon_image = image.copy()

    all_z_slices = sorted(z_lists[0] + [z for z in z_lists[1] if z not in z_lists[0]])
    img_slices = [get_image_slice(start_z, z, spacing) for z in all_z_slices]
    
    
    z_min = min(img_slices)
    z_max = max(img_slices)
    
    if cut_above:
        if reverse_z:
            z_min = 0
        else:
            z_max = len(image)

    
    mask = np.full_like(image,1)
    mask[z_min:z_max+1,0:int(y_cutoff)+1,:] = 0

    for contour in contours_to_keep:
        for c_slice in contours_to_keep[contour]:
            slice_num = get_image_slice(start_z, c_slice[2], spacing)
            if slice_num >= z_min and slice_num <= z_max:
                xi, yi, zi = c_slice[::3], c_slice[1::3], c_slice[2::3]
                X,Y,Z = xyz_to_image_coords(xi,yi,zi,spacing,origin)

                cont_xy = np.array([tuple([int(x),int(y)]) for x,y in zip(X,Y)])
                cont_xy = cont_xy.reshape((-1,1,2))
                cv2.drawContours(mask[slice_num],[cont_xy],-1,1, -1)

  # TO DO: ISODOSE BELOW      
    # print("anon*MASK************")    
    # print((anon_image*mask)[0])
    # print(mask)
    # print(mask*-1+1)
    # print(anon_image)
    return anon_image*mask++(mask*-1+1)*(-1000)


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
                print("xy",xy)
                xs, ys = zip(*xy)
                print("xs",xs)
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
    all_z_slices = sorted(z_lists[0] + [z for z in z_lists[1] if z not in z_lists[0]])
    print(body_name)
#     x_iso = X_new
#     y_iso = Y_new
#     z_iso = Z_new

#TODO:add isodoes
    contour_curves = {}
  
    for contour in contours_to_keep:
        # print(contour)
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

    return full_stack_N

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
    image_path = os.path.join(PATH,patient,CT_file)
    # CT = dcm.read_file(image_path)
    RS_file = [f for f in os.listdir(image_path) if f[0:2] == 'RS'][0]
    RS = dcm.read_file(os.path.join(image_path,RS_file))

    #TO DO:load dose and RP fiels
    return RS,RS_file

def init_data(patient_path, CT_file, RS):

    slices = [dcm.read_file(PATH+patient+'/'+CT_file + '/'+ s) for s in os.listdir(PATH+patient+'/'+CT_file) if s[0:2] =='CT']
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



def get_first_CT(patient_path):
    
    CT_list = [d for d in os.listdir(patient_path) if d[9:11] == 'CT' and len(d) == 23]
    
    CT_list.sort(key=lambda x: datetime.strptime(x[12:], "%d_%b_%Y"))
    print(CT_list[0])
    return CT_list[0]

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
    
def save_RT_struct(RS, RS_file,contour_stack,save_path,patient,CT_file):
    #TODO: make it choose roi name, putting body for now
    for i, seq in enumerate(RS.StructureSetROISequence):
        if seq.ROIName == 'BODY':
                index = i
                break
    
    dict_stack = {}
    for row in contour_stack:
        dict_stack[row[2]] = row

    for i, ROI_contour_seq in enumerate(RS.ROIContourSequence[index].ContourSequence):
        z_ref = ROI_contour_seq.ContourData[2] 
        if z_ref not in list(dict_stack.keys()):
            print("******",z_ref,"not in dict")

        print(z_ref,dict_stack(z_ref))

        RS.ROIContourSequence[index].ContourSequence[i].ContourData = dict_stack[z_ref]

    RS.save_as(os.path.join(save_path,patient,CT_file,RS_file))
        # contour_coords.append(ROI_contour_seq.ContourData)


def run_anonymization(PATH,patient,save_path,keywords_keep = [],CT_name=''):
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

    y_cutoff,z_lists,y_cutoff_roi,z_smg =  get_eye_contours(RS,start_x,start_y,start_z,z_spacing,pixel_spacing)


    if len(keywords_keep) == 0:
        dict_contours_keep = {}
    else:
        list_names_keep = []
        for keyword in keywords_keep:
            list_names_keep = list_names_keep + find_ROI_names(RS,keyword=keyword)


        dict_contours_keep,_ = get_all_ROI_contours(list_names_keep,RS)


    anon = generate_anon_image(image,z_lists,spacing,start_z,y_cutoff,reverse_z=reverse_z,contours_to_keep=dict_contours_keep,origin=origin)
    # print("anon")
    # print(anon[0])
    new_dicom = anon_dicom(anon, slices)
    save_dicom(new_dicom,save_path,patient,CT_file)

    body_names = find_ROI_names(RS,'body')
    print(body_names)
    dict_contours_body,_ = get_all_ROI_contours([body_names[0]],RS)
    # print(len(dict_contours_body['BODY']))
    # print("******************")

    full_stack_N = generate_anon_body(dict_contours_body,z_lists=z_lists,y_cutoff=y_cutoff_roi,reverse_z=reverse_z,contours_to_keep=dict_contours_keep)
    # for slice in full_stack_N:
    #     if slice[2] == z_smg:
    # #         print(slice)
    #         break
    # x_new_body = slice[::3] 
    # y_new_body = slice[1::3] 
    # z_new_body =slice[2::3] 
    # x_new, y_new, z_new = xyz_to_image_coords(x_new_body,y_new_body,z_new_body,spacing,origin)

    # print(len(full_stack_N))
    save_RT_struct(RS, RS_file,full_stack_N,save_path,patient,CT_file)
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


if __name__ == "__main__":
    start = time.time()
    
    from Anon_config import *

    PATH = config["path"]
    keywords_keep = config["contours_to_keep"]
    save_path = config["save_path"]
    CT_name = config["CT_name"]



    for patient in config["patient_list"]:
        patient = str(patient)
        if os.path.exists(PATH+patient):
            run_anonymization(PATH,patient, save_path,keywords_keep,CT_name)
        else:   
            print("Patient directory "+ PATH+patient + " does not exist.")
    


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

