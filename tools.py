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
    image_path = os.path.join(PATH,patient,CT_file)
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
