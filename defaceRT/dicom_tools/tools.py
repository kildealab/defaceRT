import numpy as np
from scipy import interpolate

debug = False
def xyz_to_image_coords(X,Y,Z,spacing,origin):
    X_new, Y_new, Z_new = [], [], []
    for x,y,z in zip(X,Y,Z):
        X_new.append((x-origin[0])/spacing[0])
        Y_new.append((y-origin[1])/spacing[1])
        Z_new.append((z-origin[2])/spacing[2])
    
    return X_new, Y_new, Z_new


def get_uni_spline(xi,yi):
    
    if len(xi) < 4:
         # Use a linear spline if there are fewer than 4 data points
        tck = interpolate.UnivariateSpline(xi, yi, k=1, s=0)
    else:
        # Fit spline with s=0 (passing through all points)
        tck = interpolate.UnivariateSpline(xi, yi, s=0)

        
#         # Evaluate spline for 1000 evenly spaced points
#         xj, yj = interpolate.splev(np.linspace(0, 1, 1000), tck)
    return tck
        
        
def reverse_pixels_hu(anon,scans):

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





def image_to_xyz_coords_single(X,spacing,origin,reverse_Z=False):
    if type(X) is not list:
        X = [X]
    X_new = []
    for x in X:
        if reverse_Z:
            X_new.append(-x*spacing+origin)
        else:
            X_new.append(x*spacing+origin)

    return X_new
    
def get_image_slice(start_z, z_smg, spacing,length=0,reverse_z=False):
    img_slice = int((abs(start_z - z_smg)/spacing[2]))
    if reverse_z:
        img_slice = length - img_slice - 1

    return img_slice



