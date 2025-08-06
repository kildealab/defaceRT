import numpy as np
import matplotlib.pyplot as plt
from tools import find_ROI_names, get_ROI_colour_dict, image_to_xyz_coords_single, get_all_ROI_contours, get_ROI_slice,get_ROI_pixel_array
from collections import OrderedDict

plt.rc('font', size=8)          # controls default text sizes
plt.rc('axes', titlesize=12)     # fontsize of the axes title
plt.rc('axes', labelsize=4)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=4)    # fontsize of the tick labels
plt.rc('ytick', labelsize=4)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    # legend fontsize
plt.rc('figure', titlesize=6)  # fontsize of the figure title

def plot_all_contours(RS,image,slice_num,origin,spacing,reverse_z=False,ignore_terms=[],legend=False):
    contours_not_on_slice = []
   
    all_ROIs = [r for r in find_ROI_names(RS)]
    for term in ignore_terms:
        all_ROIs = [r for r in all_ROIs if term.lower() not in r.lower()]
    
#     dict_contours, z_lists_b = get_all_ROI_contours(all_ROIs, RS)
    colours= get_ROI_colour_dict(RS)
    # print(colours)
    z_slice = image_to_xyz_coords_single(slice_num, spacing[2],origin[2],reverse_z)[0]
    print("Z slice:",z_slice,"Origin:",origin[2],"Spacing:",spacing[2],"Slice number:",slice_num)
#     xyz_to_image_coords_single(X,spacing,origin):
    print(z_slice)
    # 
    
    plt.imshow(image[slice_num],cmap='gray')
    for roi in all_ROIs:
        print(roi)
        dict_contours, z_lists = get_all_ROI_contours([roi], RS)
        # print("X")
#         if len(dict_contours) >1:
        for i,r in enumerate(dict_contours):
            if r == roi:
                break
        try:       
            print(z_slice, z_lists[i])
            roi_slice = get_ROI_slice(z_slice,z_lists[i])
            print("ROI",len(roi_slice))
            
            c =colours[roi]
            for s in roi_slice:
                # print(s)
                # print("**")
                roi_x, roi_y = get_ROI_pixel_array(dict_contours[roi][s],origin[0],origin[1],spacing)
                plt.plot(roi_x,roi_y,'-',color=  (c[0]/255,c[1]/255,c[2]/255),label=roi)
        except Exception as e:
            contours_not_on_slice.append(roi)
            # print(roi,"not on slice.")
        
    if legend:
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))

        plt.legend(by_label.values(), by_label.keys(),prop={'size': 6},bbox_to_anchor = (1,0.5),loc='center left')
    # print("Other contours not on slice:",contours_not_on_slice)
    
    


def plot_3_views(slices, slice_number=100,image=[],plot_cor=False,patient=''):
    ps = slices[0].PixelSpacing
    ss = slices[0].SliceThickness 
    ax_aspect = ps[1]/ps[0]
    sag_aspect = ps[1]/ss
    cor_aspect = ss/ps[0]

    # create 3D array
    img_shape = list(slices[0].pixel_array.shape)
    img_shape.append(len(slices))
    img3d = np.zeros(img_shape)
    
    # fill 3D array with the images from the files
    for i, s in enumerate(slices):
        if len(image)==0:
            img2d = s.pixel_array
        else:
            img2d = image[i]
        img3d[:, :, i] = img2d

    # plot 3 orthogonal slices
    a1 = plt.subplot(2, 2, 1)
    plt.imshow(img3d[:, :, slice_number],cmap='gray')
    a1.set_aspect(ax_aspect)

    plt.title(patient)

    a2 = plt.subplot(2, 2, 2)
    plt.imshow(img3d[:, img_shape[1]//2, :],cmap='gray')
    a2.set_aspect(sag_aspect)
    # if plot_cor:
    #     a3 = plt.subplot(2, 2, 3)
    #     plt.imshow(img3d[img_shape[0]//2, :, :].T)
    #     a3.set_aspect(cor_aspect)

    # figure_list.append(plt.figure())
    # plt.show()