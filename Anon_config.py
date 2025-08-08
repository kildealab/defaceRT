config ={
    "path": '../', # Path to patient directories 
    "patient_list": ['Example_CTs'], # Patient directory names in path
	"save_path":'results/', # Path where new dicoms will be saved
	"contours_to_keep": ["PTV","brain"], # Contour keywords to be kept in defaced image
	"eye_contour_keywords": ["eye","globe","orbit"], # Eye contour keywords 
	"body_contour_keywords": ['body','external'], # Body contour keywords
	"CT_name": '', # leave blank if you want it to find the CT
	'CT_keyword':'', # A keyword present in all CT names from the TPS, helps to differentiate between CBCT names
	'CT_dir_name_min_length':0, # minimum string length of CT name
	'CT_dir_name_max_length':100, # note maximum string length allowed in Structure Set Label tag = 16
	'ignore_keywords_in_CT':['copy', 'do not use'], # Words to ignore in CT names (ie not real planning CTs)
	"print_PDFs":True

}
