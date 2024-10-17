import pydicom as dcm
import os

path = '/data/kayla/anon_images/patients/test_id/'

for f in os.listdir(path):
    if f.endswith(".dcm"):
        file_path = os.path.join(path, f)
        ds = dcm.dcmread(file_path)

        ds.PatientID = "anon_test"
        
        ds.save_as(path+f)
