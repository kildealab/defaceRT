# defaceRT

An automated defacing algorithm specifically designed for head and neck cancer (HNC) CT scans and associated DICOM-RT data (Structure Sets and Dose maps). Our algorithm preserves critical radiotherapy structures while removing identifiable facial features. Requires a DICOM-RT Structure Set file with the eyes contoured.


## Table of Contents
- [Motivation](#Motivation)
- [Features](#Features)
- [Dependencies](#Dependencies)
- [Installation](#Installation)
- [Usage](#Usage)
- [Contributing](#Contributing)
- [Contact](#Contact)
- [Disclaimer](#Disclaimers)
- 

## Motivation

The proliferation of medical imaging datasets in the public domain has raised concerns about potential patient reidentification from CT scans of the head. Various defacing algorithms have been produced to mitigate such concerns. However, these tools overlook the need to preserve critical structures that are important for radiotherapy research. Therefore, we developed a novel automated defacing algorithm that preserves Organs at Risk (OARs) and Planning Target Volumes (PTVs) for radiotherapy planning while removing identifiable features from head and neck cancer (HNC) CT scans and associated DICOM-RT data.

#### Overview of how the defacing algorithm rowkr
<p align="center">
  <img width="500"  alt="Figure 1 - Workflow" src="https://github.com/user-attachments/assets/8e1cc9ce-d1c0-4bc9-88a7-27794454a043" />
</p>

## Features

- **Automated processing**: Uses eye contours as landmarks to deface automatically.
- **Structure preservation**: Maintains all PTV pixels even when they overlap defaced regions .
- **Multi-modal defacing**: Handles CT images, Structure Sets, and Dose maps
- **Visualization**: Generates PDF reports for quality assurance

## Requirements

### Dependencies
```
python >= 3.8
pydicom
numpy
scipy
matplotlib
opencv-python (cv2)
```

### Input Data Requirements
- DICOM CT scans
- DICOM-RT Structure Set files (with eye contours)
- DICOM-RT Dose files (optional)
- Eye contours must be present in Structure Set

## Installation

1. **Clone the main repository:**
```bash
git clone https://github.com/kildealab/defaceRT.git
cd defaceRT
```

2. **Install dependencies:**
```bash
pip install -r requirements.txt
```


## Usage
1. **Organize data** according to the following directory structure:
```
/path/to/your/data/
├── 📁 ppatient1/
│   ├── 📁 CT_folder_1/          # Auto-detected or specify in config
│   │   ├── 📄CT.*.dcm           # CT scan DICOM files
│   │   ├── 📄RS.*.dcm           # RT Structure Set (must contain eye contours)
│   │   └── 📄RD.*.dcm           # RT Dose (optional)
│   └── 📁 CT_folder_2/        
├── 📁 patient2/
│   └── ...
```

2. **Configure settings** in `Anon_config.py`:
    Please see comments in `Anon_config.py` for an explanation of the other optional configuration variables.
```python
config = {
    "path": "/path/to/your/data/",                # Root directory containing patient folders
    "save_path": "/path/to/output/",              # Path to save defaced DICOMs
    "patient_list": ["patient1", "patient2"],     # List of patient IDs to process
    "contours_to_keep": ["PTV", "CTV", "Brain"],  # ROI names to preserve in image even if within cropped bounds
`   "rint_PDFs": True                             # Generate visualization PDFs for quality assurance
  ...
    }
```


3. **Run defacing**:
```bash
python anonymization.py
```

4. **Quality assurance**: check `output.pdf` to ensure proper defacing.

*Note*: Example test CTs coming soon.
### Output
* DICOM Files: Defaced DICOM CT, RT Structure Set, and Dose files saved to `save_path` with the same initial patient directory hierarchy.
* `output.pdf`: Visual verification of defacing results.




## Testing and Validation
For researchers wanting to reproduce our validation study or validate on new datasets as per PAPER COMING SOON, please see validation tests and instructions at [github.com/kildealab/defaceRT-validation](https://github.com/kildealab/defaceRT-validation.git).

The validation repository includes:
- **Facial recognition testing**: FaceNet512-based privacy validation
- **Auto-segmentation evaluation**: Structure preservation assessment
- **PTV Location Analysis**: Target preservation assessment  

*Note: Validation requires proprietary software (LimbusAI) and manual data preparation steps not needed for basic defacing.*



## Citation

If you use DefaceRT in your research, please cite: PAPER COMING SOON


## Contributing

We welcome contributions! If you are interested in contributing, please fork the repository and create a pull request with your changes.
1. Fork the repository.
2. Create a new branch: `git checkout -b feature-name`
3. Make your changes and commit them.
4. Push your branch: `git push origin feature-name`
5. Create a pull request.

## Contact

For support or questions, please email Kayla O'Sullivan-Steben at kayla.osullivan-steben@mail.mcgill.ca.

---
## Disclaimers
* I AM IN THE PROCESS OF CLEANING THIS CODE, I KNOW IT'S A MESS!!!!
* This software is for research purposes only. Always comply with institutional policies and regulations regarding medical data sharing.
* There seems to be an issue with defacing dose maps in situations where the CT image z-values increase in the opposite direction as the dose maps'. 
* This is not the most beautifully packaged or polished code, but I hope it still proves useful. I have also only tested this code on our in-house data. That being said, I happily welcome suggestions, improvements, and contributions!
