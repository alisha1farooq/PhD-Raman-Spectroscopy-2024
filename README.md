# PhD-Raman-Spectroscopy-2024
MATLAB and Python scripts for single-cell analysis

This repository contains code used during my PhD research at the University of Leeds for the preprocessing and classification of single-cell Raman spectral data.

The goal of the project was to distinguish between healthy, pre-cancerous, and cancerous states in oesophageal cell lines using label-free Raman spectroscopy combined with multivariate analysis.

## 📂 Repository Structure

### 🔧 MATLAB Scripts – *Preprocessing Pipeline*
Located in the `MATLAB/` folder:

- `script1.m`: Sets up file paths, dataset indices, and calibration parameters for batch processing `.wdf` Raman files.
- `script2.m`: Performs preprocessing on raw spectra, including:
  - Calibration using reference peaks
  - Spline interpolation
  - Smoothing and baseline correction
  - Normalisation
  - Output of cleaned spectra and summary statistics

> These scripts prepared the data for use in Python-based classification.

---

### 🧠 Python Scripts – *PCA/LDA Classification*
Located in the `Python/` folder:

`script3.py`:  
  Performs baseline correction (Zhang, ALS, SNV) and preprocessing steps such as interpolation, standardisation, and spectral region selection.

- `script4.py`:  
  Applies PCA (Principal Component Analysis) and LDA (Linear Discriminant Analysis) for visualisation and classification of processed Raman spectra. Includes cross-validation and accuracy scoring.

- `script5.py`:  
  Additional utility functions and extensions used in the classification or preprocessing pipeline (e.g. plotting, outlier handling, advanced parameter tuning).

> These scripts were collaboratively developed. My contribution involved running and adapting the analysis for my datasets, finetuning preprocessing parameters, and interpreting PCA/LDA outputs.

---

## 📌 Notes

- This repository is intended for transparency, reproducibility, and to demonstrate my applied experience in computational biology and Raman data analysis.
- This work is linked to an open-access repository at the University of Leeds:
  🔗 [https://doi.org/10.5518/1246]

---

## 🧪 Techniques Used

- MATLAB: `.wdf` file handling, interpolation, EMSC/baseline correction, smoothing, normalisation
- Python: `numpy`, `scikit-learn`, `matplotlib`, `PyQt5`, PCA/LDA modelling, 3D plotting

