#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 21:18:45 2025

@author: csun
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import e
import mrcfile
import os, tempfile

def get_3d_map_from_file(filename):
    if filename.endswith(".gz"):
        filename_final = filename[:-3]
        import gzip, shutil
        with gzip.open(filename, 'r') as f_in, open(filename_final, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        filename_final = filename
    with mrcfile.open(filename_final, mode="r") as mrc:
        data = mrc.data.copy()
        map_crs = [int(mrc.header.mapc), int(mrc.header.mapr), int(mrc.header.maps)]
        apix = mrc.voxel_size.x.item()
        return data, map_crs, apix

def get_3d_map_from_uploaded_file(fileobj):
    original_filename = fileobj.name
    suffix = os.path.splitext(original_filename)[-1]
    with tempfile.NamedTemporaryFile(suffix=suffix) as temp:
        temp.write(fileobj.read())
        temp.flush()
        return get_3d_map_from_file(temp.name)

def plot_fsc(x, fsc, y1, y3, method):
    fig, ax = plt.subplots()
    ax.plot(x, fsc, '-b', label="FSC")
    ax.plot(x, y3, '--r', label="Upper limit")
    ax.plot(x, y1, '--k', label="Lower limit")
    xmax = max(x)
    xstep = max(round(xmax / 5 / 2, 2) * 2, 0.01)  # nice rounding
    xticks = [xstep * i for i in range(int(xmax / xstep) + 1)]
    xlabels = [r"1/$\infty$"] + ["1/%.1f" % (1 / xt) for xt in xticks[1:]]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.axhline(y=0.143, linestyle='--', color='gray', linewidth=1, label="FSC = 0.143")
    ax.set_xlabel("Spatial Frequency (1/Å)")
    ax.set_ylabel("FSC")
    ax.set_title(f"FSC with Confidence Interval ({method})")
    ax.legend()
    st.pyplot(fig)

# --- Streamlit UI ---
st.title("FSC Error Bar Calculator")
st.markdown("""
### 🧊 About This FSC Confidence Estimation Web App

Fourier Shell Correlation (FSC) is the **standard method** used in the cryo-electron microscopy (cryo-EM) field to estimate the **resolution of reconstructed 3D maps**, serving as a key indicator of overall map quality. FSC is computed as the cross-correlation coefficient between two independently reconstructed "half maps" as a function of spatial frequency.

However, a single FSC value per shell assumes uniformity in signal quality across all voxels within that shell. In practice, this assumption oversimplifies the reality: **voxel intensities vary substantially across different spatial locations**, even within the same frequency band.

---

### Error Bar Estimation Methods

This web app allows you to **visualize confidence intervals** around the FSC values using **three different statistical approaches**. You can select your preferred method from the dropdown menu:

#### 1. **Fisher Z-Transform Method** (Default)
This approach assumes the Fisher z-transformed FSC values follow a normal distribution, allowing us to compute analytical confidence bounds.

**Formula**
""")

st.latex(r"""
\text{FSC}_{\text{upper/lower}} = \frac{e^{2(z \pm \sigma/\sqrt{n-3})} - 1}{e^{2(z \pm \sigma/\sqrt{n-3})} + 1}
""")

st.markdown("Where:")

st.latex(r"""
z = \frac{1}{2} \log \left( \frac{1 + \text{FSC}}{1 - \text{FSC}} \right)
""")

st.markdown("""
- *n* is the number of voxels in the shell  
- **σ** (sigma) is the user-defined confidence level (e.g. 3 for 99.7%)
""")

st.markdown("""
#### 2. **Bootstrap Method**  
This non-parametric approach resamples voxel pairs within each shell to build a distribution of FSC values. Confidence bounds are calculated from percentiles of this distribution.  

- You can adjust the number of bootstrap samples.  
- The bounds are based on the 2.5th and 97.5th percentiles by default (95% CI).  
""")

st.markdown("""
#### 3. **Variance-Based Method**  
This method estimates FSC uncertainty analytically using a published variance formula that accounts for the number of voxels and the observed FSC value. Confidence bounds are derived as:  
""")

st.latex(r"""
\text{FSC}_{\text{upper/lower}} = \text{FSC} \pm \sigma \cdot \sqrt{\text{Var(FSC)}}
""")

st.markdown("Where Var(FSC) is computed using:")

st.latex(r"""
\text{Var(FSC)} = \frac{(1 + 6 \cdot \text{FSC} + \text{FSC}^2)(1 - \text{FSC}^2)}{(1 + \text{FSC}^4) \cdot n}
""")

st.markdown("""
**Usage Notes** 

- Sigma (σ) defines your confidence level. For example:
  - σ = 1 → ~68% confidence  
  - σ = 2 → ~95% confidence  
  - σ = 3 → ~99.7% confidence  
""")

st.markdown("""
Note: If your half map files are larger than 200 MB, please download `calcFSC.py` from [my Github repository](https://github.com/sun647/FSCerrorbar/tree/main) and run it locally on your computer. There is no size limit when using the script locally.
""")

st.markdown("""
To learn more:

- [Penczek PA (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3165049/): *Resolution measures in molecular electron microscopy*.  
- [Cardone et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837392/): *One number does not fit all: mapping local variations in resolution*.
""")

st.markdown("Upload two MRC half-maps to compute FSC with confidence intervals.")
method = st.selectbox("Select error estimation method:", ["Fisher Z-transform", "Bootstrap", "Analytical Variance"])
sigma = st.number_input("Sigma (e.g., 3 for ~99.7% confidence)", min_value=0.1, max_value=10.0, value=3.0, step=0.1)

if method == "Bootstrap":
    bootstrap_reps = st.number_input("Number of bootstrap samples", min_value=10, max_value=1000, value=100, step=10)
else:
    bootstrap_reps = 0

file1 = st.file_uploader("Upload half map 1 (.mrc)", type="mrc")
file2 = st.file_uploader("Upload half map 2 (.mrc)", type="mrc")

run_button = st.button("Run")

if run_button and file1 and file2:
    try:
        imageData1, _, apix = get_3d_map_from_uploaded_file(file1)
        imageData2, _, _ = get_3d_map_from_uploaded_file(file2)
        #print (apix)
        if imageData1.shape != imageData2.shape:
            st.error("MRC volumes must have the same shape.")
        elif len(imageData1.shape) != 3:
            st.error("MRC volumes must be 3D.")
        else:
            model_1 = np.fft.fftshift(np.fft.fftn(imageData1))
            model_2 = np.fft.fftshift(np.fft.fftn(imageData2))

            box_half = model_1.shape[0] // 2
            fsc_list, y1, y3 = [], [], []
            box_size = model_1.shape[0]  # assume cubic volume
            x_vals=[]
            for i in range(1, box_half):
                a, b, c = model_1.shape[0]//2, model_1.shape[1]//2, model_1.shape[2]//2
                x, y, z = np.ogrid[-a:model_1.shape[0]-a, -b:model_1.shape[1]-b, -c:model_1.shape[2]-c]
                mask = (x**2 + y**2 + z**2 >= i**2) & (x**2 + y**2 + z**2 < (i+1)**2)

                v1 = model_1[mask]
                v2 = model_2[mask]

                if len(v1) < 4:
                    continue

                numerator = np.sum(v1 * np.conjugate(v2))
                denominator = np.sqrt(np.sum(np.abs(v1)**2) * np.sum(np.abs(v2)**2))
                fsc_val = (numerator / denominator).real
                fsc_list.append(fsc_val)
                n = len(v1)

                if method == "Fisher Z-transform":
                    zval = 0.5 * np.log((1 + fsc_val) / (1 - fsc_val))
                    s1 = zval - sigma / np.sqrt(n - 3)
                    s2 = zval + sigma / np.sqrt(n - 3)
                    lb = (e**(2*s1) - 1) / (e**(2*s1) + 1)
                    ub = (e**(2*s2) - 1) / (e**(2*s2) + 1)

                elif method == "Bootstrap":
                    fsc_samples = []
                    indices = np.arange(n)
                    for _ in range(bootstrap_reps):
                        idx = np.random.choice(indices, size=n, replace=True)
                        num = np.sum(v1[idx] * np.conjugate(v2[idx]))
                        den = np.sqrt(np.sum(np.abs(v1[idx])**2) * np.sum(np.abs(v2[idx])**2))
                        fsc_samples.append((num / den).real)
                    lb = np.percentile(fsc_samples, 2.5)
                    ub = np.percentile(fsc_samples, 97.5)

                elif method == "Analytical Variance":
                    var = ((1 + 6*fsc_val + fsc_val**2) * (1 - fsc_val**2)) / ((1 + fsc_val**4) * n)
                    lb = fsc_val - sigma * np.sqrt(var)
                    ub = fsc_val + sigma * np.sqrt(var)

                y1.append(lb)
                y3.append(min(ub,1))
                frequency = i / (box_size * apix)  # spatial frequency in 1/Å
                x_vals.append(frequency)
            plot_fsc(x_vals, fsc_list, y1, y3, method)

            df = pd.DataFrame({
                "Spatical Frequency (1/Å)": x_vals,
                "FSC": fsc_list,
                "Lower Limit": y1,
                "Upper Limit": y3
            })
            st.dataframe(df)

            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("Download FSC Data", csv, "fsc_errorbar.csv", "text/csv")

    except Exception as e:
        st.error(f"Error processing files: {e}")
