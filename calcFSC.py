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

def plot_fsc(x, fsc, y1, y3):
    fig, ax = plt.subplots()
    ax.plot(x, fsc, '-b', label="FSC")
    ax.plot(x, y3, '--r', label="Upper 3σ")
    ax.plot(x, y1, '--k', label="Lower 3σ")
    ax.set_xlabel("Shell Index")
    ax.set_ylabel("FSC")
    ax.set_title("Fourier Shell Correlation with Confidence Interval")
    ax.legend()
    st.pyplot(fig)

# --- Streamlit UI ---
st.title("FSC Error Bar Calculator")
st.markdown("""
### 🧊 About This FSC Confidence Estimation Web App

Fourier Shell Correlation (FSC) is the **standard method** used in the cryo-electron microscopy (cryo-EM) field to estimate the **resolution of reconstructed 3D maps**, serving as a key indicator of overall map quality. FSC is computed as the cross-correlation coefficient between two independently reconstructed "half maps" as a function of spatial frequency.

However, a single FSC value per shell assumes uniformity in signal quality across all voxels within that shell. In practice, this assumption oversimplifies the reality: **voxel intensities vary substantially across different spatial locations**, even within the same frequency band.

This web app computes FSC **not just as an average**, but includes **statistical confidence intervals** using a user-defined sigma. It uses the Fisher z-transform to calculate bounds:

#### FSC confidence interval formula:
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
Note: If your half map files are larger than 200 MB, please download `calcFSC.py` from [my Github repository](https://github.com/sun647/FSCerrorbar/tree/main) and run it locally on your computer. There is no size limit when using the script locally.
""")
st.markdown("""
To learn more:

- [Penczek PA (2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3165049/): *Resolution measures in molecular electron microscopy*.  
- [Cardone et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837392/): *One number does not fit all: mapping local variations in resolution*.
""")

sigma = st.number_input("Sigma for confidence interval (default value is 3)", min_value=0.1, max_value=10.0, value=3.0, step=0.1)

file1 = st.file_uploader("Upload half map1.mrc", type="mrc")
file2 = st.file_uploader("Upload half map2.mrc", type="mrc")

if file1 and file2:
    try:
        imageData1, _, apix = get_3d_map_from_uploaded_file(file1)
        imageData2, _, _ = get_3d_map_from_uploaded_file(file2)

        if imageData1.shape != imageData2.shape:
            st.error("Error: MRC volumes must be the same shape.")
        elif len(imageData1.shape) != 3:
            st.error("Error: Input volumes must be 3D.")
        else:
            model_1 = np.fft.fftshift(np.fft.fftn(imageData1))
            model_2 = np.fft.fftshift(np.fft.fftn(imageData2))

            box_half = model_1.shape[0] // 2
            fsc_list, s1, s2 = [], [], []
            y1, y3 = [], []

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
                fsc = numerator / denominator
                fsc_list.append(fsc.real)

                n = len(v1)
                zval = 0.5 * np.log((1 + fsc.real) / (1 - fsc.real))
                s1.append(zval - sigma / np.sqrt(n - 3))
                s2.append(zval + sigma / np.sqrt(n - 3))

            for i in range(len(s1)):
                lb = (e**(2*s1[i]) - 1) / (e**(2*s1[i]) + 1)
                ub = (e**(2*s2[i]) - 1) / (e**(2*s2[i]) + 1)
                y1.append(lb)
                y3.append(ub)

            x_vals = list(range(len(fsc_list)))
            plot_fsc(x_vals, fsc_list, y1, y3)

            df = pd.DataFrame({
                "Shell Index": x_vals,
                "FSC": fsc_list,
                "Lower 3σ": y1,
                "Upper 3σ": y3
            })
            st.dataframe(df)

            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("Download FSC Data", csv, "fsc_errorbar.csv", "text/csv")

    except Exception as e:
        st.error(f"Error processing files: {e}")
