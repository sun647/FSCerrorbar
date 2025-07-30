
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import mrcz

st.title("FSC Error Bar Calculator")

st.markdown("Upload two MRC maps to compute FSC with error bars.")

file1 = st.file_uploader("Upload half map1.mrc", type="mrc")
file2 = st.file_uploader("Upload half map2.mrc", type="mrc")

if file1 and file2:
    imageData1, _ = mrcz.readMRC(file1)
    imageData2, _ = mrcz.readMRC(file2)

    if imageData1.shape != imageData2.shape:
        st.error("Error: MRC volumes must be the same shape.")
    elif len(imageData1.shape) != 3:
        st.error("Error: Input volumes must be 3D.")
    else:
        model_1 = np.fft.fftshift(np.fft.fftn(imageData1))
        model_2 = np.fft.fftshift(np.fft.fftn(imageData2))

        fsc_list = []
        err_list = []
        N_shell = []

        fr = model_1.shape[0] // 2

        for i in range(1, fr):
            a, b, c = model_1.shape[0] // 2, model_1.shape[1] // 2, model_1.shape[2] // 2
            x, y, z = np.ogrid[-a:model_1.shape[0]-a, -b:model_1.shape[1]-b, -c:model_1.shape[2]-c]
            mask = (x**2 + y**2 + z**2 >= i**2) & (x**2 + y**2 + z**2 < (i+1)**2)

            A = model_1[mask]
            B = model_2[mask]
            N = len(A)

            if N < 10:
                continue

            AB = np.real(np.sum(A * np.conjugate(B)))
            AA = np.real(np.sum(A * np.conjugate(A)))
            BB = np.real(np.sum(B * np.conjugate(B)))

            fsc = AB / np.sqrt(AA * BB)
            fsc_list.append(fsc)
            N_shell.append(N)

            err = np.sqrt((1 - fsc**2)**2 / (N - 1))
            err_list.append(err)

        fig, ax = plt.subplots()
        ax.errorbar(range(1, 1 + len(fsc_list)), fsc_list, yerr=err_list, fmt='o')
        ax.set_xlabel("Shell Index")
        ax.set_ylabel("FSC")
        ax.set_title("FSC Curve with Error Bars")
        ax.grid(True)
        st.pyplot(fig)
