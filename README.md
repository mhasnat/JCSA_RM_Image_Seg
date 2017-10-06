# JCSA-RM RGBD Image Segmentation and Analysis Method
Repository for the MATLAB implementation of the "RGB-D image Segmentation using the Joint Color-Spatial-Axial clustering and Region Merging (JCSA-RM)" method.

- The JCSA-RM method is an RGB-D (joint color+depth) image segmentation method. This repo provides demos with/without a GUI with MATLAB code to perform the following tasks: <br>
**a.** load RGB-D image data from a mat file (contains RGB, Depth and Image Normals in an structure) and display them. <br>
**b.** Generate segmented image and display it. <br>

## How to use demo (tested in Matlab2017b):
**Run** the MATLAB file name: **RGBD\_Seg\_JCSA\_RM.m** for the GUI version and **demo\_NO\_GUI.m** otherwise<br>
**Load** data/samples files name: **rgbd\_info\_1.mat** and **rgbd\_info\_2.mat** <br>

## Application:
This segmentation method has been used to segment/analyze RGB-D images captured by the Microsoft Kinect camera. For details and other possible applications please see the references.

## Results to compare:
**JCSD\_RM\_Results.zip** file contains the results of applying JCSA-RM [1,2] method on the NYU depth database (NYUD2) [3]. Each result file consists of _segmentation_ - labels of pixels and _final scores_ – VoI, BDE, PRI and GTRC for the 1449 NYUD2 [3] images in half scaled (down) image.

## Code issues:
If you encounter error with - **computeTraceTerm** then go to the directory called 'rgbd' and compile mex file as:
**mex computeTraceTerm.cpp**.

## Extensions and scopes:
- You can extend the **JCSA** method for clustering heterogeneous data. However, for now the method is limited to the 3 Dimensional data with the directional distributions [4,5]. You can also extend it to work with higher dimensional data by extending [4] or [5].

- You can use the **RM** method independently to perform segmentation. It requires the clustering labels, the color image and image normals as input.

## References:

[1] Abul Hasnat, Olivier Alata and Alain Trmeau, Joint Color-Spatial-Directional clustering and Region Merging (JCSD-RM) for unsupervised RGB-D image segmentation, In IEEE Trans. on Pattern Analysis and Machine Intelligence (**TPAMI**), Vol 38, Issue 11, pages 2255 - 2268, **2015**.

[2] M Hasnat, O Alata, A Trémeau, Joint Color-Spatial-Directional clustering and Region Merging (JCSD-RM) for unsupervised RGB-D image segmentation, in arXiv preprint arXiv:1509.01788, **2015**.

[3] N. Silberman, D. Hoiem, P. Kohli, and R. Fergus, “Indoor segmentation and support inference from RGBD images,” in **ECCV** **2012**. Springer, 2012.

[4] Abul Hasnat, Olivier Alata and Alain Trmeau, Model-Based Hierarchical Clustering with Bregman Divergence and Fisher Mixture Model: Application to Depth Image Analysis, In Statistics and Computing (**STCO**), Vol 26, Issue 4, pp 861880, **2015**.

[5] Abul Hasnat, Olivier Alata and Alain Tremeau, Unsupervised Clustering of Depth Images using Watson Mixture Model, In Int. Conference on Pattern Recognition (**ICPR**) **2014**, Stockholm, Sweden.
