# JCSA-RM RGBD Image Segmentation and Analysis Method [1,2,3]
Repository for the source code of MATLAB implementation of the "RGB-D image Segmentation using the Joint Color-Spatial-Axial clustering and Region Merging (JCSA-RM)" method.

- The JCSA-RM method is an RGB-D (joint color+depth) image segmentation method. This repo provides demos with/without a GUI with MATLAB code to perform the following tasks: <br>
**a.** load RGB-D image data from a mat file (contains RGB, Depth and Image Normals in an structure) and display them. <br>
**b.** Generate segmented image and display it. <br>

## How to use demo (tested in Matlab2017b):
**Run** the MATLAB file name: **RGBD\_Seg\_JCSA\_RM.m** for the GUI version and **demo\_NO\_GUI.m** otherwise<br>
**Load** data/samples files name: **rgbd\_info\_1.mat**, **rgbd\_info\_2.mat**, **rgbd\_info\_1\_better\_normals.mat** and **rgbd\_info\_2\_better\_normals.mat**. <br>
- Select **\_better\_normals** in order to experiment with unambiguas surface normals.
- Choose different methods for testing, among: **(a)** JCSA, **(b)** JCSD, **(c)** JCSA-RM and **(d)** JCSD-RM.

## Application:
This segmentation method has been used to segment/analyze RGB-D images captured by the Microsoft Kinect camera. For details and other possible applications please see the references.

## Results to compare:
**JCSD\_RM\_Results.zip** file contains the results of applying JCSA-RM [1,2,3] method on the NYU depth database (NYUD2) [4]. Each result file consists of _segmentation_ - labels of pixels and _final scores_ – VoI, BDE, PRI and GTRC for the 1449 NYUD2 [4] images in half scaled (down) image.

## Code running issues:
It runs on Matlab2017b. If you encounter error with - **computeTraceTerm** then go to the directory called 'rgbd' and compile mex file as: <br>
**mex computeTraceTerm.cpp**.

## Extensions and scopes:
- You can extend the **JCSA** method for clustering heterogeneous data. However, for now the method is limited to the 3 Dimensional data with the directional distributions [5,6]. You can also extend it to work with higher dimensional data by extending [5] or [6].

- You can use the **RM** method independently to perform segmentation. It requires the clustering labels, the color image and image normals as input.

## References:

[1] Abul Hasnat, Olivier Alata and Alain Trmeau, Joint Color-Spatial-Directional clustering and Region Merging (JCSD-RM) for unsupervised RGB-D image segmentation, In IEEE Trans. on Pattern Analysis and Machine Intelligence (**TPAMI**), Vol 38, Issue 11, pages 2255 - 2268, **2015**. [`pdf download`](https://arxiv.org/pdf/1509.01788.pdf)

[2] M Hasnat, O Alata, A Trémeau, Joint Color-Spatial-Directional clustering and Region Merging (JCSD-RM) for unsupervised RGB-D image segmentation, in arXiv preprint arXiv:1509.01788, **2015**. [`pdf download`](https://arxiv.org/pdf/1509.01788.pdf)

[3] M. A. Hasnat, O. Alata, and A. Tremeau, “Unsupervised RGB-D image segmentation using joint clustering and region merging,” in British Machine Vision Conference (BMVC). BMVA Press, **2014**. [`pdf download`](http://www.bmva.org/bmvc/2014/papers/paper045/index.html)

[4] N. Silberman, D. Hoiem, P. Kohli, and R. Fergus, “Indoor segmentation and support inference from RGBD images,” in **ECCV** **2012**. Springer, 2012.

[5] Abul Hasnat, Olivier Alata and Alain Trmeau, Model-Based Hierarchical Clustering with Bregman Divergence and Fisher Mixture Model: Application to Depth Image Analysis, In Statistics and Computing (**STCO**), Vol 26, Issue 4, pp 861880, **2015**. [`pdf download`](https://link.springer.com/epdf/10.1007/s11222-015-9576-3?author_access_token=guKHP_biEt1GZJH3orEjaPe4RwlQNchNByi7wbcMAY7phLRhrbt4ZOG8_2kzHa9Hk7sydZ5efCB8Saw_dieaoRCLdI5FVvwJwcLaa_D3M8GTfE1Hrg5TgMmpQisxS7GVaHibDAPx7v9GUyoqrO5GBw%3D%3D)

[6] Abul Hasnat, Olivier Alata and Alain Tremeau, Unsupervised Clustering of Depth Images using Watson Mixture Model, In Int. Conference on Pattern Recognition (**ICPR**) **2014**, Stockholm, Sweden. [`pdf download`](https://www.researchgate.net/publication/264505829_Unsupervised_Clustering_of_Depth_Images_Using_Watson_Mixture_Model)
