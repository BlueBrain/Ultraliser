# Introduction

_Ultraliser_ is an __unconditionally robust__ and __high performance__ framework dedicated primarily to **_in silico_ neuroscience** research. 
_Ultraliser_ is capable of generating high fidelity and multiscale 3D models (surface meshes and annotated volumes) of neuroscientific data, such as: nuclei, mitochondria, endoplasmic reticula, 
neurons, astrocytes, pericytes, neuronal branches with dendritic spines, minicolumns with thousands of neurons and large networks of cerebral vasculature - with realistic geometries. 

_Ultraliser_ implements an effective __voxelization-based remeshing engine__ that can rasterize non-watertight surface meshes - in the 
form of __triangular soups__ - into high resolution volumes, with which we can reconstruct topologically accurate, adaptively optimized 
and watertight surface manifolds. 

In addition to their importance for accurate quantitative analysis, the resulting models are primarily intended to automate the 
process of conducting supercomputer-based _in silico_ simulations of neuroscience experiments; complementing _in vivo_ and _in vitro_ 
techniques. 

Watertight triangular meshes are used for (i) performing 3D particle simulations, (ii) mesh-based skeletonization, in which accurate morphologies of cellular structures are obtained for performing 1D compartmental simulations and (iii) tetrahedralization, in which we can generate tetrahedral volume meshes for 3D reaction-diffusion simulations.
Annotated volumetric tissue models are also used in _in silico_ imaging studies, where we can simulate optical imaging experiments with brightfield or fluorescence microscopy<sup>10</sup>. 

## Features 

+ Reconstruction of high fidelity, optimized<sup>1</sup> and two-manifold watertight<sup>2</sup> triangular mesh models from non-watertight inputs represented by polygonal soups. 
+ Surface mesh smoothing and optimization using [Laplacian operators](https://en.wikipedia.org/wiki/Laplacian_smoothing) and feature-preserving adaptive mesh optimization<sup>1</sup>.
+ Reconstruction of large-scale volumetric models<sup>3</sup> from non-watertight input meshes using high performance [surface and solid voxelization](Mesh-Voxelization). 
+ Reconstruction of optimized and smooth surface meshes from input volumes using parallel implementations of the standard marching cubes<sup>4</sup> algorithm and the advanced dual marching cubes<sup>5</sup> algorithm. 
+ Reconstruction of optimized and smooth surface meshes from input binary masks of segmented data. 
+ Reconstruction of geometrically realistic watertight mesh models of spiny neurons from corresponding morphological skeletons<sup>6</sup> .
+ Reconstruction of geometrically realistic watertight mesh models of complete astroglial cells<sup>8</sup>  (with endfeet) from input morphological skeletons and endfeet surface patches<sup>9</sup> .        
+ Reconstruction of high-fidelity, optimized and multi-partitioned vascular meshes from fragmented and large-scale vascular network graphs<sup>7</sup>. 
+ Morphology, mesh and volume quantitative and qualitative analysis. 
+ Generation of color-coded multi-axis projections of spatial data (morphologies, meshes and volumes) for visual analytics. 

## Documentation 

Exhaustive user documentation, including step-by-step examples and detailed explanation of the 
command line options, is available on the [Wiki](https://github.com/BlueBrain/Ultraliser/wiki) of this repository.

## Installation

Installation instructions are detailed in this [page on the Wiki](https://github.com/BlueBrain/Ultraliser/wiki/Installation). 

### Software Dependencies 

* [OpenMP](https://en.wikipedia.org/wiki/OpenMP), a multi-threading library for parallel processing on multi-core CPUs. 
* [libTIFF](http://www.libtiff.org/), which gives support for the Tag Image File Format (TIFF), a widely used format for storing image data.
* [libhdf5](https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html), or the Hierarchical Data Format 5 (HDF5) library for storing data.
* [Eigen3](https://eigen.tuxfamily.org/), a template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [BZip2](https://www.sourceware.org/bzip2/), a high quality data compressor.  
* [ZLIB](https://docs.python.org/3/library/zlib.html), for data compression. 
* [FMT](https://github.com/fmtlib/fmt), a formatting library providing a fast and safe alternative to C stdio and C++ iostreams.
* [GLM](https://github.com/g-truc/glm), a header only C++ mathematics library for graphics software based on the OpenGL Shading Language (GLSL) specifications.

### Supported Operating Systems 

_Ultraliser_ has been tested on Unix-based operating systems including:

* Ubuntu 18.04, Ubuntu 20.04, Ubuntu 21.04 and Ubuntu 22.04.
* RHEL7, RHEL8. 
* macOS 10.12 Sierra, 10.13 High Sierra,  10.14 Mojave, 10.15 Catalina.


# Known Bugs or Feature Requests

Please refer to the [github issue tracker](https://github.com/BlueBrain/Ultraliser/issues?utf8=%E2%9C%93&q=) for fixed and open bugs. User can also report any bugs and request new features needed for their research. We are happy to provide direct [support](#contact) . 

# License 

_Ultraliser_ is available to download and use under the GNU General Public License, version 3 ([GPL](https://www.gnu.org/licenses/gpl.html), or “free software”). 
The code is open sourced with approval from the open sourcing committee and principal coordinators of the Blue Brain Project in March 2021. 
See file [LICENSE](https://github.com/BlueBrain/Ultraliser/blob/master/LICENSE) for the full license.

# Citation

If you use _Ultraliser_, you can use the following ${\mathrm{B{\scriptstyle{IB}} T_{\displaystyle E}X}}$ entry for citation

```
@article {Abdellah2022.07.27.501675,
    author = {Abdellah, Marwan and Garc{\'\i}a Cantero, Juan Jos{\'e} and Roman Guerrero, Nadir 
    and Foni, Alessandro and Coggan, Jay S. and Cal{\`\i}, Corrado and Agus, Marco and 
    Zisis, Eleftherios and Keller, Daniel and Hadwiger, Markus and Magistretti, Pierre and 
    Markram, Henry and Sch{\"u}rmann, Felix},
    title = {Ultraliser: a framework for creating multiscale, high-fidelity and geometrically 
    realistic 3D models for in silico neuroscience},
    elocation-id = {2022.07.27.501675},
    year = {2022},
    doi = {10.1101/2022.07.27.501675},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2022/07/29/2022.07.27.501675},
    journal = {bioRxiv}
}
```

# Publications 

The volume reconstruction algorithms in _Ultraliser_ are based on the following paper. 

```
@article{abdellah2017reconstruction,
  title={Reconstruction and visualization of large-scale volumetric models of neocortical 
  circuits for physically-plausible in silico optical studies},
  author={Abdellah, Marwan and Hernando, Juan and Antille, Nicolas and Eilemann, Stefan and 
  Markram, Henry and Sch{\"u}rmann, Felix},
  journal={BMC bioinformatics},
  volume={18},
  number={10},
  pages={402},
  year={2017},
  publisher={BioMed Central}
}
```

# Acknowledgement & Funding

The development of this software was supported by funding to the [Blue Brain Project](https://bluebrain.epfl.ch), a research center of the [École polytechnique fédérale de Lausanne (EPFL)](https://www.epfl.ch/), from the Swiss government’s ETH Board of the Swiss Federal Institutes of Technology.
Financial support was provided by competitive research funding from [King Abdullah University of Science and Technology (KAUST)](https://www.kaust.edu.sa/en).

# Attributions

* The volume reconstruction code is an extension to the work of [Marwan Abdellah's](http://marwan-abdellah.com/) [PhD (In silico Brain Imaging: Physically-plausible Methods for Visualizing Neocortical Microcircuitry)](https://infoscience.epfl.ch/record/232444?ln=en). 

* The mesh optimization code in _Ultraliser_ is based on the routiens provided by the [GAMer (Geometry-preserving Adaptive MeshER) library](http://fetk.org/codes/gamer/). GAMer is developed by _Z. Yu, M. Holst, Y. Cheng, and J.A. McCammon_, and can be redistributed and modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.

* The watertighness verification code in _Ultraliser_ is based on an extended version of the [MeshFix library](https://github.com/MarcoAttene/MeshFix-V2.1). MeshFix is developed by _Marco Attene_, Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche, Sezione di Genova, IMATI-GE / CNR, and can be redistributed and modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or any later version.    

* The mesh analysis code is implemented based on the metrics described in [The Verdict Geometric Quality Library](https://coreform.com/papers/verdict_quality_library.pdf). 

* The morphology analysis code is implemented based on the metrics described in [NeuroMorphoVis](https://github.com/BlueBrain/NeuroMorphoVis/wiki/Analysis). 

* The values of the colormaps used to generate the color-coded projections are obtained from the [matplotlib](https://matplotlib.org/) library.

* The H5 morphology samples are available with permissions from the [Blue Brain Project, Ecole Polytechnique Federale de Lausanne (EPFL)](https://www.epfl.ch/research/domains/bluebrain/).

* The SWC morphology samples of neurons and astrocytes are available from NeuroMorpho.Org. [NeuroMorpho.Org](neuromorpho.org) is a centrally curated inventory or repository of digitally reconstructed neurons associated with peer-reviewed publications.

* The SWC morphology samples of brain vasculature are available from the [Brain Vasculature (BraVa) database](http://cng.gmu.edu/brava). 
The Brain Vasculature (BraVa) database contains digital reconstructions of the human brain arterial 
arborizations from 61 healthy adult subjects along with extracted morphological measurements.
The arterial arborizations include the six major trees stemming from the circle of Willis, 
namely: the left and right Anterior Cerebral Arteries (ACAs), Middle Cerebral Arteries (MCAs), 
and Posterior Cerebral Arteries (PCAs). 
Citation: [Susan N. Wright, Peter Kochunov, Fernando Mut Maurizio Bergamino, Kerry M. Brown, John C. Mazziotta, 
Arthur W. Toga, Juan R. Cebral, Giorgio A. Ascoli. 
Digital reconstruction and morphometric analysis of human brain arterial vasculature from magnetic 
resonance angiography. NeuroImage, 82, 170-181, (2013)](http://dx.doi.org/10.1016/j.neuroimage.2013.05.089). 

* The VMV vascular morphologies are available with permissions from [Pablo Blinder, Department of Neurobiology, Faculty of Life Sciences at Tel Aviv University](https://english.tau.ac.il/profile/pb). 

* The H5 vascular morphologies are available with permissions from the 
[Blue Brain Project](https://bluebrain.epfl.ch/page-52063.html), 
[Ecole Polytechnique Federale de Lausanne (EPFL)](https://www.epfl.ch/). The original dataset 
courtesy of [Bruno Weber](https://www.neuroscience.uzh.ch/en/about/people/steering/Weber.html), 
University of Zurich, Switzerland.

* Cortical meshes (in H5 format) are courtesy of the [MICrONS Consortium](https://www.microns-explorer.org/) including 
[Seung Lab](http://seunglab.org/), 
[Brain-map.org](https://brain-map.org) - [Allen Institute for Brain Science](https://alleninstitute.org/), 
[Tolias Lab](https://toliaslab.org/) and [IARPA Microns](https://www.iarpa.gov/index.php/research-programs/microns). 
    * The cortical mm^3 datasets are available from the following pulication:
        * MICrONs Consortium et al. Functional connectomics spanning multiple areas of mouse visual cortex. bioRxiv 2021.07.28.454025; doi: https://doi.org/10.1101/2021.07.28.454025
    * Layer 2/3 datasets are available from the following publications:
        * Dorkenwald, S., Turner, N.L., Macrina, T., Lee, K., Lu, R., Wu, J., Bodor, A.L., Bleckert, A.A., Brittain, D., Kemnitz, N., et al. (2019). Binary and analog variation of synapses between cortical pyramidal neurons. bioRxiv 2019.12.29.890319; doi: https://doi.org/10.1101/2019.12.29.890319
        * Schneider-Mizell, C. Bodor, A.L., Collman, F.  Brittain,D. Bleckert, AA, Dorkenwald, S., Turner N.L. Macrina, T.  Lee, K. Lu, R.  Wu, J. et al. (2020)  Chandelier cell anatomy and function suggest a variably distributed but common signal. bioRxiv 2020.03.31.018952v1; doi: https://doi.org/10.1101/2020.03.31.018952



Full attributions and acknowledgements are available in the [ACKNOWLEDGEMENTS file](https://github.com/BlueBrain/Ultraliser/blob/master/ACKNOWLEDGEMENTS.md).

# Contact

For more information on _Ultraliser_, comments or suggestions, please contact:

__Marwan Abdellah__  
Scientific Visualiation Engineer  
Blue Brain Project  
[marwan.abdellah@epfl.ch](marwan.abdellah@epfl.ch) 
 
__Felix Schürmann__  
Co-director of the Blue Brain Project    
[felix.schuermann@epfl.ch](samuel.lapere@epfl.ch) 

Should you have any questions concerning press enquiriries, please contact:

__Kate Mullins__  
Communications  
Blue Brain Project  
[kate.mullins@epfl.ch](kate.mullins@epfl.ch)

# References

1. YU, Zeyun, HOLST, Michael J., CHENG, Yuhui, et al. Feature-preserving adaptive mesh generation for molecular shape modeling and simulation. Journal of Molecular Graphics and Modelling, 2008, vol. 26, no 8, p. 1370-1380.

2. ATTENE, Marco. A lightweight approach to repairing digitized polygon meshes. The visual computer, 2010, vol. 26, no 11, p. 1393-1406.

3. ABDELLAH, Marwan, HERNANDO, Juan, ANTILLE, Nicolas, et al. Reconstruction and visualization of large-scale volumetric models of neocortical circuits for physically-plausible in silico optical studies. BMC bioinformatics, 2017, vol. 18, no 10, p. 39-50.

4. LORENSEN, William E. et CLINE, Harvey E. Marching cubes: A high resolution 3D surface construction algorithm. ACM siggraph computer graphics, 1987, vol. 21, no 4, p. 163-169.

5. NIELSON, Gregory M. Dual marching cubes. In : IEEE visualization 2004. IEEE, 2004. p. 489-496.

6. ABDELLAH, Marwan, HERNANDO, Juan, EILEMANN, Stefan, et al. NeuroMorphoVis: a collaborative framework for analysis and visualization of neuronal morphology skeletons reconstructed from microscopy stacks. Bioinformatics, 2018, vol. 34, no 13, p. i574-i582.

7. ABDELLAH, Marwan, GUERRERO, Nadir Román, LAPERE, Samuel, et al. Interactive visualization and analysis of morphological skeletons of brain vasculature networks with VessMorphoVis. Bioinformatics, 2020, vol. 36, no Supplement_1, p. i534-i541.

8. ZISIS, Eleftherios, KELLER, Daniel, KANARI, Lida, et al. Digital reconstruction of the neuro-glia-vascular architecture. Cerebral Cortex, 2021, vol. 31, no 12, p. 5686-5703.

9. ABDELLAH, Marwan, FONI, Alessandro, ZISIS, Eleftherios, et al. Metaball skinning of synthetic astroglial morphologies into realistic mesh models for in silico simulations and visual analytics. Bioinformatics, 2021, vol. 37, no Supplement_1, p. i426-i433.

10. ABDELLAH, Marwan. In silico brain imaging: physically-plausible methods for visualizing neocortical microcircuitry. EPFL, 2017.

--- 

<p align="center">
        <img src="docs/images/logos/epfl-logo.png" width=200>
</p>
 
<p align="center">Copyright (c) 2022 Blue Brain Project/EPFL</p>
