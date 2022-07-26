# Introduction

_Ultraliser_ is an __unconditionally robust__ and __high performance__ framework dedicated primarily to **_in silico_ neuroscience** research. _Ultraliser_ is capable of generating high fidelity and multiscale 3D neuroscientific models - such as: nuclei, mitochondria, endoplasmic reticula, neurons, astrocytes, pericytes, neuronal branches with dendritic spines, minicolumns with thousands of neurons and large networks of cerebral vasculature - with realistic geometries. 

_Ultraliser_ implements an effective voxelization-based remeshing engine that can rasterize non-watertight surface meshes - in the form of triangular soups - into high resolution volumes, with which we can reconstruct topologically accurate, adaptively optimized and watertight surface manifolds. 

In addition to their importance for accurate quantitative analysis, resulting models are primarily intended to automate the process of conducting supercomputer simulations of neuroscience experiments; complementing _in vivo_ and _in vitro_ techniques. 

Watertight triangular meshes are used for (i) performing 3D particle simulations, (ii) mesh-based skeletonization, in which accurate morphologies of cellular structures are obtained for performing 1D compartmental simulations and (iii) tetrahedralization, in which we can generate tetrahedral volume meshes for 3D reaction-diffusion simulations.
Annotated volumetric tissue models are also used in _in silico_ imaging studies, where we can simulate optical imaging experiments with brightfield or fluorescence microscopy<sup>10</sup>. 
A high-level overview of _Ultraliser's_ workflow is graphically illustrated in __Figure 1__.

<figure>
<p align="center">
  <img width="1000" src="images/workflows/workflow.jpeg">
</p>
<figcaption>
<p align="justify">
<strong>Figure 1</strong> A high-level overview of Ultraliser's workflow showing its subsequent stages (surface and solid voxelization, mesh reconstruction and optimization) including the different data formats that can be processed within its workflow (masks, gray-scale volumes, polygonal surface meshes, morphological skeletons and tetrahedral volumetric meshes).
Ultraliser creates high-fidelity, adaptively optimized and watertight surface manifolds and large-scale annotated volumes. The surface meshes are used directly for molecular simulations, and can be further processed to create volumetric (tetrahedral or hexahedral) meshes for performing reaction-diffusion simulation. They can be also used by mesh-based skeletonization applications to skeletonize accurate morphological counterparts. The annotated volumes are used in <em>in silico</em> optical imaging experiments<sup>10</sup>.
</p>
</figcaption>
</figure>

## Features 

+ Reconstruction of high fidelity, optimized<sup>1</sup> and two-manifold watertight<sup>2</sup> triangular mesh models from non-watertight inputs represented by polygonal soups. 
+ Surface mesh smoothing and optimization using [Laplacian operators](https://en.wikipedia.org/wiki/Laplacian_smoothing) and feature-preserving adaptive mesh optimization<sup>1</sup>.
+ Reconstruction of large-scale volumetric models<sup>3</sup> from non-watertight input meshes using high performance [surface and solid voxelization](Mesh-Voxelization). 
+ Reconstruction of optimized and smooth surface meshes from input volumes using parallel implementations of the standard marching cubes<sup>4</sup> algorithm and the advanced dual marching cubes<sup>5</sup> algorithm. 
+ Reconstruction of optimized and smooth surface meshes from input binary masks of segmened data. 
+ Reconstruction of geometrically realistic watertight mesh models of spiny neurons from corresponding morphological skeletons<sup>6</sup> .
+ Reconstruction of geometrically realistic watertight mesh models of complete astroglial cells<sup>8</sup>  (with endfeet) from input morphological skeletons and endfeet surface patches<sup>9</sup> .        
+ Reconstruction of high-fidelity, optimized and multi-partitioned vascular meshes from fragmented and large-scale vascular network graphs<sup>7</sup>. 
+ Morphology, mesh and volume quantiative and qualitative analysis. 
+ Generation of color-coded multi-axis projections of spatial data (morphologies, meshes and volumes) for visual analytics. 


## Installation

### Software Dependencies 

* [OpenMP](https://en.wikipedia.org/wiki/OpenMP), a multi-threading library for parallel processing on multi-core CPUs. 
* [libTIFF](http://www.libtiff.org/), which gives support for the Tag Image File Format (TIFF), a widely used format for storing image data.
* [libhdf5](https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html), or the Hierarchical Data Format 5 (HDF5) library for stroring data.
* [Eigen3](https://eigen.tuxfamily.org/), a template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
* [BZip2](https://www.sourceware.org/bzip2/), a high quality data compressor.  
* [ZLIB](https://docs.python.org/3/library/zlib.html), for data compression. 
* [FMT](https://github.com/fmtlib/fmt), a formatting library providing a fast and safe alternative to C stdio and C++ iostreams.
* [GLM](https://github.com/g-truc/glm), a header only C++ mathematics library for graphics software based on the OpenGL Shading Language (GLSL) specifications.

### Supported Operating Systems 

Ultraliser has been tested on Unix-based operating systems including:

* Ubuntu 16.04, Ubuntu 18.04, Ubuntu 20.04
* RedHat 
* macOS 10.12 Sierra, 10.13 High Sierra,  10.14 Mojave, 10.15 Catalina


# Known Bugs or Feature Requests

Please refer to the [github issue tracker](https://github.com/BlueBrain/Ultraliser/issues?utf8=%E2%9C%93&q=) for fixed and open bugs. User can also report any bugs and request new features needed for their research. We are happy to provide direct [support](#contact) . 


# License 
_Ultraliser_ is available to download and use under the GNU General Public License, version 3 ([GPL](https://www.gnu.org/licenses/gpl.html), or “free software”). 
The code is open sourced with approval from the open sourcing committee and principal coordinators of the Blue Brain Project in March 2021. 

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
_Ultraliser_ is developed by the Visualization team at the [Blue Brain Project](https://bluebrain.epfl.ch/page-52063.html), [Ecole Polytechnique Federale de Lausanne (EPFL)](https://www.epfl.ch/). Financial support was provided by competitive research funding from [King Abdullah University of Science and Technology (KAUST)](https://www.kaust.edu.sa/en).

# Attributions

* The volume reconstruction code is an extension to the work of [Marwan Abdellah's](http://marwan-abdellah.com/) [PhD (In silico Brain Imaging: Physically-plausible Methods for Visualizing Neocortical Microcircuitry)](https://infoscience.epfl.ch/record/232444?ln=en). 

* The mesh optimization code in _Ultraliser_ is based on the routiens provided by the [GAMer (Geometry-preserving Adaptive MeshER) library](http://fetk.org/codes/gamer/). GAMer is developed by _Z. Yu, M. Holst, Y. Cheng, and J.A. McCammon_, and can be redistributed and modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.

* The watertighness verification code in _Ultraliser_ is based on an extended version of the [MeshFix library](https://github.com/MarcoAttene/MeshFix-V2.1). MeshFix is developed by _Marco Attene_, Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche, Sezione di Genova, IMATI-GE / CNR, and can be redistributed and modified under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or any later version.    

* The mesh analysis code is implemented based on the metrics described in the [The Verdict Geometric Quality Library](https://coreform.com/papers/verdict_quality_library.pdf). 

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

# Navigation


# Public

<p align="center">
        <img src="docs/images/logos/epfl-logo.jpg" width=200>
</p>
 
 
Copyright © 2022 Blue Brain Project/EPFL
