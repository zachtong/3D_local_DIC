# 3D-Stereo-ALDIC

## Overview
3D-Local-DIC is a MATLAB-based implementation for full-field 3D deformation measurement using stereo vision and a local, subset-based Digital Image Correlation (DIC) approach. This tool provides robust and straightforward 3D displacement and strain measurements for experimental mechanics applications.

## Features
- User-friendly tool for 3D deformation measurement.
- Fast and robust algorithm using distributed parallel computing
- Adaptive quadtree mesh refinement for optimal spatial resolution
- Support for both accumulative mode nd incremental DIC mode (for large deformation analysis)
- Multiple camera calibration format support (MATLAB, MatchID, MCC, DICe, OpenCorr)
- Comprehensive strain field computation
- Advanced visualization tools for results analysis

## Requirements
- Image Processing Toolbox
- Computer Vision Toolbox
- C++ compiler for MEX file compilation

## Code manual 
[Under development]

## Code demo videos
[Under development]


## Citation
* [1] Tong, Z, et. al. 3D Stereo Adaptive Mesh Augmented Lagrangian Digital Image Correlation. https://doi.org/10.21203/rs.3.rs-5507109/v1
* [2] Tong, Z, et, al. Machine Learning-Aided Spatial Adaptation for Improved Digital Image Correlation Analysis of Complex Geometries. https://doi.org/10.21203/rs.3.rs-5566473/v1

Previous works:
* [3] Yang, J. and Bhattacharya, K. Augmented Lagrangian Digital Image Correlation. Exp.Mech. 59: 187, 2018. https://doi.org/10.1007/s11340-018-00457-0.
* [4] Yang, J. 2D_ALDIC. https://github.com/jyang526843/2D_ALDIC
* [5] Yang, J. FEM-based Global DIC: https://github.com/jyang526843/2D_FE_Global_DIC
* [6] Yang, J. ALDVC code: https://github.com/FranckLab/ALDVC


## Usage
1. Clone this repository
2. Prepare stereo image pairs of your specimen
3. Run the main script:
```matlab
main_LocalDIC_3D
```

## Input Data Format
- Image pairs from left and right cameras
- Camera calibration data (supported formats: MATLAB, MatchID, MCC, DICe, OpenCorr)
- ROI masks for analysis region definition

## Output Results
- 3D displacement fields
- Strain fields
- Visualization plots and figures
- Results saved in MATLAB .mat format

## Contributing
We welcome contributions to improve 3D-Stereo-ALDIC. Please feel free to submit pull requests or report issues.

## Citation
If you use this code in your research, please cite:
[3D Stereo Adaptive Mesh Augmented Lagrangian Digital Image Correlation]
https://www.researchsquare.com/article/rs-5507109/v1

## License
MIT License.

## Authors
- Zixiang (Zach) Tong (zachtong@utexas.edu) @ UT-Austin
- Jin Yang (jin.yang@austin.utexas.edu) @ UT-Austin