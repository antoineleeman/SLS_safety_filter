# Predictive safety filter using system level synthesis
This repository contains the MATLAB code that accompanies the research paper:
> Leeman AP., Köhler J., Bennani S., Zeilinger MN. “Predictive safety filter using system level synthesis.” Accepted L4DC2023 (2023)

![Project Image](fig1.pdf)

The paper is freely available on [arXiv](https://arxiv.org/abs/2212.02111).

The code in this repository demonstrates the implementation of the algorithms and models discussed in the paper. To run this code, you will need MATLAB, the Multi-Parametric Toolbox 3 (MPT3) and MOSEK.

## Prerequisites
- MATLAB (tested with version R2020b)
- Multi-Parametric Toolbox 3 (MPT3)
- MOSEK

## Installation
1. Download and install MATLAB from the [official website](https://www.mathworks.com/products/matlab.html).

2. Install Mosek by following the instructions from the [official Mosek documentation](https://docs.mosek.com/9.2/install/installation.html). 

3. Install the Multi-Parametric Toolbox 3 (MPT3) by following the instructions from the [official MPT3 documentation](https://www.mpt3.org/). In summary, you will need to:


    a. Download the latest release of MPT3 from the MPT3 GitHub repository.

    b. Extract the downloaded archive and add the extracted folder to your MATLAB path by running the following command in the MATLAB Command Window:
    
        
        addpath(genpath('/path/to/mpt3/folder'));
        
    
    c. Install the required solvers by running the following command in the MATLAB Command Window:
    
        
        mpt_init;
        
    
4. Clone this repository or download the code as a ZIP archive and extract it to a folder of your choice.

5. Add the code folder to your MATLAB path by running the following command in the MATLAB Command Window:
    
        addpath('/path/to/your/code/folder');
    
## Usage

Run the main script (i.e., main.m) to execute the algorithms and models discussed in the paper. Refer to the comments and documentation within the code for further details on the implementation.

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite our paper:
  ```
  @article{leeman2022predictive,
  title={Predictive safety filter using system level synthesis},
  author={Leeman, Antoine P and K{\"o}hler, Johannes and Benanni, Samir and Zeilinger, Melanie N},
  journal={arXiv preprint arXiv:2212.02111},
  year={2022}
  }
  ```
  

## Support and Contact

For any questions or issues related to this code, please contact the author:

- Antoine Leeman: aleeman(at)ethz(dot)ch

We appreciate any feedback, bug reports, or suggestions for improvements.
