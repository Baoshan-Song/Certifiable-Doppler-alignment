# Certifiable Doppler Alignment

This repository provides the official MATLAB implementation for the paper:  
**"Certifiable Alignment of GNSS and Local Frames via Lagrangian Duality"**

The project introduces a certifiably globally optimal framework to align Global Navigation Satellite System (GNSS) frames with local robot frames (such as VIO or LIO) by leveraging Doppler measurements and Lagrangian duality.

## Prerequisites

### 1. External Dependencies (Required)
This repository relies on **MatRTKLIB** for GNSS data processing and coordinate transformations. You **must** download it and add it to your MATLAB path before running the code:

* **MatRTKLIB**: [https://github.com/taroz/MatRTKLIB](https://github.com/taroz/MatRTKLIB)

### 2. Software & Solvers
* **MATLAB**: Tested on R2022a or later.
* **CVX**: Required for solving the Semidefinite Programming (SDP) relaxation. 
* **Solvers**: Professional solvers like **MOSEK** or **SeDuMi** are highly recommended for better numerical stability and performance.

## Getting Started

### Installation
1.  Clone this repository:
    ```bash
    git clone [https://github.com/Baoshan-Song/Certifiable-Doppler-alignment.git](https://github.com/Baoshan-Song/Certifiable-Doppler-alignment.git)
    ```
2.  Download/Clone [MatRTKLIB](https://github.com/taroz/MatRTKLIB) and place it in your workspace.
3.  Open MATLAB and add both the `Certifiable-Doppler-alignment` and `MatRTKLIB` folders (including all subfolders) to your **MATLAB Path**.

### Running the Example
To verify the algorithm and observe the certifiable global optimality, we provide a comprehensive simulation script:

1.  In the MATLAB Command Window, run:
    ```matlab
    example_simulation
    ```

The script `example_simulation.m` generates a synthetic trajectory, simulates noisy Doppler measurements, and performs the certifiable alignment. It will output the rotation and translation errors alongside the **Duality Gap**, which serves as the "certificate" for global optimality.

## Methodology
The proposed method formulates the frame alignment problem as a Quadratically Constrained Quadratic Program (QCQP). By utilizing Lagrangian duality and SDP relaxation, the algorithm ensures that the obtained solution is the global optimum rather than a local minimum.

## Citation

If you find this code or the paper useful for your research, please cite our work:

```bibtex
@article{song2025certifiable,
  title={Certifiable Alignment of GNSS and Local Frames via Lagrangian Duality},
  author={Song, Baoshan and Giamou, Matthew and Yan, Penggao and Xia, Chunxi and Hsu, Li-Ta},
  journal={arXiv preprint arXiv:2512.20931},
  year={2025}
}
