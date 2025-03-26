# Monte Carlo Simulation for Fokker-Planck Equation with Metal.jl

This repository implements a simple Monte Carlo simulation framework for solving the Fokker-Planck equation using Julia and Metal.jl for GPU acceleration. The project includes tools for simulating particle trajectories, computing exit statistics, and visualizing results.

## Features

- **GPU-accelerated Monte Carlo simulations** using Metal.jl. If you happen to have a compatible GPU, you can leverage its power to speed up simulations.
- **GPU friendly random number generation** for Metal.jl. The project includes utilities for generating random numbers on the GPU.
- **Visualization tools** for generating plots and animations of simulation results.

## Example

Below is an example GIF demonstrating the simulation results:

![Simulation Results](square_boundary.png)

## Repository Structure

- `[Julia][Metal]Monte-Carlo-KFP.ipynb`: Main Jupyter notebook demonstrating the simulation workflow, including setup, execution, and visualization.
- `MCKernel.jl`: Contains GPU kernel definitions for Monte Carlo simulations.
- `plot_helper.jl`: Helper functions for generating plots and animations.
- `RandomMetal.jl`: Utilities for random number generation and Metal.jl integration.

## Usage

To run the simulation, you need to have Julia installed on your system. You can download Julia from the [official website](https://julialang.org/downloads/). To run the Jupyter notebook, you also need to have Jupyter and IJulia installed.

After installing the necessary dependencies, you can open the Jupyter notebook and follow the instructions to set up the simulation and visualize the results.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
