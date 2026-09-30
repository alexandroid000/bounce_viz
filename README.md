## Bounce Visualization
A collection of python utilities for visibility-based decomposition of 
polygons, and strategy generation for mobile "bouncing" robots.

## Installation and Setup

1. Clone this repository with submodules:

   ```bash
   git clone --recurse-submodules https://github.com/alexandroid000/bounce_viz.git
   ```

2. Make sure you have Python

3. Install python library dependencies

   ```bash
   pip3 install -r requirements.txt
   ```
    requirements.txt is a file in the root directory of this project. You can
    also use this list to install manually or convert to the python package manager
    of your choice.

## Quick Start

   ```bash
   cd $bounce_viz
   cd test
   ./run_sim.py
   ```

## Generate Documentations

   ```bash
   cd $bounce_viz
   cd docs
   make apidoc
   # to view the documentation
   open ./build/html/bounce_viz.html
   ```
