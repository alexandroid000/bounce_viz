## Bounce Visualization
A collection of python utilities for visibility-based decomposition of 
polygons, and strategy generation for mobile "bouncing" robots.

## Installation and Setup

1. Clone this repository with submodules:

   ```bash
   git clone --recurse-submodules https://github.com/alexandroid000/bounce_viz.git
   ```

2. Install python library dependencies

   ```bash
   uv sync
   ```
   or use `requirements.txt` for legacy python environment managers.

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
