# Dynamic Simulation Project

This repository collects a set of experiments for simulating vibrating mechanical systems in classical mechanics.  The focus is on modelling double and triple pendulum systems, integrating their equations of motion numerically, and visualising the resulting motion.

## Project structure

```
├── Animations/
│   ├── Autocontrolpendulul.py      # Triple-pendulum animation with damping
│   ├── Simulation-1.ipynb          # Notebook exploring the double pendulum
│   ├── solvingtripendululm.py      # SymPy derivation of triple-pendulum dynamics
│   ├── test_script.py              # Quick dependency check for SymPy
│   ├── pendulum_animation.mp4      # Sample double-pendulum animation export
│   └── tripendulum_animation.mp4   # Sample triple-pendulum animation export
├── simulation.py/
│   ├── Animations/Simulation-1.py  # Alternate double-pendulum animation script
│   ├── Autocontrolpendulul.py      # Copy of the triple-pendulum animation
│   ├── solvingtripendululm.py      # Copy of the symbolic derivation script
│   └── ...                         # Supporting assets exported from notebooks
├── Untitled.ipynb                  # Scratch notebook used during development
└── README.md
```

> The `simulation.py/` directory is an exported Jupyter notebook bundle.  The top-level `Animations/` folder contains the actively edited source files.

## Features

* **Double pendulum simulation** – integrates the coupled equations of motion using a simple Euler stepper and renders the motion with Matplotlib (`simulation.py/Animations/Simulation-1.py`).
* **Triple pendulum with damping** – extends the dynamics to three links and introduces configurable viscous damping for each joint (`Animations/Autocontrolpendulul.py`).
* **Symbolic derivation** – derives the equations of motion for a cart-triple pendulum using SymPy’s mechanics toolkit, producing the mass matrix and forcing vector (`Animations/solvingtripendululm.py`).
* **Animation exports** – example `.mp4` files demonstrate the expected results when the simulations are run locally with FFmpeg installed.

## Getting started

1. **Create and activate a virtual environment (recommended)**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows use: .venv\\Scripts\\activate
   ```
2. **Install dependencies**
   ```bash
   pip install numpy matplotlib sympy
   ```

## Running the simulations

### Double pendulum animation

```
python simulation.py/Animations/Simulation-1.py
```

The script integrates the double pendulum for 10 seconds with a fixed time step and opens an animated Matplotlib window.  If FFmpeg is available, an `.mp4` animation is also written to the path configured near the end of the script.

### Triple pendulum with damping

```
python Animations/Autocontrolpendulul.py
```

This script simulates a three-link pendulum with independent damping coefficients.  Adjust the constants at the top of the file to explore different mass distributions, link lengths, or damping ratios.

### Symbolic dynamics derivation

```
python Animations/solvingtripendululm.py
```

The script prints the mass matrix and forcing vector produced by Kane’s Method.  These expressions can be copied into other analysis pipelines or used to verify numerical implementations.

## Extending the project

* Experiment with alternative integrators (e.g. Runge–Kutta methods) to improve accuracy over Euler integration.
* Add energy diagnostics to track kinetic and potential energy exchange during the motion.
* Couple the derived equations with control inputs to investigate stabilisation of the cart–pendulum system.

## Requirements and notes

* Python 3.9+ is recommended.
* FFmpeg is required if you want to export the Matplotlib animations to video files.
* Jupyter notebooks in the repository contain intermediate derivations and visualisations.  Launch them with `jupyter notebook` if you prefer an interactive workflow.

Happy simulating!

