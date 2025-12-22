# Leapfrog Integration Simulation (Fortran)

This project implements a **3D particle simulation** using the **Leapfrog integration method**.  
The code is written in **Fortran 90** and simulates the motion of particles under Newtonian gravitational forces.

---

# Project Structure
- `ex1.f90`: main program (leapfrog integrator)
- `geometry.f90`: module defining `vector3d` and `point3d` types, and charging mathematical operations between points and vectors.  
- `particle.f90`: module defining the `particle3d` type, specified by a mass and the position and velocity components.

- `input.dat`: input file with initial conditions and parameters.
- `output.dat`: output file with the particles positions at every time step.
- `Makefile`: compile the project, generating an output file.

- `figures.py`: *Python* script that generates an animation for the particles from the `output.dat` file.
- `animation.gif`: animation generated from the simulation output
- `requirements.txt`: *Python* enviroment, with all the packages


---

# Compilation

Make sure you have a Fortran compiler installed (e.g., `gfortran`).

To compile the program:

```bash
make
```

To clean up the generated object (.o) and module (.mod) files, run:
```bash
make clean
```


The `input.dat` file must contain:

| Variable | Description |
| :--- | :--- |
| `dt` | Integration time step |
| `dt_out` | Time interval between data outputs |
| `t_end` |  Total simulation time |
| `n` | Number of particles |
| `m x y z vx vy vz` | For each particle: mass, initial position (x,y,z) and velocity (vx, vy, vz) |


The output file contains the simulation time followed by the 3D positions (x, y, z) of all particles at that time.

You can visualize the simulation using the Python script provided in Graphics/figures.py.
```bash
python Graphics/figures.py
```

This script reads output.dat and generates trajectory plots or animations (animation.gif) inside the Graphics/ folder.

To reproduce the same environment, use the requirements.txt file:
```bash
conda create -n neutron python=3.13.9
conda activate neutron
pip install -r Graphics/requirements.txt
```

