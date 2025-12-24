Angeles Moreno Guedes -- ma/course_exercise_2


This is a Barnes-Hut algorithm to do a N body simulation. There are two folders: one it is parallelise with openmp and the other with MPI. It can be run in serial executing the code of the openmp folder or in mpi writing only one processor. <br>

In each folder we have 7 files: the ex2.f90 is the main program where we use all the modules and we read the input file. The geometry module that includes several functions to compute vector and points operations. The particles module, that defines a derived type to describe the principal characteristics of a particle to compute its movement. The tree module that contains the Barnes-Hut algorithm. The N_body_generator program, is a generator of the random position of the particles (the code is the same of the notes). Finally, the folders has an out.py that you can run to make an animation from the data given by output.dat. <br>


In order to use the code, is important to first compile it with the makefile, using make in the directory where the code is: <br>

`make` <br>

However, if we use make like this in the openmp folder we are going to run the code in serial. To use openmp we have to write: <br>

`make USE_OMP=1`<br>

Then, to compile the generator of the input.dat file:<br>

`make generator` <br>

And to run it there are two options: <br>

`make gen` <br>
`./generator`<br>

Finally, to run the code there are two alternatives (in the openmp folder): <br>

`make run`<br>
`./ex2 <inputfile name>` <br>

And in the mpi folder you can select the number of processors in the make file or if you run it like in the following line instead of `make run`: <br>

`mpirun -np 8 ./ex2_mpi input.dat` <br>

It is easy to run in serial by choosing `-np 1` <br>

If you want to erased the input.dat, output.dat and the executables with (recommended): <br>

`make clean` <br>


The input.dat with the initial conditions has the following shape: <br>

dt    (timestep) <br> 
dt_out     (timestep to print the result)<br> 
t     (simulation time)<br> 
n      (number of particles)<br> 
m1 x1 y1 z1 vx1 vy1 vz1 (mass, initial position, initial velocity of particle1)<br> 
.<br>
.<br>
mn xn yn zn vxn vyn vzn        (mass, initial position, initial velocity of particlen)<br> 


Finally, the output will be an output.dat with the following shape: <br>

time    p1x     p1y     p1z     ....    pnx     pny     pnz<br>

It will appear in the same folder when the program ends. Besides, when the simulation ends in the terminal it is printed the total time it last <br>
In my computer, with 300 particles and t_end = 10s (it depends of the specifications of each computer and the set of particles) it last: <br>

Serial :  7.15 s <br>
OPMP : 3.55 s   <br>
MPI : 8.03 s  (with 4 procesors) <br>

For 1000 particles and t_end = 10s <br>

Serial : 40s   <br>
OPMP :  19s  <br>
MPI :   41.4 s  (with 4 procesors) <br>

For 3000 particles and t_end = 10s <br>

Serial : 181.30 s  <br>
OPMP :  45.04 <br>
MPI :  201 s approx 3.35 min (with 4 procesors) <br>


For 10000 particles and t_end = 10s <br>

Serial : 877.3 s aprox 15 min  <br>
OPMP :  185.42 s aprox 3 min <br>
MPI :  926.81 s aprox 15.44 min (with 4 procesors) <br>


### Performance Analysis: Why is MPI slower?

As observed in the results, **OpenMP** offers the best performance, scaling well with the number of particles. However, the **MPI** version is consistently slower than the OpenMP version and slower than the serial version. This could happen by these reasons:

1.  **The MPI_ALLREDUCE of the notes parallelize evrything not only the calculation of the forces:**
    The time spent packaging data, sending it, and synchronizing the processes (communication latency) exceeds the time saved by dividing the calculation work.

2.  **Shared Memory vs. Distributed Memory:**
    These tests were run on a single computer.
    * **OpenMP** uses *Shared Memory*: Threads access the data (positions/masses) directly without copying it. This is very efficient on a single computer.
    * **MPI** uses *Distributed Memory*: Even on a single computer, MPI simulates a network. It has to copy data from one process's memory space to another.


**Conclusion:**
OpenMP is the optimal choice for this simulation on a single workstation. MPI would only start to show advantages if the number of particles were massive or if the simulation were distributed across a cluster.
