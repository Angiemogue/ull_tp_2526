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

Serial :  7.36 s <br>
OPMP : 9.34 s   <br>
MPI : 3.37 s  (with 4 procesors) <br>

For 1000 particles and t_end = 10s <br>

Serial : 44.03s   <br>
OPMP :  22.94 s  <br>
MPI :   13.90 s  (with 4 procesors) <br>

For 3000 particles and t_end = 10s <br>

Serial : 196.58 s  <br>
OPMP :  50.72 s <br>
MPI :  21.08 s (with 4 procesors) <br>


For 10000 particles and t_end = 10s <br>

Serial : 872.08 s aprox 15 min  <br>
OPMP :  170.36 s aprox 3 min <br>
MPI :  276.71 s aprox 5 min (with 4 procesors) <br>


**Conclusion:**
As observed in the results, **MPI** (Distributed Memory) offers the best performance (in general), followed by **OpenMP** (Shared Memory) that scales well with the number of particles. In the small dataset (N=300) **OpenMP** is the slower one due to the small number of particles. However, in the largest dataset (N=10000) is the becomes the fastest because it benefits of the shared memory against the distribuited one with **MPI**. The **MPI** is harder to program but offers more options to parallelize and performance. And it could offer even better performance in a cluster of computers. However, OpenMP can be a great option in a single computer since its performance is not so bad in comparation with the serial running and in this case benefits with the single workstation comparation
