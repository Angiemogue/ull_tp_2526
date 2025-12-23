Angeles Moreno Guedes -- ma/course_exercise_2


This is a Barnes-Hut algorithm to do a N body simulation. There are two folders: one it is parallelise with openmp and the other with MPI. It can be run in serial executing the code of the openmp folder. <br>

In each folder we have 7 files: the ex2.f90 is the main program where we use all the modules and we read the input file. The geometry module that includes several functions to compute vector and points operations. The particles module, that defines a derived type to describe the principal characteristics of a particle to compute its movement. The tree module that contains the Barnes-Hut algorithm. The N_body_generator program, is a generator of the random position of the particles. Finally, the folders has an out.py that you can run to make an animation from the data given by output.dat. <br>


In order to use the code, is important to first compile it with the makefile, using make in the directory where the code is: <br>

`make` <br>

However, if we use make like this in the openmp folder we are going to run the code in serial to use openmp we have to write: <br>

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

`mpirun -np 4 ./ex2_mpi input.dat` <br>

If you want to erased the input.dat, output.dat and the executables with: <br>

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
In my computer, with 300 particles and t_end = 100s (it depends of the specifications of each computer) it last: <br>

Serial :  257.1508 s approx 4.28 min <br>
OPMP :  253.9883 s  approx 4.20 min <br>
MPI : 105.26721 s approx 2 min (with 4 proceses) <br>



