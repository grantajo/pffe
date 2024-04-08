# FEM Coding Assignment: Isoparametric Elements, Unstructured Meshes, FEM simulation


## Initial setup (same as before)

1. Check out the problem repository
    
         git clone https://github.com/solidsgroup/EM525-PS07.git
   
2. Change into the problem directory

         cd EM525-PS07
   
3. Use this command to install eigen (optional: you can skip this if you have eigen installed)

         make eigen

4. Now, compile the code

         make

   The code **will not compile** because there are missing files.
   These files must be copied as specified below before you can compile.


## Part 0: Copying previous implementation

Copy the following files directly from PS05. 

- `src/Element/Element.H`
- `src/Element/LST.H`
- `src/Element/CST.H`
- `src/Element/Q4.H`
- `src/Element/Q9.H`
- `src/Model/Isotropic.H`

If you do this correctly the code will compile.
However, it will **not** pass the tests.

## Part 1: Global stiffness matrix construction

### (1a) Update the Jacobian

Edi the following file

- `src/Element/Element.H`

In your implementation of the `Jacobian` function, change the node access from

    (*X0)[n]

to 

    (*X0)[id[n]]

The `id` array is a list of nodal indexes that each element has. 
For instance, a Q4 element might have 

    id = {0, 3, 29, 1}
	
That is, it is a Q4 element between the nodes indexed by 0, 3, 29, and 1.


### (1b) Implement the Mesh energy derivatives

Edit the following file

- `src/Mesh/Unstructured.H`

Update the `W`, `DW`, and `DDW` for Q4, LST, and Q9 elements, following the existing
implementation for CST elements.


### (1c) Run the tests

If you have implemented the above correctly, you can run the tests by

   ./bin/test

Make sure the tests pass before proceeding to the next part.

## Part 2: Finite element calculation

In this part we will run an actual FEM simulation.


### (2a) Update the Print function

Edit the following file

- `src/Mesh/Unstructured.H`

Update the `Print` function by copying in your previous implementation, and also adding the capability to output nodal displacements.

### (2b) Run the main program

Read through the code in the following file

- `src/main.cpp`

Run the finite element program by compiling and running

    ./bin/main

This will generate an output mesh called `outputfile.vtk` which you can view using Paraview or VisIt.
You should be able to visualize the displacement field to inspect the result.
If you cannot, you may need to update the Print function.

Run for all four mesh files, then report and contrast the results.


### (2c) Update the boundary conditions in the main program

Modify the boundary conditions in the `main.cpp` program to cantilever loading.
See the instructions at the bottom of `main.cpp`.
Re-run with the four meshes and report/contrast the results.


