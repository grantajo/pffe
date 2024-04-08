#define DIM 2

#include <iostream>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "Set/Set.H"
#include "Element/CST.H"
#include "Element/LST.H"
#include "Element/Q4.H"
#include "Element/Q9.H"
#include "Element/Test.H"
#include "Model/Isotropic.H"
#include "Model/Test.H"
#include "Mesh/Unstructured.H"
#include "Mesh/Test.H"


int main(int argc, char **argv)
{
    //
    // Create an unstructured mesh
    //
    Mesh::Unstructured<Model::Isotropic> mymesh("platehole_lst.vtk");

    //
    // Create a displacement vector to store nodal displacements
    //
    Eigen::VectorXd  disp = Eigen::VectorXd::Zero(mymesh.size());

    //
    // Calculate the first and second energy derivatives
    //
    Eigen::VectorXd DW = mymesh.DW(disp);
    Eigen::SparseMatrix<double> DDW = mymesh.DDW(disp);

    //
    // Modify the stiffness matrix to apply essential (displacement)
    // boundary conditions at the leftmost and rightmost nodes.
    //
    Set::Scalar du = 0.01; // the prescribed displacement
    for (int i = 0 ; i < mymesh.Points.size(); i++)
    {
        // Do this for all points on the left side (x=0)
        if (fabs(mymesh.Points[i](0) - 0.0) < 1E-8)
        {
            // x component of displacement
            DDW.row(2*i) *= 0;                      // Zero out this row of the stiffness matrix
            DDW.coeffRef(2*i,2*i) = 1.0;            // and then set the diagonal to one.
            DW(2*i) = -du;                          // Set the right-hand side to "-du"

            // y component of displacement ... same thing as z but with disp = 0
            DDW.row(2*i + 1) *= 0;
            DDW.coeffRef(2*i + 1,2*i + 1) = 1.0;
            DW(2*i + 1) = 0.0;
        }


        // Do this for all points on the right side (x1)
        // Same thing as for the right side but with positive displacement rather than negative.
        if (fabs(mymesh.Points[i](0) - 1.0) < 1E-8)
        {
            DDW.row(2*i) *= 0;
            DDW.coeffRef(2*i,2*i) = 1.0;
            DW(2*i) = du;

            DDW.row(DIM*i + 1) *= 0;
            DDW.coeffRef(DIM*i + 1,DIM*i + 1) = 1.0;
            DW(DIM*i + 1) = 0.0;
        }
    }


    // Use the Eigen library to solve DDW u = DW
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(DDW);
    disp = solver.solve(DW);

    // Use the print function to output displacements.
    mymesh.Print("outputfile.vtk", disp);


    // TODO : Run this code as-is for the following mesh files
    //        - platehole_cst.vtk
    //        - platehole_q4.vtk
    //        - platehole_lst.vtk
    //        - platehole_q9.vtk
    //        and report the results.

    // TODO : Modify the boundary conditions to apply a cantilever displacement
    //        by setting the displacements to zero along the left edge, and
    //        the y-displacement of the lower right-hand corner to -du.
    //        Plot and report the results.

}
