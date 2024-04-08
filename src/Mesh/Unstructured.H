#include <fstream>
#include <cassert>
#include <filesystem>
#include "Element/Element.H"
#include "Element/CST.H"
#include "Element/LST.H"
#include "Element/Q4.H"
#include "Element/Q9.H"
#include "Mesh/Mesh.H"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/SparseCore"
#include "eigen3/Eigen/IterativeLinearSolvers"
#include "Model/Isotropic.H"
namespace Mesh
{
template<class MODEL>
class Unstructured : Mesh
{
public:
    Unstructured()
    {};

    Unstructured(std::string vtkfile)
    {
        //
        // TODO: Copy implementation from PS6
        //
    }

    void Print(std::string vtkfile,
               Eigen::VectorXd &disp)
    {
        //
        // TODO: Copy implementation from PS6, AND modify to print out
        //       tbe components of the displacement data.
        //
        // To do this, update the code so that the following block is
        // added to the output file:
        //
        //  POINT_DATA  [[ total number of vertices ]] 
        //  VECTORS disp float
        //  [x displacement of node 0] [y displacement of node 0] 0.0
        //  [x displacement of node 1] [y displacement of node 2] 0.0
        //   ...
        // 
                
    }
    double W(Eigen::VectorXd & u)
    {
        //
        // This function calculates the total strain energy of the entire
        // mesh, including all CST, Q4, LST, and Q9 elements.
        //
        // Input:  the input vector u is a *flattened* degree of freedom vector
        //             
        //           ┌       ┐
        //           │ u^0_1 │  <-- x displacement for node 0
        //           │ u^0_2 │  <-- y displacement for node 0
        //           │ u^1_1 │  <-- x displacement for node 1
        //           │ u^1_2 │         .. etc ..
        //    u =    │ u^1_1 │
        //           │ u^1_2 │
        //           │  ...  │
        //           │ u^n_1 │
        //           │ u^n_2 │
        //           └       ┘
        //
        //         To calculated energies using the individual elements.
        //         the array must be "unflattened".
        //
        // Output: the total energy of the entire mesh of elements
        //

        double ret = 0.0;

        
        // Add up the contributions from the CST nodes.
        // Loop over all CST elements:
        for (int e = 0; e < CSTs.size(); e++)
        {
            // Get the vector of IDs from CST element e
            auto & id = CSTs[e].getid();
            
            // Construct an array of vectors ("in") to feed into
            // the element's W function.
            std::array<Set::Vector2d, 3> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            
            // Add up the energy using the CST element's energy function.
            ret += CSTs[e].W(in);
        }
        
        //  
        // TODO: Following the above, implement for Q4, LST, and Q9 elements.
        //

        // ... your code goes here ...

        return ret;
    }

    Eigen::VectorXd DW(Eigen::VectorXd & u)
    {
        //
        // This function calculates the derivative of the total strain energy of 
        // the entire mesh, including all CST, Q4, LST, and Q9 elements.
        //
        // Input:  the input vector u is a flattened degree of freedom vector, identical
        //         to that for DW.
        //             
        // Output: The output vector is also a flattened vector, with components
        //         corresponding to nodal "forces"
        //
        //                ┌       ┐
        //                │ f^0_1 │  <-- x nodal "force" for node 0
        //                │ f^0_2 │  <-- y nodal "force" for node 0
        //                │ f^1_1 │  <-- x nodal "force" for node 1
        //                │ f^1_2 │         .. etc ..
        //    f = DW(u) = │ f^1_1 │
        //                │ f^1_2 │
        //                │  ...  │
        //                │ f^n_1 │
        //                │ f^n_2 │
        //                └       ┘
        //
        //         For this function, in addition to "unflatttening" the displacement
        //         vector, the force vector must be "re-flattened" to match the above
        //         format.
        //
        

        Eigen::VectorXd ret = Eigen::VectorXd::Zero(size());

        // Calculate contributions from all CST elements by
        // looping over the CSTs array.
        for (int e = 0; e < CSTs.size(); e++)
        {
            // This part is identical to that for W.
            auto & id = CSTs[e].getid();
            std::array<Set::Vector2d, 3> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);

            // Now we calculate the dw for each element, which
            // is an array of vectors.
            std::array<Set::Vector2d,3> dw = CSTs[e].DW(in);
            
            // Now, flatten the array by storing each component in
            // the corresponding node's slot.
            //   ** note ** it is important to use += here - a node's
            //              forces come from all elements that share
            //              that node!
            ret(2*id[0]  ) += dw[0](0);
            ret(2*id[0]+1) += dw[0](1);
            ret(2*id[1]  ) += dw[1](0);
            ret(2*id[1]+1) += dw[1](1);
            ret(2*id[2]  ) += dw[2](0);
            ret(2*id[2]+1) += dw[2](1);
        }
        
        //  
        // TODO: Following the above, implement for Q4, LST, and Q9 elements.
        //

        // ... your code goes here ...


        return ret;
    }

    Eigen::SparseMatrix<double> DDW(Eigen::VectorXd & u)
    {
        //
        // This function calculates the second derivative of the total strain energy of 
        // the entire mesh, including all CST, Q4, LST, and Q9 elements, which is
        // the elastic modulus tensor.
        //
        // Input:  the input vector u is a flattened degree of freedom vector, identical
        //         to that for DW.
        //             
        // Output: The output vector is a flattened, sparse matrix.
        //         The structure of the flattened matrix is the following:
        //
        //                ┌                                                                      ┐
        //                │ ddw^{00}_{00}   ddw^{00}_{01}    ddw^{01}_{00}    ddw^{01}_{00} .... │  
        //                │ ddw^{00}_{10}   ddw^{00}_{11}    ddw^{01}_{10}    ddw^{01}_{10} .... │  
        //                │ ddw^{10}_{00}   ddw^{10}_{01}    ddw^{11}_{00}    ddw^{11}_{00} .... │  
        //    f = DW(u) = │ ddw^{10}_{10}   ddw^{10}_{11}    ddw^{11}_{10}    ddw^{11}_{10} .... │  
        //                │ ddw^{20}_{00}   ddw^{20}_{01}    ddw^{21}_{00}    ddw^{21}_{00} .... │  
        //                │ ddw^{20}_{10}   ddw^{20}_{11}    ddw^{21}_{10}    ddw^{21}_{10} .... │  
        //                │        :              :                :                  :          │
        //                │        :              :                :                  :          │
        //                └                                                                      ┘
        //
        //         In addition to being a "flattened" matrix, it is also stored using a sparse
        //         format, using an array of "triplets" which store:
        //             (first index, second index, value)
        //         of the matrix. For instance, the following matrix
        //         
        //                ┌                              ┐
        //                │ 0     val1  0     0     0    │  
        //                │ 0     0     0     0     0    │  
        //                │ 0     0     0     0     0    │  
        //                │ val2  0     0     0     0    │  
        //                │ 0     0     0     0     val3 │  
        //                └                              ┘
        //
        //        can be stored using the following list of triplets:
        //                (0,1,val1)
        //                (3,0,val2)
        //                (4,4,val3)
        //


        Eigen::SparseMatrix<double> ret(size(),size()); // < The Eigen SparseMatrix class manages sparse matrices for us
        std::vector<Eigen::Triplet<double>> coeffs;     // < This is the list of (i,j,val) triplets

        // Calculate the contribution to the stiffness matrix from the
        // CST elements by iterating over the CSTs vector.
        for (int e = 0; e < CSTs.size(); e++)
        {
            // Ths part is identical to that for W and DW.
            auto & id = CSTs[e].getid();
            std::array<Eigen::Vector2d, 3> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);

            // Calculate the LOCAL stiffness matrix for the
            // current CST element.
            std::array<std::array<Eigen::Matrix2d, 3>, 3> ddw = CSTs[e].DDW(in);
            
            // Add each component of the local stiffness matrix to the global
            // stiffness matrix by adding 
            for (int m = 0; m < 3; m++)
                for (int n = 0; n < 3; n++)
                    for (int i = 0; i < 2; i++)
                        for (int j = 0; j < 2; j++)
                        {
                            Eigen::Triplet<double> coeff(2*id[m] + i,     // first index
                                                         2*id[n] + j,     // second inde
                                                         ddw[m][n](i,j)); // value
                            // add this entry to the list of entreis.
                            coeffs.push_back(coeff);
                        }
        }
        

        //  
        // TODO: Following the above, implement for Q4, LST, and Q9 elements.
        //


        // ... your code goes here ...



        //
        // This creates a sparse matrix from the list of coefficient
        // triplest that we have created...
        //
        ret.setFromTriplets(coeffs.begin(), coeffs.end());

        // ...and returns it.
        return ret;
    }

    // nubmer of nodes
    inline const int size()
    {
        return 2 * Points.size();
    }

    inline const int nPoints()
    {
        return 2 * Points.size();
    }
    
    // Total number of elements in this mesh
    inline const int nElements()
    {
        return CSTs.size() + Q4s.size() + LSTs.size() + Q9s.size();
    }

    // Number of nodes, if not accounting for the fact that
    // some are shared between elements.
    inline const int nElementNodes()
    {
        return 4*CSTs.size() + 5*Q4s.size() + 7*LSTs.size() + 10*Q9s.size();
    }

    std::vector<Eigen::Vector2d> Points;
    std::vector<Element::CST<MODEL>> CSTs;
    std::vector<Element::Q4<MODEL>> Q4s;
    std::vector<Element::LST<MODEL>> LSTs;
    std::vector<Element::Q9<MODEL>> Q9s;

};

}
