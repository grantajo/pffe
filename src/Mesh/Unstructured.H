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
        // TODO: Complete the Unstructured Mesh file reader.
        //
        // Input: "vtkfile" is the name of a file containing a mesh in
        //        human-readable ASCII vtk format.
        //
        // Output: None - this is a constructor
        //

        // Check to see if the file exists - if not, then exit.
        if (!std::filesystem::exists(vtkfile))
            throw std::runtime_error("Could not find file " + vtkfile);
        std::ifstream in(vtkfile);

        // Create string vectors of points, cells, and cell types.
        std::vector<std::string> str_points;
        std::vector<std::string> str_cells;
        std::vector<std::string> str_cell_types;
        enum status {none, point, cell, cell_type};
        status mystatus = none;

        // Read in all of the data from the mesh file and store in these vectors.
        std::string line;
        while (getline(in,line))
        {
            if (line=="") {mystatus = none; continue;}
            std::stringstream ss(line);
            std::string first;
            ss >> first;
            if (first=="POINTS") {mystatus = point; continue;}
            if (first=="CELLS") {mystatus = cell; continue;}
            if (first=="CELL_TYPES") {mystatus = cell_type; continue;}
            
            if (mystatus == point) str_points.push_back(line);
            if (mystatus == cell) str_cells.push_back(line);
            if (mystatus == cell_type) str_cell_types.push_back(line);
        }
        
        // Make sure we did not miscount.
        assert(str_cells.size() == str_cell_types.size());
        
        // Read in all of the vertices and store in the
        // "Points" array
        for (int n = 0; n < str_points.size(); n++)
        {
            Set::Vector Point;
            std::stringstream in(str_points[n]);
            for (int d = 0; d < DIM; d++) in >> Point(d);
            Points.push_back(Point);
        }
        
        // Read in the CST elements
        for (int n = 0; n < str_cells.size(); n++) // Iterate over all cells
        {
            // CST elements have VTK ID = 5
            // (https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html)
            //
            // If this cell is of type 5, then read in and add to the CST array
            if (str_cell_types[n] == std::to_string(5)) 
            {
                // Convert the list of nodes to a stringstream
                std::stringstream in(str_cells[n]); 
                // This is an array of node IDS. It has length 3 because
                // the CST has 3 nodes.
                std::array<int,3> id;
                // We don't care about the first number but we have to read it anyway, so
                // call it "dummy"
                int dummy; 
                // Store the IDs in order - note that this order will change
                // depending on the node order we use vs that of the VTK
                // (seee the documentation!)
                in >> dummy >> id[0] >> id[1] >> id[2]; 
                // Add this new CST element to the list.
                CSTs.push_back(Element::CST<MODEL>(Points,id));
            }
        }

        // LST
        for (int n = 0; n < str_cells.size(); n++) // Iterate over all cells
        {
            // LST elements have VTK ID = 22
            // (https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html)
            //
            // If this cell is of type 22, then read in and add to the LST array
            if (str_cell_types[n] == std::to_string(22)) 
            {
                // Convert the list of nodes to a stringstream
                std::stringstream in(str_cells[n]); 
                // This is an array of node IDS. It has length 6 because
                // the LST has 6 nodes.
                std::array<int,6> id;
                // We don't care about the first number but we have to read it anyway, so
                // call it "dummy"
                int dummy; 
                // Store the IDs in order - note that this order will change
                // depending on the node order we use vs that of the VTK
                // (seee the documentation!)
                in >> dummy >> id[0] >> id[2] >> id[4] >> id[1] >> id[3] >> id[5]; 
                // Add this new LST element to the list.
                LSTs.push_back(Element::LST<MODEL>(Points,id));
            }
        }
        
        // Q4
        for (int n = 0; n < str_cells.size(); n++) // Iterate over all cells
        {
            // Q4 elements have VTK ID = 9
            // (https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html)
            //
            // If this cell is of type 9, then read in and add to the CST array
            if (str_cell_types[n] == std::to_string(9)) 
            {
                // Convert the list of nodes to a stringstream
                std::stringstream in(str_cells[n]); 
                // This is an array of node IDS. It has length 4 because
                // the Q4 has 4 nodes.
                std::array<int,4> id;
                // We don't care about the first number but we have to read it anyway, so
                // call it "dummy"
                int dummy; 
                // Store the IDs in order - note that this order will change
                // depending on the node order we use vs that of the VTK
                // (seee the documentation!)
                in >> dummy >> id[0] >> id[1] >> id[3] >> id[2]; 
                // Add this new Q4 element to the list.
                Q4s.push_back(Element::Q4<MODEL>(Points,id));
            }
        }

        for (int n = 0; n < str_cells.size(); n++) // Iterate over all cells
        {
            // Q9 elements have VTK ID = 28
            // (https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html)
            //
            // If this cell is of type 5, then read in and add to the CST array
            if (str_cell_types[n] == std::to_string(28)) 
            {
                // Convert the list of nodes to a stringstream
                std::stringstream in(str_cells[n]); 
                // This is an array of node IDS. It has length 9 because
                // the Q9 has 9 nodes.
                std::array<int,9> id;
                // We don't care about the first number but we have to read it anyway, so
                // call it "dummy"
                int dummy; 
                // Store the IDs in order - note that this order will change
                // depending on the node order we use vs that of the VTK
                // (see the documentation!)
                in >> dummy >> id[0] >> id[2] >> id[4] >> id[6] >> id[1] >> id[3] >> id[5] >> id[7] >> id[8]; 
                // Add this new Q9 element to the list.
                Q9s.push_back(Element::Q9<MODEL>(Points,id));
            }
        }

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
        //
        // Todo: Complete the mesh writer.
        //       This function writes the mesh data to VTK format.
        //
        // Input: "vtkfile" - file name to write to.

        // Open the file and write the header information.
        std::ofstream out(vtkfile);
        out << "# vtk DataFile Version 2.0" << std::endl;
        out << "created by fem" << std::endl;
        out << "ASCII" << std::endl;
        out << "DATASET UNSTRUCTURED_GRID" << std::endl;

        // Write out all of the points.
        out << "POINTS " << Points.size() << " double" << std::endl;
        for (int n = 0; n < Points.size(); n++)
            out << Points[n].transpose() << " 0.0" << std::endl;
        
        // Indicate how many total points we are about to write.
        out << std::endl;
        out << "CELLS " << nElements() << " " << nElementNodes() << std::endl;

        // Start with the CSTs. 
        for (int e = 0; e < CSTs.size(); e++)
        {
            // Get the array of node IDs.
            // The CST has three nodes, so the length of the array is 3.
            const std::array<int,3> & id = CSTs[e].getid();
            // Output the node ID order in VTK ordering
            out << "3 " << id[0] << " " << id[1] << " " << id[2] << std::endl;
        }

        // Q4
        for (int e = 0; e < Q4s.size(); e++)
        {
            // Get the array of node IDs.
            // The Q4 has four nodes, so the length of the array is 4.
            const std::array<int,4> & id = Q4s[e].getid();
            // Output the node ID order in VTK ordering
            out << "4 " << id[0] << " " << id[1] << " " << id[3] << " " << id[2] << std::endl;
        }

        // LST
        for (int e = 0; e < LSTs.size(); e++)
        {
            // Get the array of node IDs.
            // The LST has six nodes, so the length of the array is 6.
            const std::array<int,6> & id = LSTs[e].getid();
            // Output the node ID order in VTK ordering
            out << "6 " << id[0] << " " << id[2] << " " << id[4] << " " << id[1] << " " << id[3] << " " << id[5] << std::endl;
        }

        // Q9
        for (int e = 0; e < Q9s.size(); e++)
        {
            // Get the array of node IDs.
            // The Q9 has nine nodes, so the length of the array is 9.
            const std::array<int,9> & id = Q9s[e].getid();
            // Output the node ID order in VTK ordering
            out << "9 " << id[0] << " " << id[2] << " " << id[4] << " " << id[6] << " " << id[1] << " " << id[3] << " " << id[5] << " " << id[7] << " " << id[8] << std::endl;
        }

        // Now we need to specify what kind of element each of the above rows corresponds to.
        out << std::endl;
        out << "CELL_TYPES " << nElements() << std::endl;

        // CSTs are are VTK type 5, so write out one "5" for each CST in the array.
        for (int e = 0; e < CSTs.size(); e++) {
            out << "5" << std::endl;
        }

        // Q4s are VTK type 9, so write out a '9' for each Q4 in the array
        for (int e = 0; e < Q4s.size(); e++) {
            out << "9" << std::endl;
        }

        // LSTs are are VTK type 22, so write out one "22" for each CST in the array.
        for (int e = 0; e < LSTs.size(); e++) {
            out << "22" << std::endl;
        }

        // Q9s are VTK type 28, so write out a '28' for each Q4 in the array
        for (int e = 0; e < Q9s.size(); e++) {
            out << "28" << std::endl;
        }

        out << std::endl;

        out << "POINT_DATA " << Points.size() << std::endl;
        out << "VECTORS disp float" << std::endl;
        
        for (int e = 0; e < Points.size(); e++) {
            out << disp(2*e) << " " << disp(2*e + 1) << " " << 0.0 << std::endl;
        }

        //  POINT_DATA  [[ total number of vertices ]] 
        //  VECTORS disp float
        //  [x displacement of node 0] [y displacement of node 0] 0.0
        //  [x displacement of node 1] [y displacement of node 2] 0.0

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

        // Loop over all Q4 element:
        for (int e = 0; e < Q4s.size(); e++)
        {
            // Get the vector of IDs from CST element e
            auto & id = Q4s[e].getid();
            
            // Construct an array of vectors ("in") to feed into
            // the element's W function.
            std::array<Set::Vector2d, 4> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);

            
            // Add up the energy using the CST element's energy function.
            ret += Q4s[e].W(in);
        }

        // Loop over all LST element:
        for (int e = 0; e < LSTs.size(); e++)
        {
            // Get the vector of IDs from CST element e
            auto & id = LSTs[e].getid();
            
            // Construct an array of vectors ("in") to feed into
            // the element's W function.
            std::array<Set::Vector2d, 6> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            in[4] << u(2*id[4]), u(2*id[4] + 1);
            in[5] << u(2*id[5]), u(2*id[5] + 1);
        
            // Add up the energy using the CST element's energy function.
            ret += LSTs[e].W(in);
        }
        
        // Loop over all Q9 element:

        for (int e = 0; e < Q9s.size(); e++)
        {
            // Get the vector of IDs from CST element e
            auto & id = Q9s[e].getid();
            
            // Construct an array of vectors ("in") to feed into
            // the element's W function.
            std::array<Set::Vector2d, 9> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            in[4] << u(2*id[4]), u(2*id[4] + 1);
            in[5] << u(2*id[5]), u(2*id[5] + 1);
            in[6] << u(2*id[6]), u(2*id[6] + 1);
            in[7] << u(2*id[7]), u(2*id[7] + 1);
            in[8] << u(2*id[8]), u(2*id[8] + 1);
        
            // Add up the energy using the CST element's energy function.
            ret += Q9s[e].W(in);
        }

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
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            // Now we calculate the dw for each element, which
            // is an array of vectors.
            std::array<Set::Vector2d, 3> dw = CSTs[e].DW(in);
            
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
        
        // Loop over Q4s array.
        for (int e = 0; e < Q4s.size(); e++)
        {
            // This part is identical to that for W.
            auto & id = Q4s[e].getid();
            std::array<Set::Vector2d, 4> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);

            // Now we calculate the dw for each element, which
            // is an array of vectors.
            std::array<Set::Vector2d, 4> dw = Q4s[e].DW(in);
            
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
            ret(2*id[3]  ) += dw[3](0);
            ret(2*id[3]+1) += dw[3](1);
        }

        // Loop over LST array.
        for (int e = 0; e < LSTs.size(); e++)
        {
            // This part is identical to that for W.
            auto & id = LSTs[e].getid();
            std::array<Set::Vector2d, 6> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            in[4] << u(2*id[4]), u(2*id[4] + 1);
            in[5] << u(2*id[5]), u(2*id[5] + 1);

            // Now we calculate the dw for each element, which
            // is an array of vectors.
            std::array<Set::Vector2d, 6> dw = LSTs[e].DW(in);
            
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
            ret(2*id[3]  ) += dw[3](0);
            ret(2*id[3]+1) += dw[3](1);
            ret(2*id[4]  ) += dw[4](0);
            ret(2*id[4]+1) += dw[4](1);
            ret(2*id[5]  ) += dw[5](0);
            ret(2*id[5]+1) += dw[5](1);
        }

        // Loop over Q9 array.
        for (int e = 0; e < Q9s.size(); e++)
        {
            // This part is identical to that for W.
            auto & id = Q9s[e].getid();
            std::array<Set::Vector2d, 9> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            in[4] << u(2*id[4]), u(2*id[4] + 1);
            in[5] << u(2*id[5]), u(2*id[5] + 1);
            in[6] << u(2*id[6]), u(2*id[6] + 1);
            in[7] << u(2*id[7]), u(2*id[7] + 1);
            in[8] << u(2*id[8]), u(2*id[8] + 1);
        
            // Now we calculate the dw for each element, which
            // is an array of vectors.
            std::array<Set::Vector2d, 9> dw = Q9s[e].DW(in);
            
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
            ret(2*id[3]  ) += dw[3](0);
            ret(2*id[3]+1) += dw[3](1);
            ret(2*id[4]  ) += dw[4](0);
            ret(2*id[4]+1) += dw[4](1);
            ret(2*id[5]  ) += dw[5](0);
            ret(2*id[5]+1) += dw[5](1);
            ret(2*id[6]  ) += dw[6](0);
            ret(2*id[6]+1) += dw[6](1);
            ret(2*id[7]  ) += dw[7](0);
            ret(2*id[7]+1) += dw[7](1);
            ret(2*id[8]  ) += dw[8](0);
            ret(2*id[8]+1) += dw[8](1);
        }

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

        // Calculate the contribution to the stiffness matrix from the
        // Q4 elements by iterating over the Q4s vector.
        for (int e = 0; e < Q4s.size(); e++)
        {
            // Ths part is identical to that for W and DW.
            auto & id = Q4s[e].getid();
            std::array<Eigen::Vector2d, 4> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);

            // Calculate the LOCAL stiffness matrix for the
            // current Q4 element.
            std::array<std::array<Eigen::Matrix2d, 4>, 4> ddw = Q4s[e].DDW(in);
            
            // Add each component of the local stiffness matrix to the global
            // stiffness matrix by adding 
            for (int m = 0; m < 4; m++)
                for (int n = 0; n < 4; n++)
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

        // Calculate the contribution to the stiffness matrix from the
        // LST elements by iterating over the LSTs vector.
        for (int e = 0; e < LSTs.size(); e++)
        {
            // Ths part is identical to that for W and DW.
            auto & id = LSTs[e].getid();
            std::array<Eigen::Vector2d, 6> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            in[4] << u(2*id[4]), u(2*id[4] + 1);
            in[5] << u(2*id[5]), u(2*id[5] + 1);

            // Calculate the LOCAL stiffness matrix for the
            // current LST element.
            std::array<std::array<Eigen::Matrix2d, 6>, 6> ddw = LSTs[e].DDW(in);
            
            // Add each component of the local stiffness matrix to the global
            // stiffness matrix by adding 
            for (int m = 0; m < 6; m++)
                for (int n = 0; n < 6; n++)
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

        // Calculate the contribution to the stiffness matrix from the
        // Q9 elements by iterating over the Q9s vector.
        for (int e = 0; e < Q9s.size(); e++)
        {
            // Ths part is identical to that for W and DW.
            auto & id = Q9s[e].getid();
            std::array<Eigen::Vector2d, 9> in;
            in[0] << u(2*id[0]), u(2*id[0] + 1);
            in[1] << u(2*id[1]), u(2*id[1] + 1);
            in[2] << u(2*id[2]), u(2*id[2] + 1);
            in[3] << u(2*id[3]), u(2*id[3] + 1);
            in[4] << u(2*id[4]), u(2*id[4] + 1);
            in[5] << u(2*id[5]), u(2*id[5] + 1);
            in[6] << u(2*id[6]), u(2*id[6] + 1);
            in[7] << u(2*id[7]), u(2*id[7] + 1);
            in[8] << u(2*id[8]), u(2*id[8] + 1);

            // Calculate the LOCAL stiffness matrix for the
            // current Q4 element.
            std::array<std::array<Eigen::Matrix2d, 9>, 9> ddw = Q9s[e].DDW(in);
            
            // Add each component of the local stiffness matrix to the global
            // stiffness matrix by adding 
            for (int m = 0; m < 9; m++)
                for (int n = 0; n < 9; n++)
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
