#ifndef MESH_TEST_H
#define MESH_TEST_H
#include "Util/Exception.H"
#include "Mesh/Mesh.H"

namespace Mesh
{
template<class MESH>
class Test
{
public:
    static void Derivative(std::string meshfile)
    {
        MESH mesh(meshfile);

        double small = 1E-8;
        double tolerance = 1E-4;

        int N = mesh.size();

        double energysum = 0.0;

        for (int iter = 0; iter < 10; iter++)
        {
            Eigen::VectorXd u = Eigen::VectorXd::Random(N);
            
            Eigen::VectorXd dw_exact = mesh.DW(u);
            for (int n = 0; n < N; n++)
            {
                Eigen::VectorXd du = Eigen::VectorXd::Zero(N);
                du(n) = small;

                Eigen::VectorXd unew = u + du;

                energysum += mesh.W(unew);

                double dw_numerical = (mesh.W(unew) - mesh.W(u)) / small;

                double diff = dw_numerical-dw_exact(n);
                double norm = std::max(dw_numerical-dw_exact(n), 1.0);
                if (fabs(diff/norm) > tolerance)
                {
                    std::cout << "[Mesh::Test::Derivative] " << dw_numerical << ", " << dw_exact(n) << std::endl;
                    throw Util::Exception::UnitTest("numerical DW != exact DW");
                }
            }

            Eigen::MatrixXd ddw_exact = mesh.DDW(u);
            for (int n = 0; n < N; n++)
            {
                Eigen::VectorXd du = Eigen::VectorXd::Zero(N);
                du(n) = small;

                Eigen::VectorXd unew = u + du;

                Eigen::VectorXd ddw_numerical = (mesh.DW(unew) - mesh.DW(u)) / small;

                for (int m = 0; m < N; m++)
                {
                    double diff = ddw_numerical(m)-ddw_exact(n,m);
                    double norm = std::max(ddw_numerical(m)+ddw_exact(n,m), 1.0);
                    if (fabs(diff/norm) > tolerance)
                    {
                        std::cout << "[Mesh::Test::Derivative] " << ddw_numerical(m) << ", " << ddw_exact(n,m) << std::endl;
                        throw Util::Exception::UnitTest("numerical DDW != exact DDW");
                    }
                }
            }
        }

        if (fabs(energysum) < tolerance)
        {
            throw Util::Exception::UnitTest("energy returned is always zero");
        }

    }
    
};
}

#endif
