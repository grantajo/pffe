#ifndef ELEMENT_TEST_H
#define ELEMENT_TEST_H
#include "Util/Exception.H"
#include "Element/Element.H"
namespace Element
{
//
// This is a suite of tesing functions to determine whether Elements are
// coded properly.
//
template<class ELEMENT>
class Test
{
private:

    //
    // [function Element::Test<ELEMENT>::Default()]
    //
    // Static utility function (not for general use)
    // to generate a "default" element, where the 
    // embedded position is the same as the template position.
    //
    static ELEMENT Default()
    {
        ELEMENT tmp;
        std::array<Set::Vector,ELEMENT::_N> X0;
        for (int n = 0 ; n < ELEMENT::_N; n++)
        {
            X0[n] = Set::Vector::Zero();
            for (int d = 0; d < ELEMENT::_D; d++)
                X0[n](d) = tmp.getY0()[n](d);
        }
        return ELEMENT(X0);
    }

    //
    // [function RandomPointInElement]
    //
    // Static utility function (not for general use)
    // that produces a random vector in template space that is
    // guaranteed to be within the element, as long as the
    // element is convex.
    //
    static Set::Vector RandomPointInElement(ELEMENT &elem)
    {
        double sum = 0.0;
        const int N = ELEMENT::_N;
        std::array<double,N> lamb;
        for (int n = 0; n < N; n++) lamb[n] = (float)rand();
        for (int n = 0; n < N; n++) sum += lamb[n];
        for (int n = 0; n < N; n++) lamb[n] /= sum;
        auto Y0 = elem.getY0();
        Set::Vector Y = Set::Vector::Zero();
        for (int n = 0; n < N; n++) Y += lamb[n] * Y0[n];
        return Y;
    }    

    //
    // [function Random]
    //
    // Static utility function (not for general use)
    // to generate a random element with randomized
    // embedded position.
    // Element is checked to make sure the embedding is
    // nonsingular, however, there are no guarantees that
    // the points are well-ordered in embedding space.
    //
    static ELEMENT Random()
    {
        ELEMENT tmp;
        for (int iter = 0; iter < 100 ; iter++)
        {
            std::array<Eigen::Vector2d,ELEMENT::_N> X0;
            for (int n = 0 ; n < ELEMENT::_N; n++)
                X0[n] = Eigen::Vector2d::Random();
            ELEMENT elem(X0);

            bool ok = true;

            for (int iter2 = 0; iter2 < 10; iter2++)
            {
                Eigen::Vector2d Ytest = RandomPointInElement(tmp);
                if (fabs(elem.Jacobian(Ytest).determinant()) < 1E-2)
                    ok = false;
            }

            if (ok) return elem;
            else continue;
        }
        throw Util::Exception::UnitTest("Could not create non-singular element");
    }

public:

    //
    // [function Dirac]
    //
    // Test the shape function property that
    //
    // eta_i(y_j) = delta_{ij}
    //
    // where y_j are the element points and eta
    // are the shape functions.
    //
    static void Dirac()
    {
        ELEMENT elem = Default();
        auto Y0 = elem.getY0();
        int N = Y0.size();

        for (int n = 0; n < N; n++)
        {
            for (int m = 0; m < N; m++)
            {
                if (n==m && fabs(elem.eta(Y0[n])[m] - 1.0) > 1E-8)
                    throw Util::Exception::UnitTest(
                        "Dirac test failed: n="+std::to_string(n) +
                        " m=" + std::to_string(m) + 
                        "eta^m(y^n) = " + std::to_string(elem.eta(Y0[n])[m]) + " != 1.0");
                if (n!=m && fabs(elem.eta(Y0[n])[m] - 0.0) > 1E-8)
                    throw Util::Exception::UnitTest(
                        "Dirac test failed: n="+std::to_string(n) +
                        " m=" + std::to_string(m) + 
                        " eta^m(y^n) = " + std::to_string(elem.eta(Y0[n])[m]) + " != 0.0" );
            }
        }
     }

    //
    // [function SumToUnity]
    //
    // This is a test to guarantee that the shape functions all
    // sum to one at any random point within the element.
    //
    static void SumToUnity()
    {
        ELEMENT elem = Default();
        int N = ELEMENT::_N;

        for (int iter = 0; iter < 100; iter++)
        {
            Set::Vector Y = RandomPointInElement(elem);
            double sum = 0.0;
            for (int n = 0; n < N; n++) sum += elem.eta(Y)[n];
            if (fabs(sum-1.0) > 1E-8)
            {
                throw Util::Exception::UnitTest("SumToUnity test failed: sum=" + std::to_string(sum) +
                                  " for Y = [" + std::to_string(Y(0)) + ", " + std::to_string(Y(1)) + "]");
            }
        }
     }


    static void Quadrature()
    {
        ELEMENT elem = Default();
        const int N = ELEMENT::_N;
        const int Q = ELEMENT::_Q;

        double Aexact = elem.getArea();
        std::array<Set::Vector,Q> qp = elem.getQp();
        std::array<double,Q> qw = elem.getQw();

        double Anumerical = 0.0;
        for (int q = 0; q < Q; q++)
            for (int n = 0; n < N; n++)
            {
                Anumerical += qw[q] * elem.eta(qp[q])[n];
            }

        if (fabs(Aexact) < 1E-8)
            throw Util::Exception::UnitTest("Exact area of element is zero");

        if (fabs(Anumerical) < 1E-8)
            throw Util::Exception::UnitTest("Got zero computed area");

        if (fabs(Aexact-Anumerical) > 1E-8)
            throw Util::Exception::UnitTest("Quadrature error: Aexact = " + std::to_string(Aexact) +
                              ", Anumerical = " + std::to_string(Anumerical));
    }


    //
    // [function EtaDerivative]
    //
    // This function checks to make sure the eta derivatives
    // have been computed properly.
    //
    static void EtaDerivative()
    {
        ELEMENT elem = Default();
        int N = ELEMENT::_N;
        int D = ELEMENT::_D;
    
        double small = 1E-8;
        double tolerance = 1E-4;

        for (int iter = 0; iter < 100; iter++)
        {
            Set::Vector Y = RandomPointInElement(elem);

            auto deta_exact = elem.Deta(Y);
            for (int d = 0 ; d < D; d++)
            {
                Set::Vector dy = Set::Vector::Zero();
                dy(d) = small;

                for (int n = 0; n < N; n++)
                {
                    double exact = deta_exact[n](d);
                
                    double numerical = (elem.eta(Y+dy)[n] - elem.eta(Y)[n])/small;

                    double diff = exact-numerical;
                    double norm = std::max(exact+numerical,1.0);
                    if (fabs(diff/(norm+tolerance)) > tolerance)
                        throw Util::Exception::UnitTest("Element::Test::EtaDerivative failed for node " 
                                                        + std::to_string(n) + ", error = "
                                                        + std::to_string(fabs(diff)));
                }
            }
        }
    }


    //
    // [function EtaDerivative]
    //
    // This function checks to make sure the eta derivatives
    // have been computed properly.
    //
    static void IsoparametricCheck()
    {
        ELEMENT standard_element = Default();
        int N = ELEMENT::_N;
    
        Set::Scalar total_diff = 0.0;
        
        for (int iter = 0; iter < 10; iter++)
        {
            ELEMENT random_element = Random();
            for (int p = 0; p < 10; p++)
            {
                std::array<Set::Vector,ELEMENT::_N> u;
                for (int n = 0; n < N ; n++) u[n] = Set::Vector::Random();

                total_diff += fabs(standard_element.W(u) - random_element.W(u));
            }
        }

        if (total_diff < 1E-8)
            throw Util::Exception::UnitTest("Element::Test::IsoparametricCheck failed - energy does not change if embedding changes.");
    }
    
    static void EnergyDerivative()
    {
        double small = 1E-8;
        double tolerance = 1E-3;

        double energysum = 0.0;

        int N = ELEMENT::_N;

        for (int iter = 0; iter < 10; iter++)
        {
            ELEMENT elem = Random();

            std::array<Set::Vector,ELEMENT::_N> u;
            for (int n = 0; n < N ; n++) u[n] = Set::Vector::Random();
            
            double W_exact = elem.W(u);
            energysum += fabs(W_exact);

            auto DW_exact = elem.DW(u);
            for (int n = 0; n < N; n++)
            {
                for (int i = 0 ; i < 2; i++)
                {
                    std::array<Set::Vector,ELEMENT::_N> du;
                    for (int n = 0; n < N ; n++) du[n] = Set::Vector::Zero();
                    du[n](i) = small;

                    std::array<Set::Vector,ELEMENT::_N> u_plus_du;
                    for (int n = 0; n < N ; n++) u_plus_du[n] = u[n] + du[n];
                
                    double DW_numerical_ni = (elem.W(u_plus_du) - elem.W(u)) / small;

                    if (fabs((DW_numerical_ni - DW_exact[n](i))/(DW_numerical_ni + DW_exact[n](i))) > tolerance)
                    {
                        std::cout << "[Element::Test::EnergyDerivative] n = " << n << std::endl;
                        std::cout << "[Element::Test::EnergyDerivative] i = " << i << std::endl;
                        std::cout << "[Element::Test::EnergyDerivative] DW_numerical_ni = "
                                  << DW_numerical_ni << std::endl;
                        std::cout << "[Element::Test::EnergyDerivative] DW_exact[n](i) = "
                                  << DW_exact[n](i) << std::endl;
                        throw Util::Exception::UnitTest("[Element::Test::EnergyDerivative] failed (DW)");
                    }
                }
            }

            auto DDW_exact = elem.DDW(u);
            for (int n = 0; n < N; n++)
            {
                for (int i = 0 ; i < 2; i++)
                {
                    std::array<Set::Vector,ELEMENT::_N> du;
                    for (int m = 0; m < N ; m++) du[m] = Set::Vector::Zero();
                    du[n](i) = small;

                    std::array<Set::Vector,ELEMENT::_N> u_plus_du;
                    for (int m = 0; m < N ; m++) u_plus_du[m] = u[m] + du[m];
                
                    std::array<Set::Vector,ELEMENT::_N> DDW_numerical_ni;
                    for (int m = 0; m < N; m++) DDW_numerical_ni[m] = (elem.DW(u_plus_du)[m] - elem.DW(u)[m]) / small;

                    for (int m = 0; m < N; m++)
                    {
                        for (int j = 0 ; j < 2; j++)
                        {
                            if (fabs((DDW_numerical_ni[m](j) - DDW_exact[n][m](i,j)) /
                                     (DDW_numerical_ni[m](j) + DDW_exact[n][m](i,j))) > tolerance)
                            {
                                std::cout << "[Element::Test::EnergyDerivative] n = " << n << std::endl;
                                std::cout << "[Element::Test::EnergyDerivative] m = " << m << std::endl;
                                std::cout << "[Element::Test::EnergyDerivative] i = " << i << std::endl;
                                std::cout << "[Element::Test::EnergyDerivative] j = " << j << std::endl;
                                std::cout << "[Element::Test::EnergyDerivative] DDW_numerical_ni[m](j) = " << DDW_numerical_ni[m](j) << std::endl;
                                std::cout << "[Element::Test::EnergyDerivative] DDW_exact[n][m](i,j) = " << DDW_exact[n][m](i,j) << std::endl;
                                throw Util::Exception::UnitTest("Element::Test::EnergyDerivative failed (DDW)");
                            }
                        }
                    }
                }
            }
        }
        if (fabs(energysum) < 1E-8)
        {
            throw Util::Exception::UnitTest("Energy returned by element is always zero.");
        }
    }
};
}
#endif 
