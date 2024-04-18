#ifndef ELEMENT_Q9_H
#define ELEMENT_Q9_H

#include "Element.H"

namespace Element
{
template<class MODEL>
class Q9 : public Element<2,9,9,MODEL>
{
public:
    using Element<2,9,9,MODEL>::Element;

    std::array<double,9> eta(Set::Vector Y) override
    {
        // Create shape functions
        std::array<double,9> ret;

        ret[0] = 0.25*(Y(0)*Y(0)*Y(1)*Y(1) - Y(0)*Y(0)*Y(1) - Y(0)*Y(1)*Y(1) + Y(0)*Y(1));
        ret[1] = 0.5*(Y(1)*Y(1) - Y(1) - Y(0)*Y(0)*Y(1)*Y(1) + Y(0)*Y(0)*Y(1));
        ret[2] = 0.25*(Y(0)*Y(0)*Y(1)*Y(1) - Y(0)*Y(0)*Y(1) + Y(0)*Y(1)*Y(1) - Y(0)*Y(1));
        ret[3] = 0.5*(Y(0) - Y(0)*Y(1)*Y(1) + Y(0)*Y(0) - Y(0)*Y(0)*Y(1)*Y(1));
        ret[4] = 0.25*(Y(0)*Y(0)*Y(1)*Y(1) + Y(0)*Y(0)*Y(1) + Y(0)*Y(1)*Y(1) + Y(0)*Y(1));
        ret[5] = 0.5*(Y(1)*Y(1) + Y(1) - Y(0)*Y(0)*Y(1)*Y(1) - Y(0)*Y(0)*Y(1));
        ret[6] = 0.25*(Y(0)*Y(0)*Y(1)*Y(1) + Y(0)*Y(0)*Y(1) - Y(0)*Y(1)*Y(1) - Y(0)*Y(1));
        ret[7] = 0.5*(Y(0)*Y(0) - Y(0)*Y(0)*Y(1)*Y(1) - Y(0) + Y(0)*Y(1)*Y(1));
        ret[8] = 1 - Y(1)*Y(1) - Y(0)*Y(0) + Y(0)*Y(0)*Y(1)*Y(1);

        return ret;
    }

    std::array<Set::Vector,9> Deta(Set::Vector Y) override
    {
        // Create derivatives of shape functions
        std::array<Eigen::Vector2d,9> Dret;

        Dret[0](0) = 0.25*(2*Y(0)*Y(1)*Y(1) - 2*Y(0)*Y(1) - Y(1)*Y(1) + Y(1)); 
        Dret[0](1) = 0.25*(2*Y(0)*Y(0)*Y(1) - Y(0)*Y(0) - 2*Y(0)*Y(1) + Y(0));

        Dret[1](0) = 0.5*(-2*Y(0)*Y(1)*Y(1) + 2*Y(0)*Y(1));
        Dret[1](1) = 0.5*(2*Y(1) -1 - 2*Y(0)*Y(0)*Y(1) + Y(0)*Y(0));

        Dret[2](0) = 0.25*(2*Y(0)*Y(1)*Y(1) - 2*Y(0)*Y(1) + Y(1)*Y(1) - Y(1));
        Dret[2](1) = 0.25*(2*Y(0)*Y(0)*Y(1) - Y(0)*Y(0) + 2*Y(0)*Y(1) - Y(0));

        Dret[3](0) = 0.5*(1 - Y(1)*Y(1) + 2*Y(0) - 2*Y(0)*Y(1)*Y(1));
        Dret[3](1) = 0.5*(-2*Y(0)*Y(1) - 2*Y(0)*Y(0)*Y(1));

        Dret[4](0) = 0.25*(2*Y(0)*Y(1)*Y(1) + 2*Y(0)*Y(1) + Y(1)*Y(1) + Y(1));
        Dret[4](1) = 0.25*(2*Y(0)*Y(0)*Y(1) + Y(0)*Y(0) + 2*Y(0)*Y(1) + Y(0));
        
        Dret[5](0) = 0.5*(-2*Y(0)*Y(1)*Y(1) - 2*Y(0)*Y(1));
        Dret[5](1) = 0.5*(2*Y(1) + 1 - 2*Y(0)*Y(0)*Y(1) - Y(0)*Y(0));

        Dret[6](0) = 0.25*(2*Y(0)*Y(1)*Y(1) + 2*Y(0)*Y(1) - Y(1)*Y(1) - Y(1));
        Dret[6](1) = 0.25*(2*Y(0)*Y(0)*Y(1) + Y(0)*Y(0) - 2*Y(0)*Y(1) - Y(0));

        Dret[7](0) = 0.5*(2*Y(0) - 2*Y(0)*Y(1)*Y(1) - 1 + Y(1)*Y(1));
        Dret[7](1) = 0.5*(-2*Y(0)*Y(0)*Y(1) + 2*Y(0)*Y(1));

        Dret[8](0) = -2*Y(0) + 2*Y(0)*Y(1)*Y(1);
        Dret[8](1) = -2*Y(1) + 2*Y(0)*Y(0)*Y(1);

        return Dret;
    }

public:
    virtual const std::array<Set::Vector,9> & getY0() {return Y0;}
    virtual const std::array<Set::Vector,9> & getQp() {return Qp;}
    virtual const std::array<double,9> & getQw() {return Qw;}
    virtual const double getArea() {
        // Return area of template
        return 4.0;
    }

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,9> Y0 = {
        Set::Vector(-1.0,-1.0),
        Set::Vector( 0.0,-1.0),
        Set::Vector( 1.0,-1.0),
        Set::Vector( 1.0, 0.0),
        Set::Vector( 1.0, 1.0),
        Set::Vector( 0.0, 1.0),
        Set::Vector(-1.0, 1.0),
        Set::Vector(-1.0, 0.0),
        Set::Vector( 0.0, 0.0),
    };

    static constexpr double sq3 = sqrt(3.);
    static constexpr double sq5 = sqrt(5.);
    const std::array<Set::Vector,9> Qp = {
        // Create quadrature points
        Eigen::Vector2d(0.,0.),
        Eigen::Vector2d(sq3/sq5,0.),
        Eigen::Vector2d(-sq3/sq5,0.),
        Eigen::Vector2d(0.,sq3/sq5),
        Eigen::Vector2d(0.,-sq3/sq5),
        Eigen::Vector2d(sq3/sq5,sq3/sq5),
        Eigen::Vector2d(sq3/sq5,-sq3/sq5),
        Eigen::Vector2d(-sq3/sq5,sq3/sq5),
        Eigen::Vector2d(-sq3/sq5,-sq3/sq5)  
    };

    const std::array<double,9> Qw = {
        // Create quadrature weights
        64./81., 40./81., 40./81., 40./81., 40./81., 25./81., 25./81., 25./81., 25./81.
    };

public:
    
    static const int vtkid = 28;
};

}
#endif
