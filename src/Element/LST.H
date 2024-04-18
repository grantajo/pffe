#ifndef ELEMENT_LST_H
#define ELEMENT_LST_H

#include "Element.H"

namespace Element
{
template<class MODEL>
class LST : public Element<2,6,4,MODEL>
{
public:
    using Element<2,6,4,MODEL>::Element;

    std::array<double,6> eta(Set::Vector Y) override
    {
        // Create shape functions
        std::array<double,6> ret;

        ret[0] = 1 + 2*Y(0)*Y(0) - 3*Y(0) + 2*Y(1)*Y(1) - 3*Y(1) + 4*Y(0)*Y(1);
        ret[1] = 4 * (1 - Y(0) - Y(1)) * Y(0);
        ret[2] = (2*Y(0) - 1) * Y(0);
        ret[3] = 4 * Y(0) * Y(1);
        ret[4] = (2* Y(1) - 1) * Y(1);
        ret[5] = 4 * (1 - Y(0) - Y(1)) * Y(1);
        
        return ret;
    }

    std::array<Set::Vector,6> Deta(Set::Vector Y) override
    {
        // Create derivatives of shape functions
        std::array<Eigen::Vector2d,6> ret;

        ret[0](0) = 4*Y(0) + 4*Y(1) - 3; 
        ret[0](1) = 4*Y(1) + 4*Y(0) - 3;

        ret[1](0) = 4 - 8*Y(0) - 4*Y(1);
        ret[1](1) = -4*Y(0);

        ret[2](0) = 4*Y(0) - 1;
        ret[2](1) = 0;

        ret[3](0) = 4*Y(1);
        ret[3](1) = 4*Y(0);

        ret[4](0) = 0;
        ret[4](1) = 4*Y(1) - 1;
        
        ret[5](0) = -4*Y(1);
        ret[5](1) = 4 - 8*Y(1) - 4*Y(0);

        return ret;
    }

public:
    virtual const std::array<Set::Vector,6> & getY0() {return Y0;}
    virtual const std::array<Set::Vector,4> & getQp() {return Qp;}
    virtual const std::array<double,4> & getQw() {return Qw;}
    virtual const double getArea() {
        // Return the area of the template element
        return 0.5;
    }

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,6> Y0 = {
        Set::Vector(0.0,0.0),
        Set::Vector(0.5,0.0),
        Set::Vector(1.0,0.0),
        Set::Vector(0.5,0.5),
        Set::Vector(0.0,1.0),
        Set::Vector(0.0,0.5)
    };

    const std::array<Set::Vector,4> Qp = {
        // Create quadrature points
        Eigen::Vector2d(1.0/3.0,1.0/3.0),
        Eigen::Vector2d(2.0/15.0,2.0/15.0),
        Eigen::Vector2d(2.0/15.0,11.0/15.0),
        Eigen::Vector2d(11.0/15.0,2.0/15.0)
    };

    const std::array<double,4> Qw = {
        // Create quadrature weights
        -27.0/96.0, 25.0/96, 25.0/96, 25.0/96
    };


public:
    
    static const int vtkid = 22;

};

}
#endif
