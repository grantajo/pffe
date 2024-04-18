#ifndef ELEMENT_Q4_H
#define ELEMENT_Q4_H

#include "Element.H"

namespace Element
{
template<class MODEL>
class Q4 : public Element<2,4,4,MODEL>
{
public:
    using Element<2,4,4,MODEL>::Element;

    std::array<double,4> eta(Set::Vector Y) override
    {
        // Create shape functions
        std::array<double,4> ret;

        ret[0] = 0.25*(Y(0)*Y(1) - Y(0) - Y(1) + 1);
        ret[1] = -0.25*(Y(0)*Y(1) - Y(0) + Y(1) - 1);
        ret[2] = 0.25*(Y(0)*Y(1) + Y(0) + Y(1) + 1);
        ret[3] = -0.25*(Y(0)*Y(1) + Y(0) - Y(1) - 1);

        return ret;
    }

    std::array<Set::Vector,4> Deta(Set::Vector Y) override
    {
        // Create derivatives of shape functions
        std::array<Eigen::Vector2d,4> Deta;

        Deta[0](0) = 0.25*(Y(1) - 1);
        Deta[0](1) = 0.25*(Y(0) - 1);

        Deta[1](0) = -0.25*(Y(1) - 1);
        Deta[1](1) = -0.25*(Y(0) + 1);

        Deta[2](0) = 0.25*(Y(1) + 1); 
        Deta[2](1) = 0.25*(Y(0) + 1);

        Deta[3](0) = -0.25*(Y(1) + 1);
        Deta[3](1) = -0.25*(Y(0) - 1);

        return Deta;
    }

public:
    virtual const std::array<Set::Vector,4> & getY0() {return Y0;}
    virtual const std::array<Set::Vector,4> & getQp() {return Qp;}
    virtual const std::array<double,4> & getQw() {return Qw;}
    virtual const double getArea() {
        // Return area of template element
        return 4.0;
    }

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,4> Y0 = {
        Set::Vector(-1.0,-1.0),
        Set::Vector(1.0,-1.0),
        Set::Vector(1.0,1.0),
        Set::Vector(-1.0,1.0)
    };

    static constexpr double sq3 = sqrt(3.);
    const std::array<Set::Vector,4> Qp = {
        // Create quadrature points

        Eigen::Vector2d(1.0/sq3,1.0/sq3),
        Eigen::Vector2d(1.0/sq3,-1.0/sq3),
        Eigen::Vector2d(-1.0/sq3,1.0/sq3),
        Eigen::Vector2d(-1.0/sq3,-1.0/sq3)
    };

    const std::array<double,4> Qw = {
        // Create quadrature weights
        1., 1., 1., 1.
    };


public:
    
    static const int vtkid = 5;
};

}
#endif
