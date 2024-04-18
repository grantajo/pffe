#ifndef ELEMENT_CST_H
#define ELEMENT_CST_H

#include "Element.H"

namespace Element
{
template<class MODEL>
class CST : public Element<2,3,1,MODEL>
{
public:
    using Element<2,3,1,MODEL>::Element;

    std::array<double,3> eta(Set::Vector Y) override
    {
        std::array<double,3> ret;

        ret[0] = 1 - Y(0) - Y(1);
        ret[1] = Y(0);
        ret[2] = Y(1);
        return ret;
    }

    std::array<Set::Vector,3> Deta(Set::Vector Y) override
    {
        std::array<Eigen::Vector2d,3> ret;
        
        ret[0](0) = -1;
        ret[0](1) = -1;
        ret[1](0) = 1;
        ret[1](1) = 0;
        ret[2](0) = 0;
        ret[2](1) = 1;

        return ret;
    }

public:
    //
    // These are utility functions 
    //
    virtual const std::array<Set::Vector,3> & getY0() {return Y0;}
    virtual const std::array<Set::Vector,1> & getQp() {return Qp;}
    virtual const std::array<double,1> & getQw() {return Qw;}
    virtual const double getArea() {
        // TODO - copy from previous
        return 0.5;
    }

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,3> Y0 = {
        Set::Vector(0.0,0.0),
        Set::Vector(1.0,0.0),
        Set::Vector(0.0,1.0)
    };

    const std::array<Set::Vector,1> Qp = {
        Eigen::Vector2d(1.0/3.0,1.0/3.0)
    };

    const std::array<double,1> Qw = {
        0.5
    };


public:
    
    static const int vtkid = 5;
};

}
#endif
