#ifndef MODEL_MODEL_H
#define MODEL_MODEL_H

#include "eigen3/Eigen/Core"

#include "Set/Set.H"

namespace Model
{
template<int D>
class Model
{
public:
    virtual Set::Scalar   W(Set::Matrix) = 0;
    virtual Set::Matrix   DW(Set::Matrix) = 0;
    virtual Set::Tensor4  DDW(Set::Matrix) = 0;
};
}
#endif
