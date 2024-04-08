#ifndef SET_SET_H
#define SET_SET_H

namespace Set
{

using Scalar   = double;
using Vector   = Eigen::Matrix<double,DIM,1>;
using Vector2d = Eigen::Vector2d;
using Vector3d = Eigen::Vector3d;
using Matrix   = Eigen::Matrix<double,DIM,DIM>;
using Matrix2d = Eigen::Matrix2d;
using Matrix3d = Eigen::Matrix3d;


class Tensor4
{
public:
    Scalar & operator () (const int i, const int j, const int k, const int l) 
    {
        return data[i][j][k][l];
    }
    Scalar operator () (const int i, const int j, const int k, const int l) const
    {
        return data[i][j][k][l];
    }
    static Tensor4 Zero()
    {
        Tensor4 ret;
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                for (int k = 0; k < DIM; k++)
                    for (int l = 0; l < DIM; l++)
                        ret.data[i][j][k][l] = 0;
        return ret;
    }


private:
    std::array<std::array<std::array<std::array<double,DIM>,DIM>,DIM>,DIM> data;

};
}
     


#endif
