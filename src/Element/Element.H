#ifndef ELEMENT_ELEMENT_H
#define ELEMENT_ELEMENT_H
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include <time.h>
#include "Util/Exception.H"
#include "Model/Model.H"

namespace Element
{
template<int D, int N, int Q, class MODEL>
class Element
{
public:

    Element()
    {}
                
    Element(std::vector<Set::Vector> &a_X0, std::array<int,N> a_id) :
        X0(&a_X0), id(a_id)
    {}

    Element(std::array<Set::Vector,N> &a_X0) 
    {
        X0 = new std::vector<Set::Vector>();
        for (int n = 0 ; n < N ; n++)
        {
            X0->push_back(a_X0[n]);
            id[n] = n;
        }
    }

    double W(std::array<Set::Vector,N> &u)
    {
        Set::Scalar ret = 0.0;
        
        const std::array<Set::Scalar,Q> &Qw = getQw();
        const std::array<Set::Vector,Q> &Qp = getQp();
        Set::Matrix beta = Set::Matrix::Zero();

        //
        // TODO: Calculate the energy of the element as a function of
        //       nodal displacement.
        //       (See Problem 1)
        //
        // Input: the argument "u" is an array of vectors.
        //        The following expression:
        //            u[n](i)
        //        returns the ith component of the displacement vector
        //        at node n.
        //
        // Output: the total energy.
        //
        // Hint 1: use the following to access quadrature weights
        //         and points.
        //         (you will need to use quadrature to evaluate the
        //         integral)
        //
        // Hint 2: To evaluate the elastic model at a point, you can
        //         use the following expression:
        //            model[q].W(gradu)
        //         where q is a quadrature point index, gradu is a
        //         matrix containing the displacement gradient, and the
        //         result is a scalar value for the energy.
        //         (you can do the same thing for DW and DDW)
        //

        for (int q = 0; q < Q; q++) {
            Set::Matrix J = Set::Matrix::Zero();
            J = Jacobian(Qp[q]);
            Set::Matrix Jinv = Set::Matrix::Zero();
            Jinv = J.inverse();
            Set::Scalar Jdet = 0.0;
            Jdet = J.determinant();
            beta = Set::Matrix::Zero();
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < D; j++) {
                    for (int n = 0; n < N; n++) {
                        for (int k = 0; k < D; k++) {
                            beta(i,j) += u[n](i)*Deta(Qp[q])[n](k) * Jinv(k,j);
                        }
                    }
                }
            }
            ret += Qw[q] * model[q].W(beta) * Jdet;
        }

        return ret;
    }

    std::array<Set::Vector,N> DW(std::array<Set::Vector,N> &u)
    {
        std::array<Set::Vector,N> ret;
        for (int n = 0; n < N; n++) ret[n] = Set::Vector::Zero();
        const std::array<Set::Scalar,Q> &Qw = getQw();
        const std::array<Set::Vector,Q> &Qp = getQp();        
        Set::Matrix beta = Set::Matrix::Zero();

        //
        // TODO: Calculate the derivative of energy of the element as a function of
        //       nodal displacement.
        //       (See Problem 1 part a)
        //
        //
        // Input: the argument "u" is an array of vectors just as before.
        //
        // Output: an array of vectors just like u.
        //         The expression
        //              ret[n](i)
        //         gives the calculated value of the nodal force at node
        //         n in the ith direction.
        //
        // 

        for (int q = 0; q < Q; q++) {
            Set::Matrix J = Set::Matrix::Zero();
            J = Jacobian(Qp[q]);
            Set::Matrix Jinv = Set::Matrix::Zero();
            Jinv = J.inverse();
            Set::Scalar Jdet = 0.0;
            Jdet = J.determinant();
            beta = Set::Matrix::Zero();
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < D; j++) {
                    for (int n = 0; n < N; n++) {
                        for (int k = 0; k < D; k++) {
                            beta(i,j) += u[n](i) * Deta(Qp[q])[n](k) * Jinv(k,j);
                        }
                    }
                }
            }
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < D; j++) {
                    for (int n = 0; n < N; n++) {
                        for (int k = 0; k < D; k++) {
                            ret[n](i) += Qw[q] * model[q].DW(beta)(i,j) * Deta(Qp[q])[n](k) * Jinv(k,j) * Jdet;
                        }
                    }
                }
            }
        }

        return ret;
    }

    std::array<std::array<Set::Matrix,N>,N> DDW(std::array<Set::Vector,N> &u)
    {
        std::array<std::array<Set::Matrix,N>,N> ret;
        for (int n = 0; n < N; n++) for (int m = 0; m < N; m++) ret[n][m] = Set::Matrix::Zero();
        Set::Matrix beta = Set::Matrix::Zero();
        const std::array<Set::Scalar,Q> &Qw = getQw();
        const std::array<Set::Vector,Q> &Qp = getQp();

        //
        // TODO: Calculate the second derivative of energy of the element as a function of
        //       nodal displacement.
        //       (See Problem 1 part b)
        //
        // Input: the argument "u" is an array of vectors just as before.
        //
        // Output: an array of arrays of tensors.
        //         The expression
        //              ret[m][n](i,j)
        //         returns the derivative of the mth force in the ith direction
        //         with respect to the displacement at the nth node in the jth direction.
        // 
        
        for (int q = 0; q < Q; q++) {
            Set::Matrix J = Set::Matrix::Zero();
            J = Jacobian(Qp[q]);
            Set::Matrix Jinv = Set::Matrix::Zero();
            Jinv = J.inverse();
            Set::Scalar Jdet = 0.0;
            Jdet = J.determinant();
            beta = Set::Matrix::Zero();
            for (int i = 0; i < D; i++) {
                for (int j = 0; j < D; j++) {
                    for (int n = 0; n < N; n++) {
                        for (int k = 0; k < D; k++) {
                            beta(i,j) += u[n](i) * Deta(Qp[q])[n](j) * Jinv(k,j);
                        }
                    }
                }
            }
            for (int m = 0; m < N; m++) {
                for (int n = 0; n < N; n++) {
                    for (int i = 0; i < D; i++) {
                        for (int j = 0; j < D; j++) {
                            for (int k = 0; k < D; k++) {
                                for (int l = 0; l < D; l++) {
                                    for (int r = 0; r < D; r++) {
                                        for (int s = 0; s < D; s++) {
                                            ret[m][n](i,j) += Qw[q] * model[q].DDW(beta)(i,l,j,k) * Deta(Qp[q])[n](r) * Deta(Qp[q])[m](s) * Jinv(r,l) * Jinv(s,k) * Jdet;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return ret;
    }



    virtual std::array<double,N> eta(Set::Vector y) = 0;

    virtual std::array<Set::Vector,N> Deta(Set::Vector y) = 0;

    Set::Matrix Jacobian(Set::Vector Y)
    {
        // Calculate Jacobian
        Eigen::Matrix2d jacobian = Eigen::Matrix2d::Zero();
        std::array<Eigen::Vector2d,N> deta = Deta(Y);

        for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++)
                for (int n = 0; n < N ; n++)
                {
                    jacobian(i,j) += (*X0)[id[n]](i) * deta[n](j);
                }
        return jacobian;
    }

public:
    constexpr int getN() {return N;}
    virtual const std::array<Set::Vector,N> & getY0() = 0;
    virtual const std::array<Set::Vector,Q> & getQp() = 0;
    virtual const std::array<double,Q> & getQw() = 0;
    virtual const double getArea() = 0;
    virtual const std::array<int,N> & getid() {return id;};

    static const int _D = D; // The DIMENSION of the space (2d or 3d)
    static const int _N = N; // The number of NODES in the element
    static const int _Q = Q; // The number of QUADRATURE POINTS in the element

private:
    std::vector<Set::Vector> *X0;
    std::array<int,N> id;
    std::array<MODEL,Q> model; 

public:

    static const int vtkid = 0;

};



}
#endif
