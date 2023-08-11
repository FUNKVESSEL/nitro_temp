//! Copyright (c) 2019 University of Cambridge (UK), INRIA (France)
//! All rights reserved.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
//! This file is part of nitro.
//! @Authors: Xiao Xiao, Laurent Busé and Fehmi Cirak
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef CODE_SRC_CURVE_HPP_
#define CODE_SRC_CURVE_HPP_

namespace nitro {

    enum curve {
        MONOMIAL,
        LAGRANGE,
        BERNSTEIN
    };

    template< nitro::curve CURVE, typename P >
    class Curve;

    template< typename P >
    class Curve< nitro::MONOMIAL, P >;

    template< typename P >
    class Curve< nitro::LAGRANGE, P >;

}

template< typename P >
class nitro::Curve< nitro::MONOMIAL, P >
{
public:
    typedef P   Point;

public:
    Curve( ) { }

    Curve(int spaceDim) {
        spaceDim_ = spaceDim;
    }

    void setCoefficients(std::vector<Point> coeffsM) { coeffs_ = coeffsM; }

    std::vector<Point> giveCoefficients() { return coeffs_; }

    Eigen::VectorXd evaluate( double theta );

    int degree() { return coeffs_.size() - 1; }

    int spaceDim() { return spaceDim_; }

private:
    int spaceDim_;

    std::vector<Point> coeffs_;
};

template< typename P >
class nitro::Curve< nitro::LAGRANGE, P >
{
public:
    typedef P   Point;

    static const int dim = 1;

public:
    Curve( ) { }

    Curve(int spaceDim) {
        spaceDim_ = spaceDim;
    }

    void setCoefficients(std::vector<Point> coeffsL); // including the weights

    int degree() { return coeffs_.size() - 1; }

    int spaceDim() { return spaceDim_; }

    void getMonomialCoeffs(std::vector<Point>& coeffsM);

    Eigen::VectorXd evaluate( double theta );

private:
    int spaceDim_;

    std::vector<Point> coeffs_;

};

#include "CurveMonomial.ipp"
#include "CurveLagrange.ipp"

#endif /* CODE_SRC_CURVE_HPP_ */
