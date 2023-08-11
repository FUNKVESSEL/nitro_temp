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

#ifndef CODE_SRC_SURFACE_HPP_
#define CODE_SRC_SURFACE_HPP_

namespace nitro {

    enum surface {
        TENSORMONOMIAL,
        TENSORLAGRANGE,
        TENSORBERNSTEIN,
        TRIBERNSTEIN
    };

    template< nitro::surface SURF, typename P >
    class Surface;

    template< typename P >
    class Surface< nitro::TENSORMONOMIAL, P >;

    template< typename P >
    class Surface< nitro::TENSORLAGRANGE, P >;

}

template< typename P >
class nitro::Surface< nitro::TENSORMONOMIAL, P >
{
public:
    typedef P   Point;

    static const int dim = 2;

public:
    Surface( ) { }

    Surface(int du, int dv) {
        du_ = du;
        dv_ = dv;
        spaceDim_ = 3;
    }

    void setCoefficients( std::vector<Point> coeffsM ) { coeffs_ = coeffsM; }

    std::vector<Point> giveCoefficients() { return coeffs_; }

    void setIndexMap();

    void setMonomialIndexMap();

    std::vector< std::pair<int, int> > giveIndexMap() { this->setIndexMap(); return indexCoeffs_; }

    std::vector< std::pair<int, int> > giveMonomialIndexMap() { this->setMonomialIndexMap(); return indexMonomials_; }

    Eigen::VectorXd evaluate( Eigen::VectorXd theta );

public:
    std::pair<int, int> degree() {
        std::pair<int, int> deg = std::make_pair(du_, dv_);
        return deg;
    }

private:
    int du_, dv_; // bi-degree of the surface

    int spaceDim_;

    std::vector<Point> coeffs_;

    std::vector< std::pair<int, int> > indexCoeffs_; // index of monomial coefficients (points)

    std::vector< std::pair<int, int> > indexMonomials_; // index of monomial basis
};

template< typename P >
class nitro::Surface< nitro::TENSORLAGRANGE, P >
{
public:
    typedef P   Point;

    static const int dim = 2;

public:
    Surface( ) { }

    Surface(int du, int dv) {
        du_ = du;
        dv_ = dv;
        spaceDim_ = 3;
    }

    void setCoefficients( std::vector<Point> coeffsL ); // including the weights

    void setIndexMap();

    void getMonomialcoeffs( std::vector<Point>& coeffsM );

    std::vector< std::pair<int, int> > giveIndexMap() { return indexCoeffs_; }

    Eigen::VectorXd evaluate( Eigen::VectorXd theta );

public:
    std::pair<int, int> degree() {
        std::pair<int, int> deg = std::make_pair(du_, dv_);
        return deg;
    }

private:
    int du_, dv_; // bi-degree of the surface

    int spaceDim_;

    std::vector<Point> coeffs_;

    std::vector< std::pair<int, int> > indexCoeffs_; // index of the coefficient (point)
};

#include "SurfaceTensorMonomial.ipp"
#include "SurfaceTensorLagrange.ipp"

#endif /* CODE_SRC_SURFACE_HPP_ */
