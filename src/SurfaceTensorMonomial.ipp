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

#ifndef CODE_SRC_SURFACETENSORMONOMIAL_IPP_
#define CODE_SRC_SURFACETENSORMONOMIAL_IPP_

template< typename P >
void nitro::Surface< nitro::TENSORMONOMIAL, P >::setIndexMap()
{
    indexCoeffs_.clear();
    for (int j = 0; j < dv_ + 1; ++j) {
        for (int i = 0; i < du_ + 1; ++i) {

            std::pair<int, int> indexCoeff = std::make_pair(i, j);
            indexCoeffs_.push_back(indexCoeff);

        }
    }
}

template< typename P >
void nitro::Surface< nitro::TENSORMONOMIAL, P >::setMonomialIndexMap()
{
    indexMonomials_.clear();
    for (int i = 0; i < du_ + 1; ++i) {
        for (int j = 0; j < dv_ + 1; ++j) {

            std::pair<int, int> indexCoeff = std::make_pair(i, j);
            indexMonomials_.push_back(indexCoeff);

        }
    }
}

template< typename P >
Eigen::VectorXd nitro::Surface< nitro::TENSORMONOMIAL, P >::evaluate( Eigen::VectorXd theta )
{
    this->setMonomialIndexMap();

    double theta1 = theta(0), theta2 = theta(1);

    Eigen::VectorXd coordTheta = Eigen::VectorXd::Zero(spaceDim_);

    int numCoeffs = coeffs_.size();

    for (int i = 0; i < numCoeffs; ++i) {

        Point coeff = coeffs_[i];

        Eigen::VectorXd coord = coeff.coord;

        std::pair<int, int> indexCoeff = indexMonomials_[i];
        int index1 = indexCoeff.first, index2 = indexCoeff.second;

        double basisValu = 1., basisValv = 1.;

        if (index1 > 0) basisValu = std::pow(theta1, index1);

        if (index2 > 0) basisValv = std::pow(theta2, index2);

        coordTheta += basisValu*basisValv*coord;

    }

    return coordTheta;

}

#endif /* CODE_SRC_SURFACETENSORMONOMIAL_IPP_ */
