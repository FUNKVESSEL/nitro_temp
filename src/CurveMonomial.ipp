//! Copyright (c) 2019 University of Cambridge (UK), INRIA (France)
//! All rights reserved.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
//! This file is part of nitro.
//! @Author: Xiao Xiao, Laurent Busé and Fehmi Cirak
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef CODE_SRC_CURVEMONOMIAL_IPP_
#define CODE_SRC_CURVEMONOMIAL_IPP_

template< typename P >
Eigen::VectorXd nitro::Curve< nitro::MONOMIAL, P >::evaluate(double theta)
{
    Eigen::VectorXd coordTheta = Eigen::VectorXd::Zero(spaceDim_);

    int numCoeffs = coeffs_.size();

    for (int i = 0; i < numCoeffs; ++i) {
        Point coeff = coeffs_[i];

        Eigen::VectorXd coord = coeff.coord;

        double basisVal = 1.;

        if (i > 0) basisVal = std::pow(theta, i);

        coordTheta += basisVal*coord;
    }

    return coordTheta;

}


#endif /* CODE_SRC_CURVEMONOMIAL_IPP_ */
