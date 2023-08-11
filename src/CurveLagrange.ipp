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

#ifndef CODE_SRC_CURVELAGRANGE_IPP_
#define CODE_SRC_CURVELAGRANGE_IPP_

template< typename P >
void nitro::Curve< nitro::LAGRANGE, P >::setCoefficients(std::vector<Point> coeffsL)
{
    for (int i = 0; i < coeffsL.size(); ++i) {
        Point coeff = coeffsL[i];

        Eigen::VectorXd coord = coeff.coord;

        Eigen::VectorXd coordW = Eigen::VectorXd::Zero(spaceDim_ + 1);

        for (int j = 0; j < spaceDim_ + 1; ++j) {
            if (j == spaceDim_) {
                coordW(j) = 1.;
            }
            else {
                coordW(j) = coord(j);
            }
        }

        coeff.coordW = coordW;

        coeffs_.push_back(coeff);

    }
}

template< typename P >
void nitro::Curve< nitro::LAGRANGE, P >::getMonomialCoeffs(std::vector<Point>& coeffsM)
{
    int numCoeffs = coeffs_.size();
    int degree = this->degree();

    // Vondermonde matrix
    Eigen::MatrixXd matrixVM = Eigen::MatrixXd::Zero(numCoeffs, degree + 1);

    for (int i = 0; i < numCoeffs; ++i) {
        for (int j = 0; j < degree + 1; ++j) {

            double param = double(i)/degree;

            if (j == 0) {
                matrixVM(i, j) = 1.;
            }
            else {
                matrixVM(i, j) = std::pow(param, j);
            }

        }
    }

    Eigen::MatrixXd matrixVMT = matrixVM.transpose();

    Eigen::MatrixXd matrixTrans = matrixVMT.inverse();

    for (int i = 0; i < numCoeffs; ++i) {
        Point coeffM;
        coeffM.id = i;

        //-----------------------------------------------------------------------
        Eigen::VectorXd coordMW = Eigen::VectorXd::Zero(spaceDim_ + 1);

        for (int j = 0; j < numCoeffs; ++j) {
            Point coeffL = coeffs_[j];
            Eigen::VectorXd coordLW = coeffL.coordW;

            coordMW += matrixTrans(j, i)*coordLW;
        }

        coeffM.coordW = coordMW;

        //-----------------------------------------------------------------------
        Eigen::VectorXd coordM = Eigen::VectorXd::Zero(spaceDim_);

        for (int j = 0; j < spaceDim_; ++j) {
            coordM(j) = coordMW(j);
        }

        coeffM.coord = coordM;

        coeffsM.push_back(coeffM);

    }

}

template< typename P >
Eigen::VectorXd nitro::Curve< nitro::LAGRANGE, P >::evaluate(double theta)
{
    Eigen::VectorXd coordTheta = Eigen::VectorXd::Zero(spaceDim_);

    int numCoeffs = coeffs_.size();
    int degree = this->degree();

    for (int i = 0; i < numCoeffs; ++i) {
        Point coeff = coeffs_[i];

        Eigen::VectorXd coord = coeff.coord;

        double lagrangeVal = 1.;
        for (int j = 0; j < numCoeffs; ++j) {

            if (j == i) lagrangeVal *= 1.;
            else {
                lagrangeVal *= ( theta - double(j)/degree ) / ( double(i)/degree - double(j)/degree );
            }
        }

        coordTheta += lagrangeVal*coord;

    }

    return coordTheta;
}


#endif /* CODE_SRC_CURVELAGRANGE_IPP_ */
