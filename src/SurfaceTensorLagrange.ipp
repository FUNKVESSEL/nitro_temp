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

#ifndef CODE_SRC_SURFACETENSORLAGRANGE_IPP_
#define CODE_SRC_SURFACETENSORLAGRANGE_IPP_

template< typename P >
void nitro::Surface< nitro::TENSORLAGRANGE, P >::setCoefficients( std::vector<Point> coeffsL )
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
void nitro::Surface< nitro::TENSORLAGRANGE, P >::setIndexMap()
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
void nitro::Surface< nitro::TENSORLAGRANGE, P >::getMonomialcoeffs( std::vector<Point>& coeffsM )
{
    this->setIndexMap();

    int numCoeffs = coeffs_.size();

    // Vondermonde matrix (in terms of lex order starting from 1)
    Eigen::MatrixXd matrixVM = Eigen::MatrixXd::Zero(numCoeffs, (du_ + 1)*(dv_ + 1));

    for (int i = 0; i < numCoeffs; ++i) {
        std::pair<int, int> indexCoeff = indexCoeffs_[i];

        int index1 = indexCoeff.first, index2 = indexCoeff.second;

        int rowId = i;

        for (int ii = 0; ii < du_ + 1; ++ii) {

            double u = double(index1)/du_;

            if (ii == 0) u = 1.;

            for (int jj = 0; jj < dv_ + 1; ++jj) {

                double v = double(index2)/dv_;

                if (jj == 0) v = 1;

                int colId = jj + ii*(dv_ + 1);

                matrixVM(rowId, colId) = std::pow(u, ii)*std::pow(v, jj);

            }
        }
    }

    Eigen::MatrixXd matrixVMT = matrixVM.transpose();

    Eigen::MatrixXd matrixTrans = matrixVMT.inverse();

    // convert to monomial basis
    for (int i = 0; i < numCoeffs; ++i) {
        Point coeffM;

        coeffM.id = i;

        //-----------------------------------------------------------------------------------
        Eigen::VectorXd coordMW = Eigen::VectorXd::Zero(spaceDim_ + 1);

        for (int j = 0; j < numCoeffs; ++j) {
            Point coeffL = coeffs_[j];
            Eigen::VectorXd coordLW = coeffL.coordW;

            coordMW += matrixTrans(j, i)*coordLW;
        }

        coeffM.coordW = coordMW;

        //-----------------------------------------------------------------------------------
        Eigen::VectorXd coordM = Eigen::VectorXd::Zero(spaceDim_);

        for (int j = 0; j < spaceDim_; ++j) {
            coordM(j) = coordMW(j);
        }

        coeffM.coord = coordM;

        coeffsM.push_back(coeffM);

    }

}

template< typename P >
Eigen::VectorXd nitro::Surface< nitro::TENSORLAGRANGE, P >::evaluate( Eigen::VectorXd theta )
{
    this->setIndexMap();

    double theta1 = theta(0), theta2 = theta(1);

    Eigen::VectorXd coordTheta = Eigen::VectorXd::Zero(spaceDim_);

    int numCoeffs = coeffs_.size();

    for (int i = 0; i < numCoeffs; ++i) {
        Point coeff = coeffs_[i];

        Eigen::VectorXd coord = coeff.coord;

        std::pair<int, int> indexCoeff = indexCoeffs_[i];

        int index1 = indexCoeff.first, index2 = indexCoeff.second;

        double lagrangeValu = 1., lagrangeValv = 1.;
        for (int j = 0; j < du_ + 1; ++j) {
            if (j == index1) lagrangeValu *= 1.;
            else {
                lagrangeValu *= ( theta1 - double(j)/du_ ) / ( double(index1)/du_ - double(j)/du_ );
            }
        }
        for (int j = 0; j < dv_ + 1; ++j) {
            if (j == index2) lagrangeValv *= 1.;
            else {
                lagrangeValv *= ( theta2 - double(j)/dv_ ) / ( double(index2)/dv_ - double(j)/dv_ );
            }
        }

        coordTheta += lagrangeValu*lagrangeValv*coord;

    }

    return coordTheta;
}

#endif /* CODE_SRC_SURFACETENSORLAGRANGE_IPP_ */
