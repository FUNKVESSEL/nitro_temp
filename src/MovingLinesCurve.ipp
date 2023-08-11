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

#ifndef CODE_SRC_MOVINGLINESCURVE_IPP_
#define CODE_SRC_MOVINGLINESCURVE_IPP_

template< typename GEOM >
int nitro::movingLines::MovingLines< 1, GEOM >::evaluateRank( Eigen::VectorXd singularValues )
{
    int rank = 1;

    for (int i = 1; i < singularValues.size(); ++i) {
        if ( std::abs(singularValues(i) / singularValues(i - 1)) > tolRank_ ) {
            rank += 1;
        }
        else {
            break;
        }
    }

    return rank;
}

template< typename GEOM >
bool nitro::movingLines::MovingLines< 1, GEOM >::checkDegenerate( Eigen::MatrixXd matrixA, Eigen::MatrixXd matrixB )
{
    bool isDegenerate = false;

    int dim = matrixA.rows();

    int countDeficient = 0;

    for (int i = 0; i < dim + 1; ++i) {
        double randNumber = double( std::rand() ) / RAND_MAX;

        Eigen::MatrixXd randMatrix = matrixA - randNumber*matrixB;

        Eigen::JacobiSVD< Eigen::MatrixXd > matrixSVD(randMatrix, Eigen::ComputeFullV);

        Eigen::VectorXd singularValues = matrixSVD.singularValues();

        int randRank = this->evaluateRank(singularValues);

        if (randRank < dim) countDeficient += 1;

    }

    if (countDeficient == dim + 1) isDegenerate = true;

    return isDegenerate;
}

template< typename GEOM >
void nitro::movingLines::MovingLines< 1, GEOM >::computeCoefficients()
{
    int degree = curveMonomial_.degree();

    int spaceDim = curveMonomial_.spaceDim();

    if (isSpecified_ == false) {
        nu_ = degree - 1;

        if (nu_ == 0) nu_ = 1; // for the sake of inversion process
    }

    std::cout << "auxiliary degree: " << nu_ << std::endl;

    std::vector<Point> curvePoints = curveMonomial_.giveCoefficients();

    // compute the coefficient matrix from f.g = 0
    matrixCoeffs_ = Eigen::MatrixXd::Zero(nu_ + degree + 1, (spaceDim + 1)*(nu_ + 1));

    for (int j = 0; j < nu_ + 1; ++j) {
        for (int l = 0; l < degree + 1; ++l) {

            int rowId = j + l;

            Point point = curvePoints[l];
            Eigen::VectorXd coordW = point.coordW;

            for (int k = 0; k < spaceDim + 1; ++k) {

                int colId = j + k*(nu_ + 1);

                matrixCoeffs_(rowId, colId) += coordW(k);

            }

        }
    }

    // compute the null vectors of the coefficient matrix
    Eigen::JacobiSVD< Eigen::MatrixXd > matrixCoeffsSVD(matrixCoeffs_, Eigen::ComputeFullV);

    Eigen::MatrixXd matrixCoeffsV = matrixCoeffsSVD.matrixV();

    int rankMatrixCoeffs = matrixCoeffsSVD.rank();

    matrixCoeffsNull_ = matrixCoeffsV.block(0, rankMatrixCoeffs, matrixCoeffs_.cols(), matrixCoeffs_.cols() - rankMatrixCoeffs);

}

template< typename GEOM >
void nitro::movingLines::MovingLines< 1, GEOM >::computeImplicitMatrix( Eigen::VectorXd coord )
{
    int numRowsImplicit = this->numRowsImplicit();
    int numColsImplicit = this->numColsImplicit();

    int spaceDim = this->spaceDim();

    matrixImplicit_ = Eigen::MatrixXd::Zero(numRowsImplicit, numColsImplicit);

    for (int i = 0; i < spaceDim + 1; ++i) {
        Eigen::MatrixXd matrixImp = matrixCoeffsNull_.block(i*numRowsImplicit, 0, numRowsImplicit, numColsImplicit);

        if (i == spaceDim) {
            matrixImplicit_ += matrixImp;
        }
        else {
            matrixImplicit_ += coord(i)*matrixImp;
        }

    }

}

template< typename GEOM >
void nitro::movingLines::MovingLines< 1, GEOM >::inversionSVD( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol )
{
    Eigen::JacobiSVD< Eigen::MatrixXd > matrixSVD(matrix, Eigen::ComputeFullV);

    // evaluate nullity of the implicit matrix
    Eigen::VectorXd singularValues = matrixSVD.singularValues();

    int rankMatrix = this->evaluateRank(singularValues);

    if (rankMatrix < matrix.rows() && rankMatrix < matrix.cols()) {
        Eigen::MatrixXd matrixV = matrixSVD.matrixV();

        int numNullVecs = matrix.cols() - rankMatrix;

        Eigen::MatrixXd kernel = matrixV.block(0, matrix.cols() - numNullVecs, matrix.cols(), numNullVecs);

        Eigen::MatrixXd minorKernel1 = kernel.block(0, 0, numNullVecs, numNullVecs);
        Eigen::MatrixXd minorKernel2 = kernel.block(1, 0, numNullVecs, numNullVecs);

        Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolverInverse;
        eigenSolverInverse.compute(minorKernel2, minorKernel1);

        for (int i = 0; i < eigenSolverInverse.alphas().size(); ++i) {

            if ( std::abs(eigenSolverInverse.alphas()(i).imag()) < tol ) {

                double theta = eigenSolverInverse.alphas()(i).real() / eigenSolverInverse.betas()(i);

                Eigen::VectorXd param(1);

                param(0) = theta;

                params.push_back(param);

            }

        }

    }

}

template< typename GEOM >
void nitro::movingLines::MovingLines< 1, GEOM >::inversionEigenVec( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol )
{
    Eigen::MatrixXd kernel1 = matrix.block(0, 0, matrix.cols(), matrix.cols());
    Eigen::MatrixXd kernel2 = matrix.block(1, 0, matrix.cols(), matrix.cols());

    // check degenerate case (i.e. tangent in curve case using eigenvec method)
    bool isTangent = this->checkDegenerate(kernel2, kernel1);

    if (isTangent == true) std::cerr << "Warning: tangent case detected. SVD inversion is recommended." << std::endl;

    // compute generalised eigenvalue
    Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolverInverse;
    eigenSolverInverse.compute(kernel2, kernel1);

    for (int i = 0; i < eigenSolverInverse.alphas().size(); ++i) {

        if ( std::abs(eigenSolverInverse.alphas()(i).imag()) < tol ) {

            double theta = eigenSolverInverse.alphas()(i).real() / eigenSolverInverse.betas()(i);

            Eigen::VectorXd param(1);

            param(0) = theta;

            params.push_back(param);

        }

    }

}

#endif /* CODE_SRC_MOVINGLINESCURVE_IPP_ */
