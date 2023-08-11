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

#ifndef CODE_SRC_MOVINGLINESSURFACE_IPP_
#define CODE_SRC_MOVINGLINESSURFACE_IPP_

template< typename GEOM >
void nitro::movingLines::MovingLines< 2, GEOM >::setIndexMapAux()
{
    this->giveSurfaceDegree();

    int du = degree_.first;
    int dv = degree_.second;

    if (isSpecified_ == false) {
        if (du <= dv) {
            nu1_ = 2*du - 1;
            nu2_ = dv - 1;

            if (nu2_ == 0) nu2_ = 1; // for the sake of inversion process
        }
        else {
            nu1_ = du - 1;
            nu2_ = 2*dv - 1;

            if (nu1_ == 0) nu1_ = 1; // for the sake of inversion process
        }
    }

    std::cout << "auxiliary degree: (" << nu1_ << ", " << nu2_ << ")" << std::endl;

    for (int i = 0; i < nu1_ + 1; ++i) {
        for (int j = 0; j < nu2_ + 1; ++j) {

            std::pair<int, int> indexAux = std::make_pair(i, j);
            indexMapAux_.push_back(indexAux);

        }
    }

}

template< typename GEOM >
int nitro::movingLines::MovingLines< 2, GEOM >::evaluateRank( Eigen::VectorXd singularValues )
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
bool nitro::movingLines::MovingLines< 2, GEOM >::checkDegenerate( Eigen::MatrixXd matrixA, Eigen::MatrixXd matrixB )
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
void nitro::movingLines::MovingLines< 2, GEOM >::computeCoefficients()
{
    this->setIndexMapAux();

    int du = surfaceMonomial_.degree().first;
    int dv = surfaceMonomial_.degree().second;

    std::vector<Point> surfPoints = surfaceMonomial_.giveCoefficients();

    std::vector< std::pair<int, int> > indexPoints = surfaceMonomial_.giveMonomialIndexMap();

    // compute the coefficient matrix from f.g = 0
    matrixCoeffs_ = Eigen::MatrixXd::Zero((du + nu1_ + 1)*(dv + nu2_ + 1), (spaceDim_ + 1)*(nu1_ + 1)*(nu2_ + 1));

    for (int j = 0; j < (nu1_ + 1)*(nu2_ + 1); ++j) {

        int j1 = indexMapAux_[j].first, j2 = indexMapAux_[j].second;

        for (int l = 0; l < (du + 1)*(dv + 1); ++l) {

            int l1 = indexPoints[l].first, l2 = indexPoints[l].second;

            int rowId = (j2 + l2) + (j1 + l1)*(dv + nu2_ + 1);

            Point point = surfPoints[l];
            Eigen::VectorXd coordW = point.coordW;

            for (int k = 0; k < spaceDim_ + 1; ++k) {

                int colId = j + k*(nu1_ + 1)*(nu2_ + 1);

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
void nitro::movingLines::MovingLines< 2, GEOM >::inversionSVD( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol )
{
    Eigen::JacobiSVD< Eigen::MatrixXd > matrixSVD(matrix, Eigen::ComputeFullV);

    // evaluate nullity of the implicit matrix
    Eigen::VectorXd singularValues = matrixSVD.singularValues();

    int rankMatrix = this->evaluateRank(singularValues);

    if (rankMatrix < matrix.rows() && rankMatrix < matrix.cols()) {
        Eigen::MatrixXd matrixV = matrixSVD.matrixV();

        int numNullVecs = matrix.cols() - rankMatrix;

        Eigen::MatrixXd kernel = matrixV.block(0, matrix.cols() - numNullVecs, matrix.cols(), numNullVecs);

        Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolverInverse;

        //--------------------------------------------------------------------------------------------
        // solve v-param
        //--------------------------------------------------------------------------------------------
        Eigen::MatrixXd minorKernelv1 = Eigen::MatrixXd::Zero(numNullVecs, numNullVecs);
        Eigen::MatrixXd minorKernelv2 = Eigen::MatrixXd::Zero(numNullVecs, numNullVecs);

        // (highest order of v is carefully chosen in case of multiple pre-images)
        if (numNullVecs > nu2_*(nu1_ + 1)) std::cerr << "Warning: need to increase the auxiliary degree in v." << std::endl;

        int rowId  = 0;
        for (int i = 0; i < kernel.rows(); ++i) {
            if ( (i + 1) % (nu2_ + 1) != 0 ) {
                for (int j = 0; j < numNullVecs; ++j) {
                    minorKernelv1(rowId, j) = kernel(i, j);
                    minorKernelv2(rowId, j) = kernel(i + 1, j);
                }
                rowId += 1;
            }

            if (rowId == numNullVecs) break;
        }

        // compute eigenvalues v
        eigenSolverInverse.compute(minorKernelv2, minorKernelv1);

        // check degenerate case
        bool isDegeneratev = this->checkDegenerate(minorKernelv2, minorKernelv1);

        std::vector<double> thetas2;
        for (int i = 0; i < eigenSolverInverse.alphas().size(); ++i) {

            if ( std::abs( eigenSolverInverse.alphas()(i).imag() ) < tol ) {

                double theta2 = eigenSolverInverse.alphas()(i).real() / eigenSolverInverse.betas()(i);

                // effective eigenvalues (exclude the infinity value)
                if ( std::abs( eigenSolverInverse.betas()(i) ) > tol ) {
                    thetas2.push_back(theta2);
                }
                else {
                    if ( std::abs( eigenSolverInverse.alphas()(i).real() ) < tol ) {
                        thetas2.push_back(theta2);
                    }
                }

            }
        }

        //--------------------------------------------------------------------------------------------
        // solve u-param
        //--------------------------------------------------------------------------------------------
        if (numNullVecs > (nu2_ + 1)*nu1_) std::cerr << "Warning: need to increase the auxiliary degree in u." << std::endl;

        Eigen::MatrixXd minorKernelu1 = kernel.block(0, 0, numNullVecs, numNullVecs);
        Eigen::MatrixXd minorKernelu2 = kernel.block(nu2_ + 1, 0, numNullVecs, numNullVecs);

        // compute eigenvalues u
        eigenSolverInverse.compute(minorKernelu2, minorKernelu1);

        // check degenerate case
        bool isDegenerateu = this->checkDegenerate(minorKernelu2, minorKernelu1);

        std::vector<double> thetas1;
        for (int i = 0; i < eigenSolverInverse.alphas().size(); ++i) {

            if ( std::abs( eigenSolverInverse.alphas()(i).imag() ) < tol ) {

                double theta1 = eigenSolverInverse.alphas()(i).real() / eigenSolverInverse.betas()(i);

                // effective eigenvalues (exclude the infinity value)
                if ( std::abs( eigenSolverInverse.betas()(i) ) > tol ) {
                    thetas1.push_back(theta1);
                }
                else {
                    if ( std::abs( eigenSolverInverse.alphas()(i).real() ) < tol  ) {
                        thetas1.push_back(theta1);
                    }
                }

            }
        }

        //--------------------------------------------------------------------------------------------
        if (isDegeneratev == true or isDegenerateu == true) std::cerr << "Warning: degenerate case detected. Degenerate case not yet supported." << std::endl;

        //if (thetas1.size() != thetas2.size()) std::cerr << "Warning in inversion process! (Dimensions of two parameters do not match!)" << std::endl << std::endl;

        // remove duplicate u and v before combination (not a must)
        double tolerance = tolRank_; // use some small tolerance

        std::sort(thetas1.begin(), thetas1.end());
        thetas1.erase( std::unique( thetas1.begin(), thetas1.end(), unique_cmp(tolerance) ), thetas1.end() );

        std::sort(thetas2.begin(), thetas2.end());
        thetas2.erase( std::unique( thetas2.begin(), thetas2.end(), unique_cmp(tolerance) ), thetas2.end() );

        // return all possible (u, v)
        for (int i = 0; i < thetas1.size(); ++i) {
            for (int j = 0; j < thetas2.size(); ++j) {
                double theta1 = thetas1[i];
                double theta2 = thetas2[j];

                Eigen::VectorXd param(2);
                param(0) = theta1, param(1) = theta2;

                params.push_back(param);
            }

        }

    }

}

template< typename GEOM >
void nitro::movingLines::MovingLines< 2, GEOM >::inversionEigenVec( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol )
{
    int numNullVecs = matrix.cols();

    Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolverInverse;

    //------------------------------------------------------------------------------------------------
    // solve v-param
    //------------------------------------------------------------------------------------------------
    Eigen::MatrixXd kernelv1 = Eigen::MatrixXd::Zero(numNullVecs, numNullVecs);
    Eigen::MatrixXd kernelv2 = Eigen::MatrixXd::Zero(numNullVecs, numNullVecs);

    // (highest order of v is carefully chosen in case of multiple pre-images)
    if (numNullVecs > nu2_*(nu1_ + 1)) std::cerr << "Warning: need to increase the auxiliary degree in v." << std::endl;

    int rowId  = 0;
    for (int i = 0; i < matrix.rows(); ++i) {
        if ( (i + 1) % (nu2_ + 1) != 0 ) {
            for (int j = 0; j < numNullVecs; ++j) {
                kernelv1(rowId, j) = matrix(i, j);
                kernelv2(rowId, j) = matrix(i + 1, j);
            }
            rowId += 1;
        }

        if (rowId == numNullVecs) break;
    }

    // compute eigenvalues v
    eigenSolverInverse.compute(kernelv2, kernelv1);

    // check degenerate case
    bool isDegeneratev = this->checkDegenerate(kernelv2, kernelv1);

    std::vector<double> thetas2;
    for (int i = 0; i < eigenSolverInverse.alphas().size(); ++i) {

        if ( std::abs( eigenSolverInverse.alphas()(i).imag() ) < tol ) {

            double theta2 = eigenSolverInverse.alphas()(i).real() / eigenSolverInverse.betas()(i);

            // effective eigenvalues (exclude the infinity value)
            if ( std::abs( eigenSolverInverse.betas()(i) ) > tol ) {
                thetas2.push_back(theta2);
            }
            else {
                if ( std::abs( eigenSolverInverse.alphas()(i).real() ) < tol ) {
                    thetas2.push_back(theta2);
                }
            }

        }
    }

    //------------------------------------------------------------------------------------------------
    // solve u-param
    //------------------------------------------------------------------------------------------------
    if (numNullVecs > (nu2_ + 1)*nu1_) std::cerr << "Warning: need to increase the auxiliary degree in u." << std::endl;

    Eigen::MatrixXd kernelu1 = matrix.block(0, 0, matrix.cols(), matrix.cols());
    Eigen::MatrixXd kernelu2 = matrix.block(nu2_ + 1, 0, matrix.cols(), matrix.cols());

    // compute eigenvalues u
    eigenSolverInverse.compute(kernelu2, kernelu1);

    // check degenerate case
    bool isDegenerateu = this->checkDegenerate(kernelu2, kernelu1);

    std::vector<double> thetas1;
    for (int i = 0; i < eigenSolverInverse.alphas().size(); ++i) {

        if ( std::abs( eigenSolverInverse.alphas()(i).imag() ) < tol ) {

            double theta1 = eigenSolverInverse.alphas()(i).real() / eigenSolverInverse.betas()(i);

            // effective eigenvalues (exclude the infinity value)
            if ( std::abs( eigenSolverInverse.betas()(i) ) > tol ) {
                thetas1.push_back(theta1);
            }
            else {
                if ( std::abs( eigenSolverInverse.alphas()(i).real() ) < tol ) {
                    thetas1.push_back(theta1);
                }
            }

        }
    }

    //------------------------------------------------------------------------------------------------
    if (isDegeneratev == true or isDegenerateu == true) {
        std::cerr << "Warning: degenerate case or tangent case detected." << std::endl;
        std::cerr << "- SVD inversion is recommended in tangent case;" << std::endl;
        std::cerr << "- Degenerate case not yet supported." << std::endl;
    }

    //if (thetas1.size() != thetas2.size()) std::cerr << "Warning in inversion process! (Dimensions of two parameters do not match!)" << std::endl << std::endl;

    // remove duplicate u and v before combination (not a must)
    double tolerance = tolRank_; // use some small tolerance

    std::sort(thetas1.begin(), thetas1.end());
    thetas1.erase( std::unique( thetas1.begin(), thetas1.end(), unique_cmp(tolerance) ), thetas1.end() );

    std::sort(thetas2.begin(), thetas2.end());
    thetas2.erase( std::unique( thetas2.begin(), thetas2.end(), unique_cmp(tolerance) ), thetas2.end() );

    // return all possible (u, v)
    for (int i = 0; i < thetas1.size(); ++i) {
        for (int j = 0; j < thetas2.size(); ++j) {
            double theta1 = thetas1[i];
            double theta2 = thetas2[j];

            Eigen::VectorXd param(2);
            param(0) = theta1, param(1) = theta2;

            params.push_back(param);
        }

    }

}

#endif /* CODE_SRC_MOVINGLINESSURFACE_IPP_ */
