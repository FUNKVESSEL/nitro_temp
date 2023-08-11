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

#ifndef CODE_SRC_LINEINTERSECTOR_IPP_
#define CODE_SRC_LINEINTERSECTOR_IPP_

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::getConstMatrix()
{
    Eigen::VectorXd linePoint1 = line_.linePoints[0].coord;
    Eigen::VectorXd linePoint2 = line_.linePoints[1].coord;

    Eigen::MatrixXd matrixCoeffsNull = movingLines_.giveCoefficients();

    int numRowsMatrix = movingLines_.numRowsImplicit();
    int numColsMatrix = movingLines_.numColsImplicit();

    matrixA_ = Eigen::MatrixXd::Zero(numRowsMatrix, numColsMatrix);
    matrixB_ = Eigen::MatrixXd::Zero(numRowsMatrix, numColsMatrix);

    int spaceDim = movingLines_.spaceDim();

    for (int i = 0; i < spaceDim + 1; ++i) {
        Eigen::MatrixXd matrixImp = matrixCoeffsNull.block(i*numRowsMatrix, 0, numRowsMatrix, numColsMatrix);

        if (i == spaceDim) {
            matrixA_ += matrixImp;
        }
        else {
            matrixA_ += linePoint1(i)*matrixImp;
            matrixB_ += (linePoint1(i) - linePoint2(i))*matrixImp;
        }
    }

    std::cout << "implicit matrix dimension: (" << matrixA_.rows() << ", " << matrixA_.cols() << ")" << std::endl << std::endl;

}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::solverReduction()
{
    //===================================================================================================================
    // Pencil reduction process
    //===================================================================================================================
    Eigen::MatrixXd mA11 = matrixA_, mB11 = matrixB_;

    Eigen::MatrixXd mAReduced = mA11, mBReduced = mB11;

    if (mAReduced.rows() != mAReduced.cols()) {

        while ( true ) {

            int dm1 = mB11.rows(), dn1 = mB11.cols();

            // column compression of matrix B11
            Eigen::JacobiSVD< Eigen::MatrixXd > mB11SVD(mB11, Eigen::ComputeFullV);

            Eigen::VectorXd mB11Singular = mB11SVD.singularValues();

            int rankB11 = movingLines_.evaluateRank(mB11Singular);

            Eigen::MatrixXd mB11V = mB11SVD.matrixV();

            if (dm1 == dn1) { // the matrix becomes square already after row compression
                break;
            }
            else {
                if (rankB11 == dn1) { // the matrix has full column rank but not square, consider its transposed form
                    mA11.transposeInPlace();
                    mB11.transposeInPlace();

                    dm1 = mB11.rows(), dn1 = mB11.cols();

                    Eigen::JacobiSVD< Eigen::MatrixXd > mB11TSVD(mB11, Eigen::ComputeFullV);

                    Eigen::VectorXd mB11TSingular = mB11TSVD.singularValues();

                    rankB11 = movingLines_.evaluateRank(mB11TSingular);

                    mB11V = mB11TSVD.matrixV();

                }
            }

            // multiply mB11V to mA11 and mB11 and consider its block
            Eigen::MatrixXd mA11Tmp = mA11*mB11V;
            Eigen::MatrixXd mB11Tmp = mB11*mB11V;

            Eigen::MatrixXd mA21 = mA11Tmp.block(0, 0, dm1, rankB11);
            Eigen::MatrixXd mB21 = mB11Tmp.block(0, 0, dm1, rankB11);

            mA11 = mA11Tmp.block(0, rankB11, dm1, dn1 - rankB11);

            mAReduced = mA21, mBReduced = mB21;

            if (mBReduced.rows() == mBReduced.cols()) { // the matrix becomes square
                break;
            }
            else { // row number should be greater than column number; column number should be > 0
                Eigen::JacobiSVD< Eigen::MatrixXd > mA11SVD(mA11, Eigen::ComputeFullU);

                Eigen::VectorXd mA11Singular = mA11SVD.singularValues();

                int rankA11 = movingLines_.evaluateRank(mA11Singular);

                // matrix A11 can be row-compressed already
                bool rowCompressed = false;
                if (rankA11 == mA11.rows()) {
                    rowCompressed = true;
                }
                else { // consider special case of zero rows
                    int countNonZero = 0;
                    for (int i = rankA11; i < dm1; ++i) {
                        for (int j = 0; j < dn1 - rankB11; ++j) {
                            if ( std::abs(mA11(i, j)) > tol_ ) countNonZero += 1;
                        }
                    }

                    if (countNonZero == 0) rowCompressed = true;
                }

                // row compression of the new matrix A11
                if (rowCompressed == false) {
                    Eigen::MatrixXd mA11U = mA11SVD.matrixU();

                    Eigen::MatrixXd mA11UT = mA11U.transpose();

                    // multiply mA11UT to mA21 and mB21
                    Eigen::MatrixXd mA21Tmp = mA11UT*mA21;
                    Eigen::MatrixXd mB21Tmp = mA11UT*mB21;

                    Eigen::MatrixXd mA22 = mA21Tmp.block(rankA11, 0, dm1 - rankA11, rankB11);
                    Eigen::MatrixXd mB22 = mB21Tmp.block(rankA11, 0, dm1 - rankA11, rankB11);

                    mAReduced = mA22, mBReduced = mB22;

                }

                mA11 = mAReduced;
                mB11 = mBReduced;

            }

        }

    }

    // solve the generalised eigenvalue problem
    Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolver;

    eigenSolver.compute(mAReduced, mBReduced);

    std::vector<double> eigenvalues;
    for (int i = 0; i < eigenSolver.alphas().size(); ++i) {

        if ( std::abs(eigenSolver.alphas()(i).imag()) < tol_ ) {
            double eigenvalue = eigenSolver.alphas()(i).real() / eigenSolver.betas()(i);

            if ( std::abs(eigenSolver.betas()(i)) > tol_ ) { // exclude complexInfinity case
                eigenvalues.push_back(eigenvalue);
            }
        }
    }

    eigenvalues_ = eigenvalues;

    matrixA_.transposeInPlace();
    matrixB_.transposeInPlace();

}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::solverSquare()
{
    int numRowsMatrix = movingLines_.numRowsImplicit();

    Eigen::MatrixXd minorMatrixA = matrixA_.block(0, 0, numRowsMatrix, numRowsMatrix);
    Eigen::MatrixXd minorMatrixB = matrixB_.block(0, 0, numRowsMatrix, numRowsMatrix);

    // compute with transposed form
    minorMatrixA.transposeInPlace();
    minorMatrixB.transposeInPlace();

    // solve the generalised eigenvalue problem
    Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolver;

    typedef Eigen::GeneralizedEigenSolver< Eigen::MatrixXd >::EigenvectorsType EigenvectorsType;

    eigenSolver.compute(minorMatrixA, minorMatrixB);

    std::vector<double> eigenvalues;
    std::vector< Eigen::MatrixXd > eigenvectors;

    for (int i = 0; i < eigenSolver.alphas().size(); ++i) {

        if ( std::abs(eigenSolver.alphas()(i).imag()) < tol_ ) {
            double eigenvalue = eigenSolver.alphas()(i).real() / eigenSolver.betas()(i);

            if ( std::abs(eigenSolver.betas()(i)) > tol_ ) { // exclude complexInfinity case

                EigenvectorsType eigenvector = eigenSolver.eigenvectors().col(i);

                std::vector<double>::iterator iter = std::find_if(eigenvalues.begin(), eigenvalues.end(), find_cmp(eigenvalue, tol_));

                if (iter != eigenvalues.end() ) { // multiple pre-images
                    int pos = std::distance(eigenvalues.begin(), iter);

                    Eigen::MatrixXd eigenVec = eigenvectors[pos];

                    eigenVec.conservativeResize(eigenVec.rows(), eigenVec.cols() + 1);
                    eigenVec.col(eigenVec.cols() - 1) = eigenvector.real();

                    eigenvectors[pos] = eigenVec;
                }
                else {
                    eigenvalues.push_back(eigenvalue);

                    eigenvectors.push_back(eigenvector.real());
                }

            }
        }
    }

    eigenvalues_ = eigenvalues;

    eigenvectors_ = eigenvectors;

    matrixA_.resize(numRowsMatrix, numRowsMatrix);
    matrixB_.resize(numRowsMatrix, numRowsMatrix);

    matrixA_ = minorMatrixA;
    matrixB_ = minorMatrixB;

}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::solverSvd()
{
    int numRowsMatrix = movingLines_.numRowsImplicit();
	
	//concatenate vertically B and A; they have the same size:
	Eigen::MatrixXd matrixBA(2*numRowsMatrix, matrixA_.cols());
	matrixBA << matrixB_ , matrixA_ ;
	
	//Compute SVD of matrixBA, and the corresponding U matrix
	Eigen::JacobiSVD< Eigen::MatrixXd > matrixSVD(matrixBA, Eigen::ComputeFullU);
	Eigen::MatrixXd matrixU = matrixSVD.matrixU();
	
	//Extract the square pencil from U
    Eigen::MatrixXd minorMatrixA = matrixU.block(numRowsMatrix, 0, numRowsMatrix, numRowsMatrix);
    Eigen::MatrixXd minorMatrixB = matrixU.block(0, 0, numRowsMatrix, numRowsMatrix);

    // Compute with transposed form
    minorMatrixA.transposeInPlace();
    minorMatrixB.transposeInPlace();
	
    // solve the generalised eigenvalue problem
    Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolver;

    typedef Eigen::GeneralizedEigenSolver< Eigen::MatrixXd >::EigenvectorsType EigenvectorsType;

    eigenSolver.compute(minorMatrixA, minorMatrixB);

    std::vector<double> eigenvalues;
    std::vector< Eigen::MatrixXd > eigenvectors;

    for (int i = 0; i < eigenSolver.alphas().size(); ++i) {

        if ( std::abs(eigenSolver.alphas()(i).imag()) < tol_ ) {
            double eigenvalue = eigenSolver.alphas()(i).real() / eigenSolver.betas()(i);

            if ( std::abs(eigenSolver.betas()(i)) > tol_ ) { // exclude complexInfinity case

                EigenvectorsType eigenvector = eigenSolver.eigenvectors().col(i);

                std::vector<double>::iterator iter = std::find_if(eigenvalues.begin(), eigenvalues.end(), find_cmp(eigenvalue, tol_));

                if (iter != eigenvalues.end() ) { // multiple pre-images
                    int pos = std::distance(eigenvalues.begin(), iter);

                    Eigen::MatrixXd eigenVec = eigenvectors[pos];

                    eigenVec.conservativeResize(eigenVec.rows(), eigenVec.cols() + 1);
                    eigenVec.col(eigenVec.cols() - 1) = eigenvector.real();

                    eigenvectors[pos] = eigenVec;
                }
                else {
                    eigenvalues.push_back(eigenvalue);

                    eigenvectors.push_back(eigenvector.real());
                }
                
            }
        }
    }

    eigenvalues_ = eigenvalues;

    eigenvectors_ = eigenvectors;

    matrixA_.resize(numRowsMatrix, numRowsMatrix);
    matrixB_.resize(numRowsMatrix, numRowsMatrix);

    matrixA_ = minorMatrixA;
    matrixB_ = minorMatrixB;
    
}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::solverQr()
{
    int numRowsMatrix = movingLines_.numRowsImplicit();
		
	//Compute a QR decomposition of the transpose of B
	Eigen::HouseholderQR< Eigen::MatrixXd > qr(matrixB_.transpose());
	qr.compute(matrixB_.transpose());
	Eigen::MatrixXd q = qr.householderQ();
	Eigen::MatrixXd r = qr.matrixQR().triangularView<Eigen::Upper>();
	
	//Extract the square pencil
    Eigen::MatrixXd minorMatrixA = (matrixA_*q).block(0, 0, numRowsMatrix, numRowsMatrix);
    Eigen::MatrixXd minorMatrixB = r.block(0, 0, numRowsMatrix, numRowsMatrix).transpose();

    // compute with transposed form
    minorMatrixA.transposeInPlace();
    minorMatrixB.transposeInPlace();
	
    // solve the generalised eigenvalue problem
    Eigen::GeneralizedEigenSolver< Eigen::MatrixXd > eigenSolver;

    typedef Eigen::GeneralizedEigenSolver< Eigen::MatrixXd >::EigenvectorsType EigenvectorsType;

    eigenSolver.compute(minorMatrixA, minorMatrixB);

    std::vector<double> eigenvalues;
    std::vector< Eigen::MatrixXd > eigenvectors;

    for (int i = 0; i < eigenSolver.alphas().size(); ++i) {

        if ( std::abs(eigenSolver.alphas()(i).imag()) < tol_ ) {
            double eigenvalue = eigenSolver.alphas()(i).real() / eigenSolver.betas()(i);

            if ( std::abs(eigenSolver.betas()(i)) > tol_ ) { // exclude complexInfinity case

                EigenvectorsType eigenvector = eigenSolver.eigenvectors().col(i);

                std::vector<double>::iterator iter = std::find_if(eigenvalues.begin(), eigenvalues.end(), find_cmp(eigenvalue, tol_));

                if (iter != eigenvalues.end() ) { // multiple pre-images
                    int pos = std::distance(eigenvalues.begin(), iter);

                    Eigen::MatrixXd eigenVec = eigenvectors[pos];

                    eigenVec.conservativeResize(eigenVec.rows(), eigenVec.cols() + 1);
                    eigenVec.col(eigenVec.cols() - 1) = eigenvector.real();

                    eigenvectors[pos] = eigenVec;
                }
                else {
                    eigenvalues.push_back(eigenvalue);

                    eigenvectors.push_back(eigenvector.real());
                }
                
            }
        }
    }

    eigenvalues_ = eigenvalues;

    eigenvectors_ = eigenvectors;

    matrixA_.resize(numRowsMatrix, numRowsMatrix);
    matrixB_.resize(numRowsMatrix, numRowsMatrix);

    matrixA_ = minorMatrixA;
    matrixB_ = minorMatrixB;
    
}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::computeIntersectionSVD()
{
    Eigen::VectorXd linePoint1 = line_.linePoints[0].coord;
    Eigen::VectorXd linePoint2 = line_.linePoints[1].coord;

    int numIntersectPoints = 0;

    // inversion process
    for (int i = 0; i < eigenvalues_.size(); ++i) {
        double lineParam = eigenvalues_[i];

        Eigen::VectorXd coordIntersect = linePoint1 + lineParam*(linePoint2 - linePoint1);

        // intersection point candidate
        IntersectPoint intersectPoint;

        intersectPoint.id = numIntersectPoints;
        intersectPoint.coord = coordIntersect;
        intersectPoint.paramLine = lineParam;

        bool isIntersect = false;

        // implicit matrix at the point (or truncated square one if "SQUARE" method is used)
        Eigen::MatrixXd matrixIntersect = matrixA_ - lineParam*matrixB_; // M = A - xi*B

        std::vector< Eigen::VectorXd > paramsIntersect;

        movingLines_.inversionSVD(matrixIntersect, paramsIntersect, tol_);

        for (int ii = 0; ii < paramsIntersect.size(); ++ii) {
            Eigen::VectorXd param = paramsIntersect[ii];

            Eigen::VectorXd coordParam = movingLines_.evaluateGeometry(param);

            Eigen::VectorXd coordDiff = coordParam - coordIntersect;

            double dist = coordDiff.norm();

            if ( std::abs(dist) < tol_) {
                intersectPoint.paramGeom.push_back(param);
                isIntersect = true;
            }

        }

        if (isIntersect == true) intersectPoints_.push_back(intersectPoint);

    }

}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::computeIntersectionEigenVec()
{
    if (method_ == "REDUCTION") {
        std::cerr << "ERROR: EIGENVEC inversion NOT valid for REDUCTION method! Select SVD inversion." << std::endl;
    }

    Eigen::VectorXd linePoint1 = line_.linePoints[0].coord;
    Eigen::VectorXd linePoint2 = line_.linePoints[1].coord;

    int numIntersectPoints = 0;

    // inversion process
    for (int i = 0; i < eigenvalues_.size(); ++i) {
        double lineParam = eigenvalues_[i];

        Eigen::VectorXd coordIntersect = linePoint1 + lineParam*(linePoint2 - linePoint1);

        // intersection point candidate
        IntersectPoint intersectPoint;

        intersectPoint.id = numIntersectPoints;
        intersectPoint.coord = coordIntersect;
        intersectPoint.paramLine = lineParam;

        bool isIntersect = false;

        std::vector< Eigen::VectorXd > paramsIntersect;

        movingLines_.inversionEigenVec(eigenvectors_[i], paramsIntersect, tol_);

        for (int ii = 0; ii < paramsIntersect.size(); ++ii) {
            Eigen::VectorXd param = paramsIntersect[ii];

            Eigen::VectorXd coordParam = movingLines_.evaluateGeometry(param);

            Eigen::VectorXd coordDiff = coordParam - coordIntersect;

            double dist = coordDiff.norm();

            if ( std::abs(dist) < tol_) {
                intersectPoint.paramGeom.push_back(param);
                isIntersect = true;
            }

        }

        if (isIntersect == true) intersectPoints_.push_back(intersectPoint);

    }

}

template< typename MLINE, typename LINE >
void nitro::lineIntersector::LineIntersector<MLINE, LINE>::giveIntersection( std::vector<IntersectPoint>& intersectPoints )
{
    int numIntersects = intersectPoints_.size();

    int numIntersectRemain = 0;

    for (int i = 0; i < intersectPoints_.size(); ++i) {
        IntersectPoint ip = intersectPoints_[i];
        Eigen::VectorXd coord = ip.coord;

        if (i == 0) {
            ip.id = numIntersectRemain;
            intersectPoints.push_back(ip);
            numIntersectRemain += 1;
        }
        else {

            int countDuplicate = 0;

            for (int j = 0; j < intersectPoints.size(); ++j) {
                IntersectPoint ipRemain = intersectPoints[j];
                Eigen::VectorXd coordRemain = ipRemain.coord;

                Eigen::VectorXd coordDiff = coordRemain - coord;

                double dist = coordDiff.norm();

                if ( std::abs(dist) < tol_ ) {
                    countDuplicate += 1;
                }

            }

            if (countDuplicate == 0) {
                ip.id = numIntersectRemain;
                intersectPoints.push_back(ip);
                numIntersectRemain += 1;
            }

        }

    }
}

#endif /* CODE_SRC_LINEINTERSECTOR_IPP_ */
