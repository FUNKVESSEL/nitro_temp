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

#ifndef CODE_SRC_LINEINTERSECTOR_HPP_
#define CODE_SRC_LINEINTERSECTOR_HPP_

namespace nitro {

    namespace lineIntersector {

        template< typename MLINE, typename LINE >
        class LineIntersector;

    }

}

template< typename MLINE, typename LINE >
class nitro::lineIntersector::LineIntersector
{
public:
    typedef MLINE   MLine;

    typedef LINE    Line;

    struct IntersectPoint {
        int id;
        Eigen::VectorXd coord;
        double paramLine;
        std::vector< Eigen::VectorXd > paramGeom;
    };

public:
    LineIntersector(MLine movingLines, Line line, std::string method, double tol) {
        movingLines_ = movingLines;
        line_ = line;
        method_ = method;
        tol_ = tol;
    }

    void getConstMatrix();

    void solveEigenValues() {

        this->getConstMatrix();

        if (method_ == "REDUCTION") {
            std::cout << "use REDUCTION method" << std::endl << std::endl;
        	this->solverReduction();
        } else if (method_ == "SQUARE") {
            std::cout << "use SQUARE method" << std::endl << std::endl;
        	this->solverSquare();
        } else if (method_ == "SVD") {
            std::cout << "use SVD method" << std::endl << std::endl;
        	this->solverSvd();
        } else if (method_ == "QR") {
            std::cout << "use QR method" << std::endl << std::endl;
        	this->solverQr();
        } else {
            std::cerr << "ERROR: method not available!" << std::endl;
        	assert(false);
        }
    }

    // Pencil reduction method
    void solverReduction();

    // Largest square matrix method
    void solverSquare();

    // square matrix extracted via SVD
    void solverSvd();
	
    // square matrix extracted via QR
    void solverQr();

    void computeIntersection(std::string inversionMethod)
    {
        if (inversionMethod == "SVD") {
            std::cout << "inversion using SVD" << std::endl;
            this->computeIntersectionSVD();
        }
        else if (inversionMethod == "EIGENVEC") {
            std::cout << "inversion using eigenvectors" << std::endl;
            this->computeIntersectionEigenVec();
        }
    }
	
    // inversion process with SVD
    void computeIntersectionSVD();

    // inversion process with eigenvectors
    void computeIntersectionEigenVec();

    void giveIntersection( std::vector<IntersectPoint>& intersectPoints );

public:
    void printConstMatrixA() { std::cout << matrixA_ << std::endl; }
    void printConstMatrixB() { std::cout << matrixB_ << std::endl; }

private:
    struct find_cmp {
        find_cmp(double v, double epsilon) {
            val_ = v;;
            epsilon_ = epsilon;
        }
        inline bool operator()(const double &x) const {
            return std::abs( x - val_ ) < epsilon_;
        }
    private:
        double val_;
        double epsilon_;
    };

private:
    MLine   movingLines_;

    Line    line_;

    std::string     method_; // method to solve the generalised eigenvalue problem

    double  tol_; // tolerance for numerical comparison

    Eigen::MatrixXd     matrixA_, matrixB_;

    std::vector<double>     eigenvalues_; // all intersection eigenvalues (including fictitious ones if "SQUARE" method is selected)

    std::vector< Eigen::MatrixXd >  eigenvectors_;

    std::vector< IntersectPoint >   intersectPoints_;

};

#include "LineIntersector.ipp"

#endif /* CODE_SRC_LINEINTERSECTOR_HPP_ */
