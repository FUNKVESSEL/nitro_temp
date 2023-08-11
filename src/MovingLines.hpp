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

#ifndef CODE_SRC_MOVINGLINES_HPP_
#define CODE_SRC_MOVINGLINES_HPP_

namespace nitro {

    namespace movingLines {

        template<int DIM, typename GEOM>
        class MovingLines;

        // moving lines for curves
        template<typename GEOM>
        class MovingLines<1, GEOM>;

        // moving lines (planes) for surfaces
        template<typename GEOM>
        class MovingLines<2, GEOM>;

    }

}

template< typename GEOM >
class nitro::movingLines::MovingLines< 1, GEOM >
{
public:
    typedef GEOM    CurveMonomial;

    typedef typename CurveMonomial::Point   Point;

public:
    MovingLines( ) {
        nu_ = 0; // Initialisation of auxiliary degree;
        isSpecified_ = false;

        tolRank_ = 1.e-6;
    }

    MovingLines( int nu ) {
        nu_ = nu;
        isSpecified_ = true;

        tolRank_ = 1.e-6;

        if (nu_ == 0) nu_ = 1; // for the sake of inversion process
    }

    void registerCurve(GEOM curveMonomial) { curveMonomial_ = curveMonomial; } ;

    void setRankTolerance(double tol) { tolRank_ = tol; }

    void computeCoefficients();

    Eigen::MatrixXd giveCoefficients() { return matrixCoeffsNull_; }

    void computeImplicitMatrix( Eigen::VectorXd coord );

    Eigen::MatrixXd giveImplicitMatrix() { return matrixImplicit_; }

    Eigen::VectorXd evaluateGeometry( Eigen::VectorXd param ) {

        double theta = param(0);

        return curveMonomial_.evaluate(theta);

    }

    void inversionSVD( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol );

    void inversionEigenVec( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol );

public:
    int evaluateRank( Eigen::VectorXd singularValues );

    bool checkDegenerate( Eigen::MatrixXd matrixA, Eigen::MatrixXd matrixB ); // check if A - x.B is always rank-deficient

    int numRowsImplicit() { return nu_ + 1; }

    int numColsImplicit() { return matrixCoeffsNull_.cols(); }

    int spaceDim() { return curveMonomial_.spaceDim(); }

    void printMatrixCoeffs() { std::cout << matrixCoeffs_ << std::endl; }

private:
    GEOM    curveMonomial_;

    int     nu_; // the auxiliary degree

    bool    isSpecified_; // if or not the auxiliary degree is specified

    Eigen::MatrixXd     matrixCoeffs_;

    Eigen::MatrixXd     matrixCoeffsNull_;

    Eigen::MatrixXd     matrixImplicit_;

    double tolRank_; // tolerance for rank estimation
};

template< typename GEOM >
class nitro::movingLines::MovingLines< 2, GEOM >
{
public:
    typedef GEOM    SurfaceMonomial;

    typedef typename SurfaceMonomial::Point     Point;

public:
    MovingLines( ) {
        nu1_ = 0; // Initialisation of auxiliary degrees
        nu2_ = 0;
        spaceDim_ = 3;
        isSpecified_ = false;

        tolRank_ = 1.e-6;
    }

    MovingLines(int nu1, int nu2) {
        nu1_ = nu1;
        nu2_ = nu2;
        spaceDim_ = 3;
        isSpecified_ = true;

        tolRank_ = 1.e-6;

        if (nu1_ == 0) nu1_ = 1; // for the sake of inversion process
        if (nu2_ == 0) nu2_ = 1; // for the sake of inversion process
    }

    void registerSurface(GEOM surfaceMonomial) { surfaceMonomial_ = surfaceMonomial; }

    void setRankTolerance(double tol) { tolRank_ = tol; }

    void giveSurfaceDegree() { degree_ = surfaceMonomial_.degree(); }

    void setIndexMapAux();

    void computeCoefficients();

    Eigen::MatrixXd giveCoefficients() { return matrixCoeffsNull_; }

    Eigen::VectorXd evaluateGeometry( Eigen::VectorXd param ) { return surfaceMonomial_.evaluate(param); }

    void inversionSVD( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol );

    void inversionEigenVec( Eigen::MatrixXd matrix, std::vector< Eigen::VectorXd >& params, double tol );

public:
    int evaluateRank( Eigen::VectorXd singularValues );

    bool checkDegenerate( Eigen::MatrixXd matrixA, Eigen::MatrixXd matrixB ); // check if A - x.B is always rank-deficient

    int numRowsImplicit() { return (nu1_ + 1)*(nu2_ + 1); }

    int numColsImplicit() { return matrixCoeffsNull_.cols(); }

    int spaceDim() { return spaceDim_; }

    void printMatrixCoeffs() { std::cout << matrixCoeffs_ << std::endl; }

private:
    struct unique_cmp {
        unique_cmp(double epsilon) {
            epsilon_ = epsilon;
        }
        inline bool operator()(const double &x, const double &y) const {
            return std::abs(x - y) < epsilon_;
        }

    private:
        double epsilon_;

    };

private:
    std::pair<int, int> degree_; // Surface degree

    int nu1_, nu2_; // Auxiliary bi-degree

    int spaceDim_;

    bool isSpecified_; // whether or not the auxiliary degee is specified

    GEOM surfaceMonomial_;

    std::vector< std::pair<int, int> > indexMapAux_;

    Eigen::MatrixXd matrixCoeffs_;

    Eigen::MatrixXd matrixCoeffsNull_;

    double tolRank_; // tolerance for rank estimation

};

#include "MovingLinesCurve.ipp"
#include "MovingLinesSurface.ipp"

#endif /* CODE_SRC_MOVINGLINES_HPP_ */
