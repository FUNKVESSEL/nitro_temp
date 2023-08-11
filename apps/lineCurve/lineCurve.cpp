// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
//! This file is part of nitro.
//! @Authors: Xiao Xiao, Laurent Busé and Fehmi Cirak
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

#include <eigen/Eigen/Dense>

#include "../../src/PropertiesParser.hpp"
#include "../../src/Input.hpp"
#include "../../src/Curve.hpp"
#include "../../src/MovingLines.hpp"
#include "../../src/LineIntersector.hpp"

#include "Parameters.hpp"

#define LAGRANGE_

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: lineCurve input.txt" << std::endl;
        return 0;
    }

    //==================================================================================
    // read input parameters
    const std::string inputData(argv[1]);
    std::ifstream inp(inputData.c_str());

    if (inp.is_open() == false) std::cerr << "Cannot open input.txt" << std::endl;

    lineCurve::Parameters parameters(inp);
    inp.close();

    std::ifstream curvef(parameters.geometryFileName.c_str());

    std::string basisType = parameters.basisType;

    const int spaceDim = parameters.spaceDim;

    // read geometry
    typedef nitro::input::Point CurvePoint;

    std::vector<CurvePoint> curvePoints;

    nitro::input::readGeometry(curvef, curvePoints, spaceDim);

    int numCurvePoints = curvePoints.size();
    int degree = numCurvePoints - 1;

    std::cout << "number of curve points: " << numCurvePoints << std::endl << std::endl;

    // create curve
#if defined LAGRANGE_
    const nitro::curve curveBasis = nitro::LAGRANGE;
#elif defined MONOMIAL_
    const lineCurve::curve curveBasis = lineCurve::MONOMIAL;
#endif

    typedef nitro::Curve< curveBasis, CurvePoint > Curve;

    Curve curve(spaceDim);

    curve.setCoefficients(curvePoints);

    std::vector<CurvePoint> coeffsMonomial;
    curve.getMonomialCoeffs(coeffsMonomial);

    // convert to monomial curve
    typedef nitro::Curve< nitro::MONOMIAL, CurvePoint > CurveMonomial;

    CurveMonomial curveMonomial(spaceDim);

    curveMonomial.setCoefficients(coeffsMonomial);

    //==================================================================================
    // compute coefficient matrix
    typedef nitro::movingLines::MovingLines< Curve::dim, CurveMonomial > MovingLines;

    MovingLines movingLines;

    movingLines.registerCurve(curveMonomial);

    movingLines.setRankTolerance(1.e-6); // tolerance for numerical rank estimation; the default is 1.e-6.

    movingLines.computeCoefficients();

    //==================================================================================
    // create a line
    std::ifstream linef(parameters.lineFileName.c_str());

    typedef nitro::input::Line Line;

    Line line;

    nitro::input::readLine(linef, line, spaceDim);

    // compute line / curve intersection
    typedef nitro::lineIntersector::LineIntersector< MovingLines, Line > LineIntersector;

    std::string method = parameters.method; // method to solve the generalised eigenvalue problem
    const double tol = 1.e-8; // tolerance for numerical comparison in intersection computation

    LineIntersector lineIntersector(movingLines, line, method, tol);

    lineIntersector.solveEigenValues();

    std::string inversion = parameters.inversion;

    lineIntersector.computeIntersection(inversion);

    //==================================================================================
    // output intersection points
    typedef LineIntersector::IntersectPoint IntersectPoint;

    std::vector<IntersectPoint> intersectPoints;

    lineIntersector.giveIntersection(intersectPoints);

    Eigen::IOFormat outFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " [", "]");

    std::cout << "\n" << "number of intersection points: " << intersectPoints.size() << std::endl;
    for (int i = 0; i < intersectPoints.size(); ++i) {
        IntersectPoint ip = intersectPoints[i];

        std::cout << "***** point " << i + 1 << " *****" << std::endl;

        std::cout << "coordinate: " << std::endl;
        std::cout << ip.coord.format(outFormat) << std::endl;

        std::cout << "line parameter: " << std::endl;
        std::cout << " " << ip.paramLine << std::endl;

        std::cout << "curve parameter: ";

        if (ip.paramGeom.size() > 1) std::cout << "(multiple pre-images)";

        std::cout << std::endl;

        for (int j = 0; j < ip.paramGeom.size(); ++j) {
            std::cout << " " << ip.paramGeom[j] << std::endl;
        }

        std::cout << "verify: " << std::endl;
        for (int j = 0; j < ip.paramGeom.size(); ++j) {
            double theta = ip.paramGeom[j](0);
            Eigen::VectorXd coordTheta = curve.evaluate(theta);
            std::cout << coordTheta.format(outFormat) << std::endl;
        }

        std::cout << std::endl;

    }

    return 0;
}

