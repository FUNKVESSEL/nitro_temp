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
#include "../../src/Surface.hpp"
#include "../../src/MovingLines.hpp"
#include "../../src/LineIntersector.hpp"

#include "Parameters.hpp"

#define TENSORLAGRANGE_

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: lineSurface input.txt" << std::endl;
        return 0;
    }

    //================================================================================
    // read input parameters
    const std::string inputData(argv[1]);
    std::ifstream inp(inputData.c_str());

    if (inp.is_open() == false) std::cerr << "Cannot open input.txt" << std::endl;

    lineSurface::Parameters parameters(inp);
    inp.close();

    std::ifstream surfacef(parameters.geometryFileName.c_str());

    std::string basisType = parameters.basisType;

    const int spaceDim = parameters.spaceDim;
    const int uDeg = parameters.uDeg, vDeg = parameters.vDeg;

    // read geometry
    typedef nitro::input::Point SurfPoint;

    std::vector<SurfPoint> surfPoints;

    nitro::input::readGeometry(surfacef, surfPoints, spaceDim);

    int numSurfPoints = surfPoints.size();

    if ((uDeg + 1)*(vDeg + 1) != numSurfPoints) std::cerr << "Surface points not compatible with degrees!" << std::endl;

    std::cout << "number of surface points: " << numSurfPoints << std::endl << std::endl;

    // create surface
#if defined TENSORLAGRANGE_
    const nitro::surface surfaceBasis = nitro::TENSORLAGRANGE;
#elif defined TENSORMONOMIAL_
    const lineSurface::surface surfaceBasis = lineSurface::TENSORMONOMIAL;
#endif

    typedef nitro::Surface< surfaceBasis, SurfPoint > Surface;

    Surface surface(uDeg, vDeg);

    surface.setCoefficients(surfPoints);

    std::vector<SurfPoint> coeffsMonomial;
    surface.getMonomialcoeffs(coeffsMonomial);

    // convert to monomial surface
    typedef nitro::Surface< nitro::TENSORMONOMIAL, SurfPoint > SurfaceMonomial;

    SurfaceMonomial surfaceMonomial(uDeg, vDeg);

    surfaceMonomial.setCoefficients(coeffsMonomial);

    //==============================================================================================
    // compute coefficient matrix
    typedef nitro::movingLines::MovingLines< Surface::dim, SurfaceMonomial > MovingLines;

    MovingLines movingLines;

    movingLines.registerSurface(surfaceMonomial);

    movingLines.setRankTolerance(1.e-6); // tolerance for numerical rank estimation; the default is 1.e-6.

    movingLines.computeCoefficients();

    //=============================================================================================
    // create a line
    std::ifstream linef(parameters.lineFileName.c_str());

    typedef nitro::input::Line Line;

    Line line;

    nitro::input::readLine(linef, line, spaceDim);

    // compute line / surface intersection
    typedef nitro::lineIntersector::LineIntersector< MovingLines, Line > LineIntersector;

    std::string method = parameters.method; // method to solve the generalised eigenvalue problem
    const double tol = 1.e-7; // tolerance for numerical comparison in intersection computation

    LineIntersector lineIntersector(movingLines, line, method, tol);

    lineIntersector.solveEigenValues();

    std::string inversion = parameters.inversion;

    lineIntersector.computeIntersection(inversion);

    //=============================================================================================
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

        std::cout << "surface parameter: ";

        if (ip.paramGeom.size() > 1) std::cout << "(multiple pre-images)";

        std::cout << std::endl;

        for (int j = 0; j < ip.paramGeom.size(); ++j) {
            std::cout << ip.paramGeom[j].format(outFormat) << std::endl;
        }

        std::cout << "verify: " << std::endl;
        for (int j = 0; j < ip.paramGeom.size(); ++j) {
            Eigen::VectorXd param = ip.paramGeom[j];
            Eigen::VectorXd coordParam = surface.evaluate(param);
            std::cout << coordParam.format(outFormat) << std::endl;
        }

        std::cout << std::endl;

    }

    return 0;
}

