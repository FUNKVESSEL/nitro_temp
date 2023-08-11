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

#ifndef CODE_SRC_INPUT_HPP_
#define CODE_SRC_INPUT_HPP_

namespace nitro {

    namespace input {

        struct Point {
            int id;
            Eigen::VectorXd coord;
            Eigen::VectorXd coordW;
        };

        struct Line {
            int id;
            std::vector<Point> linePoints;
        };

        // read .obj file of the geometry
        void readGeometry(std::ifstream& inp, std::vector< Point >& points, int spaceDim);

        // read .obj file of the line
        void readLine(std::ifstream& inp, Line& line, int spaceDim);

    }

}

void nitro::input::readGeometry(std::ifstream& inp, std::vector< Point >& points, int spaceDim)
{

    std::cout << "spaceDim = " << spaceDim << std::endl;

    while (true) {
        std::string lineString;
        std::getline(inp, lineString);

        int numPoints = 0;

        if (lineString.find("v ") != std::string::npos) {
            std::istringstream lineStringStream(lineString);

            Point point;

            point.id = numPoints;

            char tag;
            lineStringStream >> tag;

            Eigen::VectorXd coord(spaceDim);
            for (int i = 0; i < spaceDim; ++i) {
                double coordVal = 0.;
                lineStringStream >> coordVal;
                coord(i) = coordVal;
            }

            point.coord = coord;
            points.push_back(point);

            numPoints += 1;

            char c = inp.peek();
            if (c != 'v') break;

        }
    }
}

void nitro::input::readLine(std::ifstream& inp, Line& line, int spaceDim)
{
    while (true) {
        std::string lineString;
        std::getline(inp, lineString);

        int numPoints = 0;

        if (lineString.find("v ") != std::string::npos) {
            std::istringstream lineStringStream(lineString);

            Point point;
            point.id = numPoints;

            char tag;
            lineStringStream >> tag;

            Eigen::VectorXd coord(spaceDim);
            for (int i = 0; i < spaceDim; ++i) {
                double coordVal = 0.;
                lineStringStream >> coordVal;
                coord(i) = coordVal;
            }

            point.coord = coord;
            line.linePoints.push_back(point);

            numPoints += 1;

            char c = inp.peek();
            if (c != 'v') break;

        }

    }
}

#endif /* CODE_SRC_INPUT_HPP_ */
