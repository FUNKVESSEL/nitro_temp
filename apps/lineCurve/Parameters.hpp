// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
//! This file is part of nitro.
//! @Authors: Xiao Xiao, Laurent Busé and Fehmi Cirak
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef CODE_APPS_LINECURVEM_PARAMETERS_HPP_
#define CODE_APPS_LINECURVEM_PARAMETERS_HPP_

namespace lineCurve {

    class Parameters;

}

class lineCurve::Parameters
{
public:
    Parameters(std::istream& inp)
    {
        nitro::PropertiesParser* prop = new nitro::PropertiesParser;

        prop->registerPropertiesVar( "geometryFileName",    geometryFileName );
        prop->registerPropertiesVar( "lineFileName",        lineFileName );
        prop->registerPropertiesVar( "basisType",           basisType );
        prop->registerPropertiesVar( "method",              method );
        prop->registerPropertiesVar( "inversion",           inversion );

        prop->registerPropertiesVar( "spaceDim",            spaceDim);

        prop->readValues(inp);

        delete prop;
    }

public:
    std::string geometryFileName;
    std::string lineFileName;
    std::string basisType;
    std::string method;
    std::string inversion;

    int spaceDim;

};

#endif /* CODE_APPS_LINECURVEM_PARAMETERS_HPP_ */
