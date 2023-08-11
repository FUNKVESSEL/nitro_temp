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

#ifndef CODE_SRC_PROPERTIESPARSER_HPP_
#define CODE_SRC_PROPERTIESPARSER_HPP_

#include <map>
#include <iostream>

namespace nitro {

    class Mutator;

    template<typename T>
    class TypedMutator;

    class PropertiesParser;

}

class nitro::Mutator
{
public:
    virtual void read(std::istream& inp) = 0;

    virtual ~Mutator( ) { }

};

template<typename T>
class nitro::TypedMutator : public nitro::Mutator
{
public:
    TypedMutator(T& v) : v_(v) {}

    void read(std::istream& inp) { inp >> v_; }

private:
    T&	v_;
};

class nitro::PropertiesParser
{
public:
    PropertiesParser( ) { }

public:
    template<typename T>
    void registerPropertiesVar(const std::string& name, T& t)
    {
        // If name is already in the table, create a warning message
        if (table_.find(name) != table_.end()) {
            std::cerr << "Duplicate variable name!" << std::endl;
        }
        else {
            nitro::TypedMutator<T>* v = new nitro::TypedMutator<T>(t);
            table_[name] = v;
        }
    }

    void readValues(std::istream& inp)
    {
        while (inp) {
            inp >> std::ws;

            if (inp) {
                // get name of variable
                std::string name;
                inp >> name;

                if (table_.find(name) != table_.end()) {
                    table_[name] -> read(inp);
                }
                else if (!name.empty() && name != "\n") {
                    std::cerr << "Variable not recognised: " << name << std::endl;
                }

            }

        }
    }

private:
    typedef std::map< std::string, nitro::Mutator* > TableType_;

    TableType_	table_;

};

#endif /* CODE_SRC_PROPERTIESPARSER_HPP_ */
