/**
 * @file                Vect4D.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               A floating point vector in R4 that inherits from Vect
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VECT4D_H
#define VECT4D_H

#include "VectND.h"
#include <array>

namespace TPP_VMath
{
    #pragma region Vect4D_Declaration
    //========================================================================//
    // START OF VECT4D CLASS DECLARATION
    //========================================================================//

    class Vect4D : public VectND
    {
    public:

        /**
         * @brief       A vector in R4 that defaults to zero
         * @param       void
         */
        Vect4D(); 

        /**
         * @brief       A vector in R4 that inits to (x, y, z, w)
         * @param       a an array of floats representing the value(s) x,y,z,w
         */
        Vect4D(std::array<float, 4>);

        /**
         * @brief       Copy Constructor for a vector in R4
         * @param       original a const Vect4D reference
         */
        Vect4D(const Vect4D&);

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect4D(); 
         * @param       components the entries in the 4D vector
         */
        Vect4D(float*); 
    };

    //========================================================================//
    // END OF VECT4D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect4D_Declaration
}

#endif