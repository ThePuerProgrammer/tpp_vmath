/**
 * @file                Vect1D.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               A one dimensional vector can be thought of as a scalar 
 *                      but there are some unique edge cases where this is a
 *                      very interesting and useful property. A floating point
 *                      vector in R1 that inherits from VectND
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VECT1D_H
#define VECT1D_H

#include "VectND.h"
#include <array>

namespace TPP_VMath
{
    #pragma region Vect1D_Declaration
    //========================================================================//
    // START OF VECT1D CLASS DECLARATION
    //========================================================================//

    class Vect1D : public VectND
    {
    public:

        /**
         * @brief       constructor for a vector in R1 that inits to zero
         * @param       void
         */ 
        Vect1D();

        /**
         * @brief       constructor for a vector in R1 that inits to (x)
         * @param       a an array of floats representing the value(s) x
         */ 
        Vect1D(std::array<float, 1>);

        /**
         * @brief       Copy constructor for a vector in R1
         * @param       original a const Vect1D reference
         */
        Vect1D(const Vect1D&); 

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect1D(); 
         * @param       components the entries in the 1D vector
         */
        Vect1D(float*);
    }; 
    //========================================================================//
    // END OF VECT1D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect1D_Declaration
}

#endif