/**
 * @file                Vect2D.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               A floating point vector in R2 that inherits from Vect
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VECT2D_H
#define VECT2D_H

#include "VectND.h"
#include <array>

namespace TPP_VMath
{
    #pragma region Vect2D_Declaration
    //========================================================================//
    // START OF VECT2D CLASS DECLARATION
    //========================================================================//

    class Vect2D : public VectND
    {
    public:

        /**
         * @brief       constructor for a vector in R2 that inits to zero
         * @param       void
         */ 
        Vect2D();

        /**
         * @brief       constructor for a vector in R2 that inits to (x,y)
         * @param       a an array of floats representing the value(s) x,y
         */ 
        Vect2D(std::array<float, 2>);

        /**
         * @brief       Copy constructor for a vector in R2
         * @param       original a const Vect2D reference
         */
        Vect2D(const Vect2D&); 

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect2D(); 
         * @param       components the entries in the 2D vector
         */
        Vect2D(float*);
    }; 
    //========================================================================//
    // END OF VECT2D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect2D_Declaration
}

#endif