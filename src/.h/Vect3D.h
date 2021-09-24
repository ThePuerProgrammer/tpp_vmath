/**
 * @file                Vect3D.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               A floating point vector in R3 that inherits from Vect
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VECT3D_H
#define VECT3D_H

#include "VectND.h"
#include <array>

namespace TPP_VMath
{
    #pragma region Vect3D_Declaration
    //========================================================================//
    // START OF VECT3D CLASS DECLARATION
    //========================================================================//

    class Vect3D : public VectND
    {
    public:
    
        /**
         * @brief       Constructor for a vector in R3 that inits to zero
         * @param       void
         */ 
        Vect3D();

        /**
         * @brief       Constructor for a vector in R3 that inits to (x,y,z)
         * @param       a an array of floats representing the value(s) x,y,z
         */
        Vect3D(std::array<float, 3>);

        /**
         * @brief       Copy Constructor for a vector in R3
         * @param       original a const Vect3D reference
         */ 
        Vect3D(const Vect3D&);

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect3D(); 
         * @param       components the entries in the 3D vector
         */
        Vect3D(float*);

        /**
         * @brief       The cross product of two Vect3Ds is a Vect3D that is
         *              right angled to both.          
         * @param       that a Vect3D
         * @return      A Vect3D that is the cross product of this and that
         */
        Vect3D          cross_product(Vect3D&);

        /**
         * @brief       The cross product of two Vect3Ds is a Vect3D that is
         *              right angled to both.
         * @param       a a Vect3D          
         * @param       b a Vect3D
         * @return      A Vect3D that is the cross product of a and b
         */
        static Vect3D   cross_product(Vect3D&, Vect3D&);
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect3D_Declaration
}

#endif