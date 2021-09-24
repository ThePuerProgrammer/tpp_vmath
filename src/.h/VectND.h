/**
 * @file                VectND.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               A floating point vector in RN (variable length) that
 *                      inherits from Vect
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VECTND_H
#define VECTND_H

#include "Vect.h"

namespace TPP_VMath
{
    #pragma region VectND_Declaration
    //========================================================================//
    // START OF VECTND CLASS DECLARATION
    //========================================================================//

    class VectND : public Vect
    {
    public:

        /**
         * @brief       Default constructor for a vector in RN that inherits
         *              from Vect
         * @param       void
         */
        VectND();

        /**
         * @brief       Constructor for a vector in RN where N is provided that 
         *              defaults to the zero vector
         * @param       n the dimension of the vector
         */ 
        VectND(int);

        /**
         * @brief       Constructor for a vector in RN where N is provided that
         *              initializes to the provided array of floats
         * @param       n the dimension of the vector
         * @param       components the provided array of floats
         */
        VectND(int, float*);

        /**
         * @brief       Constructor for a vector in RN where N is the size of a
         *              std::vector<float>
         * @param       v a std::vector of floats
         */
        VectND(std::vector<float>); 

        /**
         * @brief       Destructor deletes the coordinate array
         */
        ~VectND();
    }; 

    //========================================================================//
    // END OF VECTND CLASS DECLARATION
    //========================================================================//
    #pragma endregion VectND_Declaration
}

#endif