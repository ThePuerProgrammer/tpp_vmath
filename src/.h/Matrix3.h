/**
 * @file                Matrix3.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               A 3x3 matrix that inherits from Matrix.h
 * @version             0.1
 * @date                2021-09-28
 * @copyright           Copyright (c) 2021
 */

#include "Matrix.h"

namespace TPP_VMath
{
    class Matrix3 : public Matrix
    {
    public:

        /**
         * @brief       Default constructor for the zero 3x3 matrix 
         * @param       vects a std::vector<Vect*> representing the columns of a
         *              3x3 matrix
         */
        Matrix3(std::vector<Vect*>);
    };
}