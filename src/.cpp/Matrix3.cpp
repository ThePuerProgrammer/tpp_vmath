/**
 * @file                Matrix3.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               Implementation file for a 3x3 matrix
 * @version             0.1
 * @date                2021-09-28
 * @copyright           Copyright (c) 2021
 */

#include "../.h/Matrix3.h"

namespace TPP_VMath
{
    #pragma region Matrix3_Implementation
    //========================================================================//
    // START OF MATRIX3 CLASS IMPLEMENTATION
    //========================================================================//

    Matrix3::Matrix3(std::vector<Vect*> vects) : Matrix(vects)
    {   
        if (vects.size() != 3) throw std::runtime_error("Matrix3 must be 3x3");

        for (int i = 0; i < 3; ++i)
        {
            if (vects[i]->get_dimension() != 3)
            {
                throw std::runtime_error("Matrix3 must be 3x3");
            }
        }
    }

    //========================================================================//
    // END OF MATRIX3 CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Matrix3_Implementation
}