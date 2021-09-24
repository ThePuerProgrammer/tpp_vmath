/**
 * @file                Vect1D.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/VectND.h"
#include "../.h/Vect1D.h"

namespace TPP_VMath
{
    #pragma region Vect1D_Implementation
    //========================================================================//
    // START OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//

    Vect1D::Vect1D() : VectND()
    {   }

    Vect1D::Vect1D(std::array<float, 1> a) : VectND(1, a.begin())
    {   }

    Vect1D::Vect1D(const Vect1D& original) : VectND(1, original.components)
    {   }

    Vect1D::Vect1D(float* components) : VectND(1, components)
    {   }

    //========================================================================//
    // END OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect1D_Implementation
}