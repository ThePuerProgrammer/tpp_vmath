/**
 * @file                TPP_VMath_Exception.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/TPP_VMath_Exception.h"
#include <string>

namespace TPP_VMath
{
    #pragma region TPP_VMath_Exception_Implementation
    //========================================================================//
    // START OF TPP_VMATH_EXCEPTION CLASS IMPLEMENTATION
    //========================================================================//

    TPP_VMath_Exception::TPP_VMath_Exception(std::string msg, int num)
    : errCode(num), errMsg(msg), std::runtime_error(msg)
    {   }

    int     TPP_VMath_Exception::get_error_code() const
    {
        return errCode;
    }

    //========================================================================//
    // END OF TPP_VMATH_EXCEPTION CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion TPP_VMath_Exception_Implementation
}