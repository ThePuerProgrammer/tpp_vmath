/**
 * @file                TPP_VMath_Exception.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               Class is used to handle runtime errors associated with
 *                      undefined vector and matrix operations, such as 
 *                      transformations outside of the defined range of T
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef TPP_VMATH_EXCEPTION_H
#define TPP_VMATH_EXCEPTION_H

#include <string>
#include <stdexcept>

namespace TPP_VMath
{
    #pragma region TPP_VMath_Exception_Declaration
    //========================================================================//
    // START OF TPP_VMATH_EXCEPTION CLASS DECLARATION
    //========================================================================//

    class TPP_VMath_Exception : public std::runtime_error
    {
    public:
        /** 
         *  @brief      Exception handler for runtime errors associated with
         *              undefined vector and matrix operations.
         *  @param      msg a description of the exception thrown
         *  @param      num a value attributed to the error type
         */
        TPP_VMath_Exception(std::string, int);

        /**
         * @brief       return the integer value of the error code
         */
        int             get_error_code() const; 

    private:

        ///             the error code
        int             errCode;     

        ///             a description of the exception, e.what()
        std::string     errMsg;     
    };

    //========================================================================//
    // START OF TPP_VMATH_EXCEPTION CLASS DECLARATION
    //========================================================================//
    #pragma endregion TPP_VMath_Exception_Declaration
}

#endif