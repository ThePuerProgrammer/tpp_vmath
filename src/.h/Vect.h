/**
 * @file                Vect.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               Parent class for all vectors in the library. Children 
 *                      that inherit from Vect should be thought of as having a 
 *                      vertical orientation when considering their relationship
 *                      to matrices. 
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VECT_H
#define VECT_H

#include <vector>

namespace TPP_VMath
{
    #pragma region Vect_Declaration
    //========================================================================//
    // START OF VECT CLASS DECLARATION
    //========================================================================//

    class Vect
    {
    public:

        /**
         * @brief       The magnitude of a Vect is its distance from the origin
         * @param       void
         * @return      std::sqrt(v₁^2+...+vᵢ^2) as a float
         */
        float           get_magnitude();

        /**
         * @brief       The magnitude of a Vect is its distance from the origin
         * @param       vect a reference to a child of Vect 
         * @return      std::sqrt(v₁^2+...+vᵢ^2) as a float
         */
        static float    get_magnitude(Vect&);

        /**
         * @brief       Scales each entry in the vector by integer scalar value
         * @param       c an integer scalar value
         */ 
        void            scale_by(int);

        /**
         * @brief       Scales each entry in the vector by integer scalar value
         * @param       c an integer scalar value
         * @param       vect a reference to a child of Vect
         */ 
        static void     scale_by(int, Vect&);

        /**
         * @brief       Scales each entry in the vector by float scalar value
         * @param       c a floating point scalar value
         */ 
        void            scale_by(float);

        /**
         * @brief       Scales each entry in the vector by float scalar value
         * @param       c an float scalar value
         * @param       vect a reference to a child of Vect
         */ 
        static void     scale_by(float, Vect&);

        /**
         * @brief       Faulty dimensions throws TPP_VMath_Exception
         * @param       that a const Vect reference
         * @return      sum of the entrywise products of two Vects
         */ 
        float           dot_product(const Vect&);

        /**
         * @brief       Faulty dimensions throws TPP_VMath_Excpetion
         * @param       a a const Vect reference
         * @param       b a const Vect reference
         * @return      sum of the entrywise products of two Vects
         */ 
        static float    dot_product(const Vect&, const Vect&);

        /**
         * @brief       Normalization scales the vector by the inverse of its
         *              magnitude, which has the effect of setting the magnitude
         *              of the vector to 1. If called upon the zero vector, a
         *              TPP_VMath_Exception is thrown to prevent divide by zero
         * @param       void
         */
        void            normalize(); 

        /**
         * @brief       Normalization scales the vector by the inverse of its
         *              magnitude, which has the effect of setting the magnitude
         *              of the vector to 1. If called upon the zero vector, a
         *              TPP_VMath_Exception is thrown to prevent divide by zero
         * @param       vect a reference to a child of Vect
         */
        static void     normalize(Vect&); 

        /**
         * @brief       Component-wise addition of two Vects
         * @param       vect a Vect reference of matching dimensions 
         */
        void            vect_addition(Vect&);

        /**
         * @brief       Component-wise addition of two Vects of matching
         *              dimensions.
         * @param       a a Vect reference
         * @param       b a Vect reference
         */
        static void     vect_addition(Vect&, Vect&);

        /**
         * @brief       Component wise subtraction of two Vects
         * @param       vect a Vect reference of matching dimensions
         */
        void            vect_subtraction(Vect&);

        /**
         * @brief       Component wise subtraction of two Vects of matching
         *              dimensions
         * @param       a a Vect reference
         * @param       b a Vect reference
         */
        static void     vect_subtraction(Vect&, Vect&);

        /**
         * @brief       Manual assignment of entires in the vector. Faulty
         *              dimensions throws TPP_VMath_Exception
         * @param       v the provided array of floats
         */ 
        void            set_components(std::vector<float>);

        /**
         * @param       i the ith coordinate
         * @return      a pointer to the ith coordinate in the vector
         */
        float*          operator[](int);

        /**
         * @brief       Overloaded *= calls the appropriate scale_by() function.
         *              Faulty dimensions throws TPP_VMath_Exception
         * @param       c a scalar
         */
        void            operator*=(float c); 

        /**
         * @brief       Overloaded /= calls the appropriate scale_by() function,
         *              passing in the inverse of the right side argument.
         *              Faulty dimensions throws TPP_VMath_Exception
         *              Zero throws TPP_VMath_Exception to prevent divide by 0
         * @param       c a scalar
         */
        void            operator/=(float c); 

        /**
         * @brief       Component wise addition of two Vects of matching
         *              dimensions.
         * @param       right a Vect reference
         */
        void            operator+=(Vect&);

        /**
         * @brief       Component wise subtraction of two Vects of matching
         *              dimensions.
         * @param       right a Vect reference
         */
        void            operator-=(Vect&);

        /**
         * @brief       Component wise addition of two Vects of matching
         *              dimensions. Result should be wrapped with a VWrap object
         *              in order to prevent memory leaks. 
         * @param       right a Vect reference
         * @return      Vect* to new Vect child of matching dimensions that is
         *              the sum of the left and right Vect.
         */
        Vect*           operator+(Vect&);

        /**
         * @brief       Component wise subtraction of two Vects of matching
         *              dimensions. Result should be wrapped with a VWrap object
         *              in order to prevent memory leaks.
         * @param       right a Vect reference
         * return       Vect* to a new Vect child of matching dimensions that
         *              is the difference between the left and right Vect.
         */
        Vect*           operator-(Vect&);

        /**
         * @brief       Copies components from right vector to left vector. 
         *              Mismatched dimensions throws TPP_VMath_Exception
         * @param       right a const Vect reference to the right side of =
         */ 
        Vect&           operator=(const Vect&);

        /**
         * @brief       Polymorphic resolution for copy construction
         * @param       vect the Vect being copied
         * @return      a Vect* to a new child vect of m dimensions
         */ 
        static Vect*    get_vect_of_valid_dimensions(Vect*);

        /**
         * @param       void
         * @return      the coordinate array as a pointer to components[0]
         */
        float*          get_components() const;

        /**
         * @param       void
         * @return      the dimension of the vector (the number of entries)
         */ 
        unsigned int    get_dimension() const;

        /**
         * @brief       A print function for simple tests in the console
         * @param       void
         */
        void            print();

        /**
         * @brief       Pure virtual destructor
         */ 
        virtual         ~Vect() = 0;

    protected:

        ///             The number of entries in components
        unsigned int    dimension;

        ///             The array representing the vector
        float*          components;
    };

    //========================================================================//
    // END OF VECT CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect_Declaration
}

#endif