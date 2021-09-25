/**
 * @file                Vect.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/Vect.h"
#include "../.h/Vect1D.h"
#include "../.h/Vect2D.h"
#include "../.h/Vect3D.h"
#include "../.h/Vect4D.h"
#include "../.h/VectND.h"
#include <cmath>
#include <string>
#include <iostream>

namespace TPP_VMath
{
    #pragma region Vect_Implementation
    //========================================================================//
    // START OF VECT CLASS IMPLEMENTATION
    //========================================================================//

    Vect::~Vect()
    {   }

    float   Vect::get_magnitude()
    {
        float sumOfSquares = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfSquares += components[i] * components[i];
        }

        return std::sqrt(sumOfSquares);
    }

    float   Vect::get_magnitude(Vect& vect)
    {
        float sumOfSquares = 0;
        float* components = vect.get_components(); 
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            sumOfSquares += components[i] * components[i];
        }

        return std::sqrt(sumOfSquares);
    }

    void    Vect::scale_by(int c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    void    Vect::scale_by(int c, Vect& vect)
    {
        float* components = vect.get_components();
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    void    Vect::scale_by(float c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    void    Vect::scale_by(float c, Vect& vect)
    {
        float* components = vect.get_components();
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    float   Vect::dot_product(const Vect& that)
    {
        if (this->dimension != that.get_dimension())
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "The dot product of a vector in R";
            error += std::to_string(this->dimension) + " ";
            error += "and a vector in R";
            error += std::to_string(that.dimension) + " ";
            error += "is undefined.";
            throw std::runtime_error(error);
        }

        float sumOfProducts = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += this->components[i] * that.get_components()[i];
        }

        return sumOfProducts;
    }

    float   Vect::dot_product(const Vect& a, const Vect& b)
    {
        int dimension = a.get_dimension();

        if (dimension != b.get_dimension())
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "The dot product of a vector in R";
            error += std::to_string(a.dimension) + " ";
            error += "and a vector in R";
            error += std::to_string(b.dimension) + " ";
            error += "is undefined.";
            throw std::runtime_error(error);
        }

        float sumOfProducts = 0;
        float* acomponents = a.get_components();
        float* bcomponents = b.get_components();

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += acomponents[i] * bcomponents[i];
        }

        return sumOfProducts;
    }

    void    Vect::normalize()
    {
        float magnitude = this->get_magnitude();

        if (magnitude == 0)
        {
            std::string error = "TPP_VMath_Exception: Divide by zero ";
            error += "while attempting to normalize Vect";
            throw std::runtime_error(error);
        }

        float inverseMagnitude = 1 / magnitude;

        this->scale_by(inverseMagnitude);
    }

    void    Vect::normalize(Vect& vect)
    {
        float magnitude = vect.get_magnitude();

        if (magnitude == 0)
        {
            std::string error = "TPP_VMath_Exception: Divide by zero ";
            error += "while attempting to normalize Vect";
            throw std::runtime_error(error);
        }

        float inverseMagnitude = 1 / magnitude;

        vect.scale_by(inverseMagnitude);
    }

    void    Vect::vect_addition(Vect& vect)
    {
        if (this->dimension != vect.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw std::runtime_error(error);
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->components[i] += vect.components[i];
        }
    }

    void    Vect::vect_addition(Vect& a, Vect& b)
    {
        if (a.dimension != b.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw std::runtime_error(error);
        }

        for (int i = 0; i < a.dimension; ++i)
        {
            a.components[i] += b.components[i];
        }
    }

    void    Vect::vect_subtraction(Vect& vect)
    {
        if (this->dimension != vect.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw std::runtime_error(error);
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->components[i] -= vect.components[i];
        }
    }

    void    Vect::vect_subtraction(Vect& a, Vect& b)
    {
        if (a.dimension != b.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw std::runtime_error(error);
        }

        for (int i = 0; i < a.dimension; ++i)
        {
            a.components[i] -= b.components[i];
        }
    }

    void    Vect::set_components(std::vector<float> v)
    {
        if (dimension != 0 && v.size() != dimension)
        {
            std::string error = "TPP_VMath_Exception: Provided entries of ";
            error += "set_components(std::vector<float>) do not match the ";
            error += "dimensions of the Vect from which is was called.";
            throw std::runtime_error(error);
        }

        if (dimension == 0)
        {
            dimension = v.size();
            components = new float[v.size()];
        }

        for (int i = 0; i < dimension; ++i)
        {
            components[i] = v[i];
        }
    }

    float*  Vect::operator[](int i)
    {
        return &components[i];
    }

    void    Vect::operator*=(float c)
    {
        this->scale_by(c);
    }

    void    Vect::operator/=(float c)
    {
        if (c == 0)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "operator/=(float c) divide by zero.";
            throw std::runtime_error(error);
        }

        this->scale_by(1/c);
    }

    void    Vect::operator+=(Vect& right)
    {
        this->vect_addition(right);
    }

    void    Vect::operator-=(Vect& right)
    {
        this->vect_subtraction(right);
    }

    Vect*   Vect::operator+(Vect& right)
    {
        Vect* result = get_vect_of_valid_dimensions(this);
        result->vect_addition(right);
        return result;
    }

    Vect*   Vect::operator-(Vect& right)
    {
        Vect* result = get_vect_of_valid_dimensions(this);
        result->vect_subtraction(right);
        return result;
    }

    Vect&   Vect::operator=(const Vect& right)
    {
        if (this->dimension != right.get_dimension())
        {
            std::string error = "TPP_VMath_Exception: Dimensions of ";
            error += "operator=(const Vect& right) do not match the ";
            error += "dimensions of the left side Vect.";
            throw std::runtime_error(error);
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->components[i] = right.get_components()[i];
        }

        return *this;
    }

    Vect*   Vect::get_vect_of_valid_dimensions(Vect* vect)
    {
        float* c = vect->get_components();
        int    m = vect->get_dimension();

        switch (m)
        {
            case  1: return new Vect1D(c);
            case  2: return new Vect2D(c);
            case  3: return new Vect3D(c);
            case  4: return new Vect4D(c);
            default: return new VectND(m, c);
        }
    }

    float*  Vect::get_components() const
    {
        return components;
    }

    unsigned int Vect::get_dimension() const
    {
        return dimension;
    }

    void    Vect::print()
    {
        std::cout << "------\n";
        for (int i = 0; i < dimension; ++i)
        {
            std::cout << "|" << components[i] << "|\n";
        }
        std::cout << "------\n";
    }

    //========================================================================//
    // END OF VECT CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect_Implementation
}