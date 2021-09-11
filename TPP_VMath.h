// TPP_VMath.h
// written by Jesse Rankins 2021
#pragma once

namespace TPP_VMath
{

    //========================================================================//
    // START OF VECT CLASS DECLARATION
    //========================================================================//

    class Vect
    {
    public:

        // Returns sqrtf(v₁^2+...+vᵢ^2)
        float           get_magnitude();

        // Scales each entry in the vector by integer scalar value
        void            scale_by(int);

        // Scales each entry in the vector by float scalar value
        void            scale_by(float);

        // Returns pointer to the ith coordinate in the vector
        float*          operator[](int);

        // Returns the coordinate array as a pointer to coordinates[0]
        float*          get_coordinates() const;

        // Returns the dimension of the vector (the number of entries)
        unsigned int    get_dimension() const;

        // Pure virtual destructor
        virtual         ~Vect() = 0;

    protected:

        // The number of entries in coordinates
        unsigned int    dimension;

        // The array representing the vector
        float*          coordinates;
    };

    //========================================================================//
    // END OF VECT CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VWRAPPER CLASS DECLARATION
    //========================================================================//

    /*
    This class should be used to wrap any function that returns a Vect* to a 
    newly allocated resource, such as the function
    
    Vect* Matrix::get_matrix_vector_product()

    Wrapping a new Vect allows the destructor to handle deallocation of heap 
    resources. If two VWrappers are used to wrap the same Vect*, this will 
    result in undefined behavior. 
    
    Therefore, a single address should only ever be wrapped once. 
    */
    class VWrapper
    {
    public:

        // Set vect to nullptr and wrap as a seperate step
        VWrapper();

        // Wrap upon construction
        VWrapper(Vect*);

        // Destructor handles deallocation
        ~VWrapper();

        // If vect is null, wrap a Vect*
        void    wrap(Vect*);

        // Explicit deallocation of vect
        void    unwrap();

        // Assignment for access of the wrapped pointer
        Vect*   get_vect();

    protected:

        // The wrapped Vect*
        Vect*   vect;
    };

    //========================================================================//
    // END OF VWRAPPER CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT1D CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECT1D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT2D CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECT2D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT3D CLASS DECLARATION
    //========================================================================//

    class Vect3D : public Vect
    {
    public:
    
        // Default constructor for a vector in R3 that inits to the zero vector
        Vect3D();

        /*
        Overloaded constructor for a vector in R3 that accepts x,y,z coordinates
        */
        Vect3D(const float, const float, const float);

        // Copy constructor
        Vect3D(const Vect3D&);

        /*
        Constructor useful for polymorphic copies of Vect* v = new Vect3D();
        Error prone!! Call should ALWAYS be wrapped with guard in the form
        if (v->get_dimension() == 3)
        {
            Vect3D(v->get_components());
        }
        */
        Vect3D(const float*);

        // Destructor deletes coordinate array
        ~Vect3D();

        // Copies coordinates from the right side Vect3D to the left side Vect3D
        Vect3D& operator=(const Vect3D&);

        // Manual assignment of x, y, z vector entries
        void    set_coordinates(float, float, float);

        // Calculates the sum of the entrywise products of two Vect3Ds
        float   dot_product(const Vect3D&);
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT4D CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECT4D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECTND CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECTND CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF SET CLASS DECLARATION
    //========================================================================//

    class VSet
    {
    public:

        // Default constructor is a null set
        VSet();

        // Overloaded constructor for one Vect3D
        VSet(const Vect3D&);

        // Overloaded constructor for an array of Vect3D
        VSet(const int&, Vect3D[]);

        // Copy constructor
        VSet(const VSet&);

        // Destructor deletes the set
        ~VSet();

        // Returns the number of vectors in the set
        unsigned int    get_n();
        
        // Returns the number of entries in each vector
        unsigned int    get_m();

        // Returns the set as an array of Vects
        Vect**          get_set_of_vectors();

        /*
        Add a new vect of *matching dimension* to the set
        If the vector dimension doesn't == m, a runtime exception is thrown
        It is best to wrap the call to this funtion in a guard in the form
        if (Vect*->get_dimension() == Set.get_m())
        {
            add_vect_to_set(Vect*);
        }
        This does not apply if the set is currently null
        */
        void            add_vect_to_set(Vect*);

    private:

        // The number of entries in each vector in the set
        unsigned int    m;

        // The number of vectors in the set
        unsigned int    n;

        // The set of vectors as an array in pointer notation
        Vect**          setOfVectors;
    };

    //========================================================================//
    // END OF SET CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF MATRIX CLASS DECLARATION
    //========================================================================//

    class Matrix
    {
    public:

        // Overloaded constructor accepts a set of Vects and generates a matrix
        Matrix(VSet&);

        // Matrix copy constructor
        Matrix(const Matrix&);

        // Destructor deletes the matrix A
        ~Matrix();

        // Reduce this matrix to reduced echelon form
        void            reduce_matrix();

        // Produce a new reduced echelon matrix from this matrix
        Matrix          get_reduced();

        // Returns a pointer to b where Ax = b
        Vect*           get_matrix_vector_product(Vect*);

        // Console representation of coefficient matrix for testing
        void            print_matrix();

        // Returns an nxn identity matrix
        static Matrix   get_identity_matrix_of_size(int);

    private:

        // The number of rows in A
        unsigned int    m;

        // The number of columns in A
        unsigned int    n;

        /*
        2D array using pointer notation representing a matrix
        */
        float**         A;
    };

    //========================================================================//
    // END OF MATRIX CLASS DECLARATION
    //========================================================================//

};