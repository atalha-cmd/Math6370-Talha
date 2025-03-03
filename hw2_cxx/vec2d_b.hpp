/* Abu Talha
   SMU Mathematics
   Math 6370 */

   #ifndef VEC2D_B_DEFINED__
   #define VEC2D_B_DEFINED__
   //    #define idx(i,j,cols) ((i*this->cols)+j)
   
   // Inclusions
   #include <stdlib.h>
   #include <iostream>
   #include <fstream>
   #include <iomanip>
   #include <string.h>
   #include <algorithm>
   #include <numeric>
   #include <cmath>
   #include <stdexcept>
   
   
   // This defines a simple arithmetic vector class
   class vec2d_b {
   
    private:
   
     ///// Contents /////
     long int rows;
     long int cols;
     long int length; // total number of elements
     double *data;
   
    public:
   
     ///// General class routines /////
     // constructor (initializes values to 0.0)
     vec2d_b(long int m, long int n);
   
     // destructor
     ~vec2d_b();
   
     // write myself to stdout
     void Write() const;
   
     // write myself to a file
     void Write(const char* outfile) const;
   
     // returns my number of rows, columns, and overall length
     long int Rows() const {return rows;};
     long int Cols() const {return cols;};
     long int Length() const {return length;};
   
     // access my data array
     double* GetData() const {return data;};
   
     // Overloaded () operator for (row, col) access
      double& operator() (long int i, long int j ) {return data[i*cols + j];};
      double operator() (long int i, long int j) const {return data[i*cols + j];};
      
      // Overload [] operator for direct 1D access
      double& operator[](long int index) { return data[index]; }
      double operator[](long int index) const { return data[index];}
   
     ///// Arithmetic operations defined on a vec1d /////
   
     // in-place operations (x is the calling vec1d) -- 0=success, 1=fail
     void LinearSum(double a, const vec2d_b &y,      // x = a*y + b*z
                    double b, const vec2d_b &z);
     void Scale(double a);                         // x = x*a
     void Copy(const vec2d_b &y);                    // x = y
     void Constant(double a);                      // x = a
   
     // scalar quantites derived from vectors
     double Min();                                 // min x_i
     double Max();                                 // max x_i
   
   };  // end vec2d_b
   
   
   // independent constructor routines
   vec2d_b Linspace(double a, double b, long int m, long int n);
   vec2d_b Random(long int m, long int n);
   
   // independent arithmetic routines
   double Dot(const vec2d_b &x, const vec2d_b &y);     // sum_i (x_i * y_i)
   double TwoNorm(const vec2d_b &x);                 // sqrt(sum_i x_i^2)
   double RmsNorm(const vec2d_b &x);                 // sqrt(sum_i x_i^2 / n)
   double MaxNorm(const vec2d_b &x);                 // max_i |x_i|
   
   #endif

