/* Abu Talha
   SMU Mathematics
   Math 6370 */

// Inclusions
#include "vec2d.hpp"


// This file implements the operations defined in the vec2d class


///// General class routines /////

// constructor (initializes values to 0.0)
vec2d::vec2d(long int m, long int n) {
  // if m or n is illegal, create an empty vector
  if (m < 1 || n < 1) {
    this->rows = 0;
    this->cols = 0;
    this->length = 0;
    this->data = NULL;
  } else {
    this->rows = m;
    this->cols = n;
    this->length = m*n;
    this->data = new double*[m];
    for (long int i=0; i<m; i++)
      this->data[i] = new double[n]();
  }
}


// destructor (frees space associated with a vec2d)
vec2d::~vec2d() {
  if (this->data!=NULL)
    for (long int i=0; i<this->rows; i++)
      delete[] this->data[i];
    delete[] this->data;
    this->rows = 0;
    this->cols = 0;
    this->length = 0;
}


// write myself to stdout
void vec2d::Write() const {
  // throw exception if data array isn't allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec1d::Write error, empty data array" );

  // print data to screen
  for (long int i=0; i<this->rows; i++){
    for (long int j=0; j<this->cols; j++)
        std::cout << std::setprecision(17) << "  " << this->data[i][j] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
}


// write myself to a file
void vec2d::Write(const char *outfile) const {
  // throw exception if data array isn't allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d::Write error, empty data array" );

  // throw exception if 'outfile' is empty
  if (strlen(outfile) < 1)
    throw std::invalid_argument( "vec2d::Write error, empty outfile" );

  // open output file
  std::ofstream fstr;
  fstr.open(outfile);
  if (!fstr.is_open())
    throw std::invalid_argument( "vec2d::Write error, unable to open file for writing" );

  // print data to file
  for (long int i=0; i<this->rows; i++){
    for (long int j=0; j<this->cols; j++)
        fstr << std::setprecision(17) << "  " << this->data[i][j];
      fstr << std::endl;
  }
  fstr << std::endl;

  // close output file
  fstr.close();
}

///// Arithmetic operations defined on a given vec2d /////

// x = a*y + b*z
void vec2d::LinearSum(double a, const vec2d &y, double b, const vec2d &z) {
  // check that array sizes match
  if (y.rows != this->rows || z.rows != this ->rows || y.cols != this->cols || z.cols != this->cols)
    throw std::invalid_argument( "vec2d::LinearSum error, vector sizes do not match" );

  // check that data is not NULL
  if (this->data == NULL || y.data == NULL || z.data == NULL)
    throw std::invalid_argument( "vec2d::LinearSum error: empty data array" );

  // perform operation
  for (long int i=0; i<this->rows; i++)
    for (long int j=0; j<this->cols; j++)
        this->data[i][j] = a*y.data[i][j] + b*z.data[i][j];
}


//   x = x*a  (scales my data by scalar a)
void vec2d::Scale(double a) {
  // check that data is not NULL
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d::Scale error: empty data array" );

  // perform operation
  for (long int i=0; i<this->rows; i++)
    for (long int j=0; j<this->cols; j++)
        this->data[i][j] *= a;
}


//   x = y  (copies y into x)
void vec2d::Copy(const vec2d& y) {
  // check that array sizes match
  if (y.rows != this->rows || y.cols != this->cols)
    throw std::invalid_argument( "vec2d::Copy error, vector sizes do not match" );

  // check that data is not NULL
  if (this->data == NULL || y.data == NULL)
    throw std::invalid_argument( "vec2d::Copy error: empty data array" );

  // perform operation
  for (long int i=0; i<this->rows; i++)
    for (long int j=0; j<this->cols; j++)
      this->data[i][j] = y.data[i][j];
}


//   x = a  (sets all entries of x to the scalar a)
void vec2d::Constant(double a) {
  // check that data is not NULL
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d::Copy error: empty data array" );

  // perform operation and return
  for (long int i=0; i<this->rows; i++){
    for (long int j=0; j<this->cols; j++)
        this->data[i][j] = a;
  }
}


///// scalar quantities derived from vectors /////

// min x_i
double vec2d::Min() {
  // check that my data is allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d::Min error: empty data array" );

  // perform operation and return
  double mn = this->data[0][0];
  for (long int i=0; i<this->rows; i++)
    for (long int j=0; j<this->cols; j++)
        mn = std::min(mn, this->data[i][j]);
  return mn;
}


// max x_i
double vec2d::Max() {
  // check that my data is allocated
  if (this->data == NULL)
    throw std::invalid_argument( "vec2d::Max error: empty data array" );

  // perform operation and return
  double mx = this->data[0][0];
  for (long int i=0; i<this->rows; i++)
    for (long int j=0; j<this->cols; j++)
        mx = std::max(mx, this->data[i][j]);
  return mx;
}

///// independent constructor routines /////

// create a vector of linearly spaced data
vec2d Linspace(double a, double b, long int m, long int n) {
  vec2d *x = new vec2d(m,n);
  double **xd = x->GetData();
  double h = (b-a)/(m*n-1);
  for (long int i=0; i<m; i++)
    for (long int j=0; j<n; j++)
        xd[i][j] = a + h*(i*n + j);
  return *x;
}


// create a vector of uniformly-distributed random data
vec2d Random(long int m, long int n) {
  vec2d *x = new vec2d(m,n);
  double **xd = x->GetData();
  for (long int i=0; i<m; i++)
    for (long int j=0; j<n; j++)
        xd[i][j] = random() / (pow(2.0,31.0) - 1.0);
  return *x;
}


///// independent arithmetic routines /////

// dot-product of x and y
double Dot(const vec2d &x, const vec2d &y) {
  // check that array sizes match
    if (y.Rows() != x.Rows() || y.Cols() != x.Cols())
        throw std::invalid_argument( "vec2d::Dot error, vector sizes do not match" );
  // perform operation and return
  double sum = 0.0;

  for (long int i=0; i<x.Rows(); i++)
      for (long int j=0; j<x.Cols(); j++)
          sum += x[i][j]*y[i][j];
  return sum;
}



// ||x||_2
double TwoNorm(const vec2d &x) {
  double sum = 0.0;
  for (long int i=0; i<x.Rows(); i++)
    for (long int j=0; j<x.Cols(); j++)
        sum += x[i][j]*x[i][j];
  return sqrt(sum);
}


// ||x||_RMS
double RmsNorm(const vec2d &x) {
  double sum = 0.0;
  for (long int i=0; i<x.Rows(); i++)
    for (long int j=0; j<x.Cols(); j++)
        sum += x[i][j]*x[i][j];
    return sqrt(sum/x.Length());
}


// ||x||_infty
double MaxNorm(const vec2d &x) {
  double mx = 0.0;
  for (long int i=0; i<x.Rows(); i++)
    for (long int j=0; j<x.Cols(); j++)
        mx = std::max(mx, std::abs(x[i][j]));
  return mx;
}
