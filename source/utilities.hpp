/*! \file utilities.hpp
  \author Almog Yalinewich
  \brief A collection of general useful functions
 */

#include <vector>
#include <cmath>

using std::vector;
using std::size_t;

/*! \brief Minimum term in a vector
  \param list A vector
  \return Minimum term
 */
template<class T> T min_term(const vector<T>& list)
{
  T res = list[0];
  for(size_t i=1, endp=list.size();i<endp;++i)
    res = fmin(res,list[i]);
  return res;
}

//! \brief Lazy list
template<class T> class Index2Member
{
public:

  /*! \brief returns the length of the list
    \return Length of the list
   */
  virtual size_t getLength(void) const = 0;

  /*! \brief Calculates list member
    \param i Index
    \return i'th list member
   */
  virtual T operator()(size_t i) const = 0;

  //! \brief Class destructor
  virtual ~Index2Member(void) {}
};

/*! \brief Evaluates all members of a lazy list
  \param i2m Lazy list
  \return A vector
 */
template<class T> vector<T> serial_generate(const Index2Member<T>& i2m)
{
  vector<T> res(i2m.getLength());
  const size_t endp = res.size();
  //  #pragma omp parallel for
  for(size_t i=0;i<endp;++i)
    res[i] = i2m(i);
  return res;
}

/*! \brief Calculates the difference between consecutive vector members
  \param arg Vector
  \return Differences
 */
template<class T> vector<T> diff(const vector<T>& arg)
{
  class Differ: public Index2Member<T>
  {
  public:

    Differ(const vector<T>& arg_i):
      arg_(arg_i) {}

    size_t getLength(void) const
    {
      return arg_.size() - 1;
    }

    T operator()(size_t i) const
    {
      return arg_[i+1] - arg_[i];
    }

  private:
    const vector<T>& arg_;
  };

  return serial_generate(Differ(arg));
}

/*! \brief Increases a list by specified amounts
  \param change List of changes
  \param target Array to be changed
 */
template<class T> void serial_increment(const vector<T>& change,
				   vector<T>& target)
{
  for(size_t i=0, endp=target.size();i<endp;++i)
    target[i] += change[i];
}

/*! \brief Decreases a list by specified amounts
  \param change List of changes
  \param target Array to be changed
 */
template<class T> void serial_decrement(const vector<T>& change,
				   vector<T>& target)
{
  for(size_t i=0, endp=target.size();i<endp;++i)
    target[i] -= change[i];
}

/*! \brief Termwise product
  \param v1 Left argument
  \param v2 Right argument
  \return List of products
 */
template<class T> vector<T> serial_product(const vector<T>& v1,
					   const vector<T>& v2)
{
  assert(v1.size()==v2.size() && "Vectors must have same length");
  class Multiplier: public Index2Member<T>
  {
  public:

    Multiplier(const vector<T>& v1_i,
	       const vector<T>& v2_i):
      v1_(v1_i), v2_(v2_i) {}

    size_t getLength(void) const
    {
      return v1_.size();
    }

    T operator()(size_t i) const
    {
      return v1_[i]*v2_[i];
    }

  private:
    const vector<T>& v1_;
    const vector<T>& v2_;
  };

  return serial_generate(Multiplier(v1,v2));
}

/*! \brief Multiplies a list by a scalar
  \param v List
  \param s Scalar
  \return List multiplied by a scalar
 */
template<class T> vector<T> serial_product(const vector<T>& v,
					   double s)
{
  class Multiplier: public Index2Member<T>
  {
  public:

    Multiplier(const vector<T>& v_i, double s_i):
      v_(v_i), s_(s_i) {}

    size_t getLength(void) const
    {
      return v_.size();
    }

    double operator()(size_t i) const
    {
      return v_[i]*s_;
    }

  private:
    const vector<T>& v_;
    const double s_;
  } multiplier(v,s);

  return serial_generate(multiplier);
}
