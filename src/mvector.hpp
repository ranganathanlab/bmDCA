/* mvector.hpp
 *
 * Copyright (C) 2011 Carlo Baldassi (the "Author") <carlobaldassi@gmail.com>.
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the Licence, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org.licences/>.
 */

#ifndef MVECTOR_HPP
#define MVECTOR_HPP

#include <cstdarg>
#include <ostream>
#include <vector>

#ifdef XSTD_DEBUG
#include <iostream>
#define XSTD_DBGOUT(x) (std::cerr << x << std::endl)
#else
#define XSTD_DBGOUT(x)
#endif

namespace xstd {

template<size_t d>
class mshape
{
public:
  explicit mshape<d>(size_t n0, ...)
    : data(new size_t[d])
    , refs(new size_t[1])
  {
    XSTD_DBGOUT("constructor data=" << data << " refs=" << refs);
    data[0] = n0;
    va_list ap;
    va_start(ap, n0);
    for (size_t i = 1; i < d; ++i) {
      data[i] = va_arg(ap, size_t);
    }
    va_end(ap);
    *refs = 1;
  }

  ~mshape<d>()
  {
    XSTD_DBGOUT("destructor data=" << data << " refs=" << refs);
    if (refs) {
      XSTD_DBGOUT("  *refs=" << *refs);
      if (--(*refs) == 0) {
        XSTD_DBGOUT("  DELETE");
        delete[] data;
        delete[] refs;
      }
    }
  }

  mshape<d>(mshape<d> const& other)
    : data(other.data)
    , refs(other.refs)
  {
    XSTD_DBGOUT("copy data=" << data << " refs=" << refs);
    if (refs) {
      ++(*refs);
      XSTD_DBGOUT("  *refs=" << *refs);
    }
  }

  mshape<d>& operator=(mshape<d> const& other)
  {
    XSTD_DBGOUT("assign new data=" << other.data << " refs=" << other.refs);
    XSTD_DBGOUT("       old data=" << data << " refs=" << refs);
    this->~mshape<d>();
    data = other.data;
    refs = other.refs;
    if (refs) {
      ++(*refs);
      XSTD_DBGOUT("  *refs=" << *refs);
    }
    return *this;
  }

private:
  mshape<d>(size_t* o_d)
    : data(o_d)
    , refs(0)
  {
    XSTD_DBGOUT("lazy copy data=" << data);
  }

public:
  size_t current() const { return *data; }

  mshape<d - 1> next() const { return mshape<d - 1>(data + 1); }

  friend class mshape<d + 1>;

private:
  size_t* data;
  size_t* refs;
};

template<>
class mshape<0>
{
  size_t current() const { return 0; }
  void next() const {}
};

template<size_t d, typename T>
class mvector : public std::vector<mvector<d - 1, T>>
{
public:
  mvector<d, T>(size_t n = 0, T const& t = T())
    : std::vector<mvector<d - 1, T>>(n, mvector<d - 1, T>(n, t))
  {}
  mvector<d, T>(size_t v[], T const& t = T())
    : std::vector<mvector<d - 1, T>>(v[0], mvector<d - 1, T>(v + 1, t))
  {}
  mvector<d, T>(mshape<d> v, T const& t = T())
    : std::vector<mvector<d - 1, T>>(v.current(),
                                     mvector<d - 1, T>(v.next(), t))
  {}

  void reshape(size_t n, T const& t = T())
  {
    XSTD_DBGOUT("reshape d=" << d << " n=" << n);
    this->resize(n, mvector<d - 1, T>(n, t));
    for (typename mvector<d, T>::iterator i = this->begin(); i != this->end();
         ++i) {
      i->reshape(n, t);
    }
  }

  void reshape(size_t v[], T const& t = T())
  {
    XSTD_DBGOUT("reshape d=" << d << " v=" << v[0]);
    this->resize(v[0], mvector<d - 1, T>(v + 1, t));
    for (typename mvector<d, T>::iterator i = this->begin(); i != this->end();
         ++i) {
      i->reshape(v + 1, t);
    }
  }

  void reshape(mshape<d> v, T const& t = T())
  {
    XSTD_DBGOUT("reshape d=" << d << " v=" << v.current());
    this->resize(v.current(), mvector<d - 1, T>(v.next(), t));
    for (typename mvector<d, T>::iterator i = this->begin(); i != this->end();
         ++i) {
      i->reshape(v.next(), t);
    }
  }
};

template<typename T>
class mvector<1, T> : public std::vector<T>
{
public:
  mvector<1, T>(size_t n = 0, T const& t = T())
    : std::vector<T>(n, t)
  {}
  mvector<1, T>(size_t v[], T const& t = T())
    : std::vector<T>(v[0], t)
  {}
  mvector<1, T>(mshape<1> const& v, T const& t = T())
    : std::vector<T>(v.current(), t)
  {}
  void reshape(size_t n, T const& t = T()) { this->resize(n, t); }

  void reshape(size_t v[], T const& t = T()) { this->resize(v[0], t); }

  void reshape(mshape<1> v, T const& t = T()) { this->resize(v.current(), t); }
};

// This is defined mainly to avoid infinte recursion;
// it should not be used.
template<typename T>
class mvector<0, T> : public T
{
public:
  mvector<0, T>(T const& t = T())
    : T(t)
  {}
  mvector<0, T>(mshape<0> v, T const& t = T())
    : T(t)
  {}
  void reshape(mshape<0> v, T const& t = T()) {}
};

} // namespace xtl end

template<size_t d, typename T>
std::ostream&
operator<<(std::ostream& os, xstd::mvector<d, T> const& mv)
{
  os << "[";
  for (size_t i = 0, n = mv.size(); i < n; ++i) {
    os << mv[i];
    if (i < n - 1) {
      // os << std::endl;
      os << ",";
    }
  }
  os << "]";
  return os;
}

template<typename T>
std::ostream&
operator<<(std::ostream& os, xstd::mvector<0, T> const& mv)
{
  os << dynamic_cast<T const&>(mv);
  return os;
}

#endif // MVECTOR_HPP
