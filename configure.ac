AC_PREREQ([2.68])

AC_INIT([bmDCA],
        [0.7.0],
        [https://github.com/ranganathanlab/bmDCA/issues],
        [bmDCA])

AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([enable])
AC_PREFIX_DEFAULT(/usr/local)
AM_SILENT_RULES([yes])

AC_LANG([C++])
AC_PROG_CXX

AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_cxx_compile_stdcxx.m4])
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])

PKG_PROG_PKG_CONFIG

# Check for required libraries and headers.
PKG_CHECK_MODULES(
  [ARMADILLO],
  [armadillo],
  [],
  [AC_CHECK_HEADERS(
    [armadillo],
    [],
    [AC_MSG_ERROR("missing armadillo")],
    [])
  ])
AC_OPENMP
AC_CHECK_HEADERS(
  [stdio.h assert.h],
  [],
  [AC_MSG_ERROR("missing headers")],
  [])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_UINT16_T

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
