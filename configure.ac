AC_PREREQ([2.68])

AC_INIT([bmDCA],
        [0.8.1],
        [https://github.com/ranganathanlab/bmDCA/issues],
        [bmDCA])

AC_CANONICAL_HOST

build_unix=no
build_windows=no

# Detect the target system
case "${host_os}" in
  linux*|darwin*)
    build_unix=yes
    ;;
  cygwin*|mingw*|msys*)
    build_windows=yes
    ;;
  *)
    AC_MSG_ERROR(["OS $host_os is not supported"])
    ;;
esac

# Pass the conditionals to automake
AM_CONDITIONAL([UNIX], [test "$build_unix" = "yes"])
AM_CONDITIONAL([WINDOWS], [test "$build_windows" = "yes"])

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
