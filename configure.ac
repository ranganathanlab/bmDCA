AC_PREREQ([2.68])

AC_INIT([bmDCA],
        [0.6.2],
        [https://github.com/sudorook/bmDCA/issues],
        [bmDCA])

AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([enable])
AC_PREFIX_DEFAULT(/usr/local)
AM_SILENT_RULES([yes])

AC_PROG_CXX
AC_LANG([C++])

# Check for required libraries and headers.
PKG_CHECK_MODULES([ARMADILLO], [armadillo])
AC_SEARCH_LIBS([cos], [m], [], [AC_MSG_ERROR([-lm missing])])
AC_OPENMP
AC_CHECK_HEADERS([stdio.h assert.h], [], [AC_MSG_ERROR("missing headers")], [])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_UINT16_T

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
