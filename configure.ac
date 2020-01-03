AC_PREREQ([2.68])

AC_INIT([bmDCA],
        [0.3.0],
        [https://github.com/sudorook/bmDCA/issues],
        [bmDCA])

AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([enable])
AC_PREFIX_DEFAULT(/usr/local)
AM_SILENT_RULES([yes])

# AC_PROG_CC_C99
# AC_LANG([C])
AC_PROG_CXX
AC_LANG([C++])

AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([-lm is missing, so your life has no meaning.])
])

# Checks for header files.
AC_CHECK_HEADERS([locale.h stdint.h stdlib.h string.h unistd.h signal.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_UINT16_T

AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
