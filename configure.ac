AC_PREREQ([2.68])

AC_INIT([bmDCA],
        [0.1.0],
        [https://github.com/sudorook/bmDCA/issues],
        [bmDCA])

AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE([enable])
AC_PREFIX_DEFAULT(/usr/local)
AM_SILENT_RULES([yes])

AC_PROG_CC_C99
AC_LANG([C])
AC_PROG_CXX

# PKG_CHECK_MODULES([GLIB], [glib-2.0 gobject-2.0], [have_glib=yes], [have_glib=no])
AC_CHECK_PROG(BC_CHECK, bc, yes)
if test x"$BC_CHECK" != x"yes"; then
  AC_MSG_ERROR(["BC" executable is not installed.])
fi
AC_CHECK_PROG(AWK_CHECK, awk, yes)
if test x"$AWK_CHECK" != x"yes"; then
  AC_MSG_ERROR(["AWK" executable is not installed.])
fi

# Checks for header files.
AC_CHECK_HEADERS([locale.h stdint.h stdlib.h string.h unistd.h signal.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_UINT16_T

AC_CONFIG_FILES([Makefile scripts/Makefile src/Makefile src/mcmc/Makefile])

AC_OUTPUT
