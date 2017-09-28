#                                               -*- Autoconf -*-

AC_INIT([npstat], [1.0.0], [lucaferretti@users-no-reply.github.com])

AM_INIT_AUTOMAKE([-Wall foreign])

AC_CONFIG_MACRO_DIR([aclocal])

AC_PROG_CC
AM_PROG_CC_C_O
GSL_MIN_VERSION=1.8
AX_PATH_GSL([$GSL_MIN_VERSION], [],
  [AC_MSG_ERROR([(Could not find GSL version >=$GSL_MIN_VERSION.])])

AC_CHECK_HEADER([stdio.h])

## From: https://www.gnu.org/software/gsl/manual/html_node/Autoconf-Macros.html#Autoconf-Macros
#AC_CHECK_LIB([m],[cos])
#AC_CHECK_LIB([gslcblas],[cblas_dgemm], [], [AC_MSG_ERROR(["libgslcblas is required"])], [])
#AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [], [AC_MSG_ERROR(["libgsl is required"])])

## autoscan output
AC_CHECK_FUNCS([getdelim])
AC_CHECK_FUNCS([sqrt])
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_PREREQ
AC_TYPE_SIZE_T


AC_CONFIG_FILES([Makefile])

AC_OUTPUT
