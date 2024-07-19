AC_INIT([tcenum], [0.1])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
dnl AC_CONFIG_FILES([bin/tcenum], [chmod +x bin/tcenum])
AC_REQUIRE_AUX_FILE([tap-driver.sh])
AC_PROG_CC
AC_PROG_F77([gfortran egfortran])

AC_ARG_WITH([maxmib], AS_HELP_STRING([--with-maxmib=X],
	              [allocate X megabytes to the TC. Default is X=6000.]))
AS_IF([test x$with_maxmib = x], [maxmib=6000], [maxmib=$with_maxmib])

AC_ARG_WITH([intbytes], AS_HELP_STRING([--with-intbytes=Y],
	              [use Y-byte integers. Default is Y=4, the only other possibility is 8.]))
AS_IF([test x$with_intbytes = x], [intbytes=4], [intbytes=$with_intbytes])
AS_IF([test $intbytes = 4 -o $intbytes = 8], [], [AC_MSG_ERROR([intbytes = $intbytes is not allowed])])

FFLAGS="$FFFLAGS -DMAXMIB=$maxmib -DINTBYTES=$intbytes"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
