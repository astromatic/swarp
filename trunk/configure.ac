#
#				configure.ac
#
# Process this file with autoconf to produce a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SWarp
#
#	Copyright:		(C) 2000-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
#
#	SWarp is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SWarp is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SWarp.  If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		17/12/2010
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT(swarp, 2.21, [bertin@iap.fr])
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_AUX_DIR(autoconf)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "Emmanuel Bertin")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_pthread.m4)
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_urbi_resolve_dir.m4)

# Display pakage and version number
AC_MSG_RESULT([*********** Configuring:  $PACKAGE $VERSION ($date) **********])

# Initialize the list of compilers to consider
cclist="cc gcc"

# Provide special option for choosing automatically the compilation flags.
AC_MSG_CHECKING([if compilation flags are set automatically])
AC_ARG_ENABLE(auto-flags,
	[AC_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags \
(off by default)])],
	use_autoflag="yes"
	AC_MSG_RESULT([yes]),
	use_autoflag="no"
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([for Linux Intel C compiler mode])
AC_ARG_ENABLE(icc,
	[AC_HELP_STRING([--enable-icc],
	[Enable special mode for compilation with the Intel compiler \
(off by default)])],
        use_icc="yes"
        cclist="icc $cclist"
	AC_MSG_RESULT([yes]),
        use_icc="no"
	AC_MSG_RESULT([no]))

# Checks for programs.
# GCC is chosen last because it is likely to yield less optimized code
AC_PROG_CC([$cclist])
# C Compiler: Check that it is ANSI C
AM_PROG_CC_STDC
# C Compiler: Check that it is POSIX-compliant
AC_ISC_POSIX
# Override automatic CFLAGS and LDFLAGS with those of user
if test "$use_autoflag" = "yes"; then
  CFLAGS=""
  LDFLAGS =""
  ACX_PROG_CC_OPTIM
fi
AC_PROG_RANLIB
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h unistd.h sys/mman.h)
if test $use_icc = "yes" -a $CC = "icc"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_TYPE_SIGNAL
AC_TYPE_UID_T
AC_CHECK_TYPES([unsigned long long, long long])
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MMAP
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv gethostname memcpy memmove memset mmap strstr])
AC_CHECK_FUNCS([cosd sind tand acosd asind atand atan2d sincos])

# Check support for large files
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Set the data directory to a true absolute path
datadir2=$(URBI_RESOLVE_DIR([$datadir]))

# Provide a special option for the default XSLT URL
AC_ARG_WITH(xsl_url,
        [AC_HELP_STRING([--with-xsl_url=<default URL for XSLT filter>],
        [Provide an alternative default URL of the XSLT filter])],
    [xsl_url=$withval],
    [xsl_url="file://$datadir2/$PACKAGE_NAME/$PACKAGE_NAME.xsl"]
    )

AC_DEFINE_UNQUOTED([XSL_URL], "$xsl_url",[Default URL of the XSLT filter])

# Set flags for multithreading
n_pthreads=16
AC_ARG_ENABLE(threads,
        [AC_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 16 threads by default)])],
    if test "$enableval" = "no"; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test "$enableval" != "yes"; then
        n_pthreads=$enableval
      fi
    fi,
    use_pthreads=yes
    )

# Provide special option for gprof profiling
AC_MSG_CHECKING([for gprof profiler mode])
AC_ARG_ENABLE(gprof,
	[AC_HELP_STRING([--enable-gprof],
	[Enable special mode for compilation with the gprof profiler \
(off by default)])],
        use_gprof="yes"
	AC_MSG_RESULT([yes]),
        use_gprof="no"
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AC_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (off by default)])],
	use_best="yes"
	AC_MSG_RESULT([yes]),
	use_best="no"
	AC_MSG_RESULT([no]))

################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$PTHREAD_LIBS $LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

# Link with gprof option
if test "$use_gprof" = "yes"; then
  if test "$use_icc" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  use_best="no"
fi

# "Best" linking option
if test "$use_best" = "yes"; then
  AM_LDFLAGS="-shared-libgcc -static-libtool-libs $AM_LDFLAGS"
fi

AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $AM_CPPFLAGS $CPPFLAGS $AM_CFLAGS $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $AM_LDFLAGS $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
man/Makefile
src/Makefile
src/fits/Makefile
src/wcs/Makefile
swarp.spec
man/swarp.1])
AC_OUTPUT
