#
#				configure.ac
#
# Process this file with autoconf to produce a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SWarp
#
#	Copyright:		(C) 2000-2019 IAP/CNRS/SorbonneU
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
#	Last modified:		20/12/2019
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT(SWarp, 2.41.0, [astromatic@astromatic.net],
	swarp, [http://astromatic.net/software/swarp])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date1=`date -R`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "AstrOmatic")
AC_SUBST(DATE1, "$date1")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_cfitsio.m4)
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_pthread.m4)
sinclude(acx_urbi_resolve_dir.m4)

# Provide a special option for setting the package release number
AC_ARG_WITH(release,
	[AS_HELP_STRING([--with-release=<release number>],
	[set the package release number (default = 1)])],
	[],
	[with_release="no"])
if test "$with_release" = "no"; then
  with_release=1
fi
AC_SUBST(PACKAGE_RELEASE, "$with_release")

# Display pakage and version number
AC_MSG_RESULT([******** Configuring:  $PACKAGE_NAME $PACKAGE_VERSION - $PACKAGE_RELEASE ($date) ********])

# Provide special option for choosing automatically the compilation flags.
AC_MSG_CHECKING([if compilation flags are set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([whether the INTEL compiler is enabled])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Use the Intel compiler (default = no)])],
        CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Checks for programs.
AC_LANG(C)

AC_SEARCH_LIBS([strerror],[cposix])
if test "$enable_auto_flags" = "yes"; then
  CFLAGS=""
  LDFLAGS=""
  ACX_PROG_CC_OPTIM
fi
LT_INIT
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h unistd.h sys/mman.h)
if test "$enable_icc" = "yes"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv gethostname memcpy memmove memset mmap strstr \
		getrlimit])
AC_CHECK_FUNCS([cosd sind tand acosd asind atand atan2d sincos])
AC_CHECK_FUNC([isnan], AC_DEFINE_UNQUOTED([HAVE_ISNAN2], 1,
		[Second isnan check]))

# Check support for large files
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Provide special options for CFITSIO
AC_ARG_WITH(cfitsio,
        [Link with the CFITSIO library])
AC_ARG_WITH(cfitsio-libdir,
        [AS_HELP_STRING([--with-cfitsio-libdir=<CFITSIO library path>],
        [Provide an alternative path to the CFITSIO library])])
AC_ARG_WITH(cfitsio-incdir,
        [AS_HELP_STRING([--with-cfitsio-incdir=<CFITSIO include dir>],
        [Provide an alternative path to the CFITSIO include directory])])

# Provide a special option for the default XSLT URL
with_xsl_url="file://"$(URBI_RESOLVE_DIR([$datadir]))"/$PACKAGE_TARNAME/$PACKAGE_TARNAME.xsl"
AC_ARG_WITH(xsl_url,
	[AS_HELP_STRING([--with-xsl_url=<default URL for XSLT filter>],
	[Provide an alternative default URL of the XSLT filter])])

AC_DEFINE_UNQUOTED([XSL_URL], "$with_xsl_url",[Default URL of the XSLT filter])

# Set flags for multithreading
n_pthreads=1024
AC_ARG_ENABLE(threads,
	[AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 1024 threads by default)])],
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

# Provide special option for profiling
AC_MSG_CHECKING([for profiler mode])
AC_ARG_ENABLE(profiling,
	[AS_HELP_STRING([--enable-profiling],
	[Enable special mode for profiling (default = no)])]
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (default = no)])],
	AC_MSG_RESULT([yes]),
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
  LIBS="$LIBS $PTHREAD_LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

########################## handle the CFITSIO library ########################
if test "$with_cfitsio" != "no"; then
  ACX_CFITSIO($with_cfitsio_libdir, $with_cfitsio_incdir,
	[LIBS="$CFITSIO_LIBS $LIBS"]
	if test "$CFITSIO_WARN" != ""; then
	  AC_MSG_WARN([$CFITSIO_WARN])
	fi,
	AC_MSG_ERROR([$CFITSIO_ERROR Exiting.]))
fi

# Compile with profiling option
if test "$enable_profiling" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  enable_best_link="no"
fi

# "Best" linking option
if test "$enable_best_link" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_LDFLAGS="-static-intel -qopenmp-link static -shared-libgcc \
	-static-libtool-libs -no-intel-extensions -avoid-version $AM_LDFLAGS"
  else
    AM_LDFLAGS="-shared-libgcc -static-libtool-libs -avoid-version $AM_LDFLAGS"
  fi
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
swarp.spec
debian/Makefile
debian/changelog
doc/Makefile
doc/src/conf.py
man/Makefile
man/swarp.1
src/Makefile
src/fits/Makefile
src/wcs/Makefile
])
AC_OUTPUT
