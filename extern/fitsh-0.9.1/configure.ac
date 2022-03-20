# Initialize autoconf:
AC_INIT(fitsh, 0.9.1, apal@szofi.net)

RELEASE_DATE=2011.11.11

###############################################################################

# ===========================================================================
#      http://www.nongnu.org/autoconf-archive/ax_cflags_gcc_option.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CFLAGS_GCC_OPTION (optionflag [,[shellvar][,[A][,[NA]]])
#
# DESCRIPTION
#
#   AX_CFLAGS_GCC_OPTION(-fvomit-frame) would show a message as like
#   "checking CFLAGS for gcc -fvomit-frame ... yes" and adds the optionflag
#   to CFLAGS if it is understood. You can override the shellvar-default of
#   CFLAGS of course. The order of arguments stems from the explicit macros
#   like AX_CFLAGS_WARN_ALL.
#
#   The cousin AX_CXXFLAGS_GCC_OPTION would check for an option to add to
#   CXXFLAGS - and it uses the autoconf setup for C++ instead of C (since it
#   is possible to use different compilers for C and C++).
#
#   The macro is a lot simpler than any special AX_CFLAGS_* macro (or
#   ac_cxx_rtti.m4 macro) but allows to check for arbitrary options.
#   However, if you use this macro in a few places, it would be great if you
#   would make up a new function-macro and submit it to the ac-archive.
#
#     - $1 option-to-check-for : required ("-option" as non-value)
#     - $2 shell-variable-to-add-to : CFLAGS (or CXXFLAGS in the other case)
#     - $3 action-if-found : add value to shellvariable
#     - $4 action-if-not-found : nothing
#
#   Note: in earlier versions, $1-$2 were swapped. We try to detect the
#   situation and accept a $2=~/-/ as being the old option-to-check-for.
#
#   There are other variants that emerged from the original macro variant
#   which did just test an option to be possibly added. However, some
#   compilers accept an option silently, or possibly for just another option
#   that was not intended. Therefore, we have to do a generic test for a
#   compiler family. For gcc we check "-pedantic" being accepted which is
#   also understood by compilers who just want to be compatible with gcc
#   even when not being made from gcc sources.
#
#   See also: AX_CFLAGS_SUN_OPTION, AX_CFLAGS_HPUX_OPTION,
#   AX_CFLAGS_AIX_OPTION, and AX_CFLAGS_IRIX_OPTION.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_CFLAGS_GCC_OPTION_OLD], [dnl
AS_VAR_PUSHDEF([FLAGS],[CFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cflags_gcc_option_$2])dnl
AC_CACHE_CHECK([m4_ifval($1,$1,FLAGS) for gcc m4_ifval($2,$2,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_SAVE
 AC_LANG_C
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic -Werror % m4_ifval($2,$2,-option)"  dnl   GCC
   "-pedantic % m4_ifval($2,$2,-option) %% no, obsolete"  dnl new GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_RESTORE
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($1,$1,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($1,$1,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"])
                      m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])


dnl the only difference - the LANG selection... and the default FLAGS

AC_DEFUN([AX_CXXFLAGS_GCC_OPTION_OLD], [dnl
AS_VAR_PUSHDEF([FLAGS],[CXXFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cxxflags_gcc_option_$2])dnl
AC_CACHE_CHECK([m4_ifval($1,$1,FLAGS) for gcc m4_ifval($2,$2,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic -Werror % m4_ifval($2,$2,-option)"  dnl   GCC
   "-pedantic % m4_ifval($2,$2,-option) %% no, obsolete"  dnl new GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_RESTORE
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($1,$1,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($1,$1,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"])
                      m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])

dnl -------------------------------------------------------------------------

AC_DEFUN([AX_CFLAGS_GCC_OPTION_NEW], [dnl
AS_VAR_PUSHDEF([FLAGS],[CFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cflags_gcc_option_$1])dnl
AC_CACHE_CHECK([m4_ifval($2,$2,FLAGS) for gcc m4_ifval($1,$1,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_SAVE
 AC_LANG_C
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic -Werror % m4_ifval($1,$1,-option)"  dnl   GCC
   "-pedantic % m4_ifval($1,$1,-option) %% no, obsolete"  dnl new GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_RESTORE
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($2,$2,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($2,$2,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"])
                      m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])


dnl the only difference - the LANG selection... and the default FLAGS

AC_DEFUN([AX_CXXFLAGS_GCC_OPTION_NEW], [dnl
AS_VAR_PUSHDEF([FLAGS],[CXXFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cxxflags_gcc_option_$1])dnl
AC_CACHE_CHECK([m4_ifval($2,$2,FLAGS) for gcc m4_ifval($1,$1,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic -Werror % m4_ifval($1,$1,-option)"  dnl   GCC
   "-pedantic % m4_ifval($1,$1,-option) %% no, obsolete"  dnl new GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_RESTORE
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($2,$2,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($2,$2,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"])
                      m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])

AC_DEFUN([AX_CFLAGS_GCC_OPTION],[ifelse(m4_bregexp([$2],[-]),-1,
[AX_CFLAGS_GCC_OPTION_NEW($@)],[AX_CFLAGS_GCC_OPTION_OLD($@)])])

AC_DEFUN([AX_CXXFLAGS_GCC_OPTION],[ifelse(m4_bregexp([$2],[-]),-1,
[AX_CXXFLAGS_GCC_OPTION_NEW($@)],[AX_CXXFLAGS_GCC_OPTION_OLD($@)])])

###############################################################################

# our default $(CFLAGS) are:
if test -n "$CFLAGS" ; then
	cflags_default=no
else
	CFLAGS="-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE"
	cflags_default=yes
fi

# some additional: 
CFLAGS="${CFLAGS} -D_FI_SOURCE"

dynamic_extensions=yes

AC_ARG_ENABLE(dynamic-extensions,
 [  --enable-dynamic-extensions		`lfit` may use dynamic extensions at run-time (default: yes)],
 [ if test "x$enableval" = "xno" ; then dynamic_extensions=no; else dynamic_extensions=yes; fi ])

# check compiler:
AC_PROG_CC

# check archiver:
AC_CHECK_PROG(ac_prog_ar,ar,ar,false)

# check linker:
AC_CHECK_PROG(ac_prog_ld,ld,ld,false)

# check installer:
AC_CHECK_PROG(ac_prog_install,install,install,false)

# Output of: cat `find . -name "*.[c]"` | grep "^#include <.*>$" | sort | uniq
#include <bzlib.h>
#include <ctype.h>
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <fnmatch.h>
#include <libgen.h>
#include <libintl.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

# help2man utility:
AC_MSG_CHECKING([help2man])
if help2man --version >/dev/null 2>&1 ; then
	AC_MSG_RESULT([help2man])
	ac_cv_help2man="yes"
	HELP2MAN=help2man
else
	AC_MSG_RESULT([no])
	ac_cv_help2man="no"
	HELP2MAN=true
fi

AC_MSG_CHECKING([for long string extension])
tmpdir=`mktemp -d`
echo	'char * somestring = __extension__ "QQRIQ";' > $tmpdir/test.c 
if $CC -c $tmpdir/test.c -o $tmpdir/test.o >/dev/null 2>&1 ; then
	AC_MSG_RESULT([yes])
	ac_cv_gcc_extension="yes"
else
	AC_MSG_RESULT([no])
	ac_cv_gcc_extension="no"
fi
rm -f -r $tmpdir

if test "x$ac_cv_gcc_extension" = "xno" ; then
	CFLAGS="${CFLAGS} -DHAVE_NO_CC_EXTENSION"
fi

# check used headers:
AC_CHECK_HEADERS(ctype.h)
AC_CHECK_HEADERS(dlfcn.h)
AC_CHECK_HEADERS(fcntl.h)
AC_CHECK_HEADERS(fnmatch.h)
AC_CHECK_HEADERS(math.h)
AC_CHECK_HEADERS(time.h)
AC_CHECK_HEADERS(errno.h)
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(signal.h)
AC_CHECK_HEADERS(sys/time.h)
AC_CHECK_HEADERS(sys/ioctl.h)
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(sys/stat.h)
AC_CHECK_HEADERS(sys/mman.h)
AC_CHECK_HEADERS(sys/select.h)
AC_CHECK_HEADERS(sys/socket.h)
AC_CHECK_HEADERS(sys/wait.h)

# Some OS-dependent functions:
AC_CHECK_LIB(c,fnmatch)

# dynamic library management/handling:
AC_CHECK_LIB(dl,dlopen)
AC_CHECK_LIB(dl,dlsym)
AC_CHECK_LIB(dl,dlclose)

ax_cflags_gcc_no_overlength_strings=no
AX_CFLAGS_GCC_OPTION(-Wno-overlength-strings,[],[ax_cflags_gcc_no_overlength_strings=yes])
if test $cflags_default == yes && test $ax_cflags_gcc_no_overlength_strings == yes ; then
	CFLAGS="$CFLAGS -Wno-overlength-strings"
fi

ax_cflags_gcc_no_long_long=no
AX_CFLAGS_GCC_OPTION(-Wno-long-long,[],[ax_cflags_gcc_no_long_long=yes])
if test $cflags_default == yes && test $ax_cflags_gcc_no_long_long == yes ; then
	CFLAGS="$CFLAGS -Wno-long-long"
fi

ax_cflags_gcc_no_unused_result=no
AX_CFLAGS_GCC_OPTION(-Wno-unused-result,[],[ax_cflags_gcc_no_unused_result=yes])
if test $cflags_default == yes && test $ax_cflags_gcc_no_unused_result == yes ; then
	CFLAGS="$CFLAGS -Wno-unused-result"
fi

# gzip:
AC_CHECK_HEADERS(zlib.h)
AC_CHECK_LIB(z,zlibVersion)
AC_CHECK_LIB(z,deflateInit2)
AC_CHECK_LIB(z,deflateInit2_)
AC_CHECK_LIB(z,deflate)
AC_CHECK_LIB(z,deflateEnd)
AC_CHECK_LIB(z,inflateInit2)
AC_CHECK_LIB(z,inflateInit2_)
AC_CHECK_LIB(z,inflate)
AC_CHECK_LIB(z,inflateEnd)
AC_CHECK_LIB(z,crc32)

# bzip2:
AC_CHECK_HEADERS(bzlib.h)
AC_CHECK_LIB(bz2,BZ2_bzlibVersion)
AC_CHECK_LIB(bz2,BZ2_bzBuffToBuffCompress)
AC_CHECK_LIB(bz2,BZ2_bzBuffToBuffDecompress)

FIZIP_CFLAGS=""
FIZIP_LIBS=""
fizip_supported="dummy"
if  test "x$ac_cv_header_zlib_h" == "xyes" && \
    test "x$ac_cv_lib_z_zlibVersion" == "xyes" && \
    test "x$ac_cv_lib_z_deflate" == "xyes" && \
    test "x$ac_cv_lib_z_inflate" == "xyes" ; then
	FIZIP_CFLAGS="${FIZIP_CFLAGS} -DENABLE_GZ"
	FIZIP_LIBS="${FIZIP_LIBS} -lz"
	fizip_supported="$fizip_supported gz"
fi
if  test "x$ac_cv_header_bzlib_h" == "xyes" && \
    test "x$ac_cv_lib_bz2_BZ2_bzlibVersion" == "xyes" && \
    test "x$ac_cv_lib_bz2_BZ2_bzBuffToBuffCompress" == "xyes" && \
    test "x$ac_cv_lib_bz2_BZ2_bzBuffToBuffDecompress" == "xyes" ; then
	FIZIP_CFLAGS="${FIZIP_CFLAGS} -DENABLE_BZ2"
	FIZIP_LIBS="${FIZIP_LIBS} -lbz2"
	fizip_supported="$fizip_supported bz2"
fi
	

# some additional functions:

# basic math functions:
AC_CHECK_LIB(m,sin)
AC_CHECK_LIB(m,cos)
AC_CHECK_LIB(m,tan)
AC_CHECK_LIB(m,asin)
AC_CHECK_LIB(m,acos)
AC_CHECK_LIB(m,atan)
AC_CHECK_LIB(m,atan2)
AC_CHECK_LIB(m,sqrt)
AC_CHECK_LIB(m,log)
AC_CHECK_LIB(m,exp)
AC_CHECK_LIB(m,log10)
AC_CHECK_LIB(m,exp10)
AC_CHECK_LIB(m,pow)
AC_CHECK_LIB(m,floor)
AC_CHECK_LIB(m,fmod)

AC_CHECK_LIB(m,finite)
AC_CHECK_LIB(m,isnan)
AC_CHECK_LIB(m,isinf)

# Debian archiver:
AC_MSG_CHECKING([for a debian archiver])
if dpkg-deb --version >/dev/null 2>&1 ; then
	AC_MSG_RESULT([dpkg-deb])
	ac_cv_debian_archiver="yes"
	DPKG_DEB=dpkg-deb
else
	AC_MSG_RESULT([no])
	ac_cv_debian_archiver="no"
	DPKG_DEB=false
fi

# Groff:
AC_MSG_CHECKING([for a roff formatter])
if groff --version >/dev/null 2>&1 ; then
	AC_MSG_RESULT([groff])
	ac_cv_roff_groff="yes"
	ROFF=groff
else
	AC_MSG_RESULT([no])
	ac_cv_roff_groff="no"
	ROFF=false
fi

AC_CHECK_PROG(ac_prog_man,man,man,false)
AC_CHECK_PROG(ac_prog_dvips,dvips,dvips,false)
AC_CHECK_PROG(ac_prog_gzip,gzip,gzip,false)

# Check paper size:
AC_MSG_CHECKING([papersize])
if test -x /usr/bin/paperconf ; then 
	PAPERSIZE=`/usr/bin/paperconf`
	PAPER_WIDTH=`/usr/bin/paperconf -w`
	PAPER_HEIGHT=`/usr/bin/paperconf -h`
elif test -f /etc/papersize ; then 
	PAPERSIZE=`grep -v "^\#" /etc/papersize | head -1`
	PAPER_WIDTH=0
	PAPER_HEIGHT=0
else
	PAPERSIZE=letter
	PAPER_WIDTH=612 
	PAPER_HEIGHT=792
fi
AC_MSG_RESULT([$PAPERSIZE: ${PAPER_WIDTH}x${PAPER_HEIGHT}])

AC_MSG_CHECKING([dynamic library extensions])
if test `uname` == Darwin ; then
	DLEXT=dylib
	DLSWC=-dylib
	DLDYN=-dynamic
	AC_MSG_RESULT([dylib])
else
	DLEXT=so
	DLSWC=-shared
	DLDYN=-rdynamic
	AC_MSG_RESULT([so])
fi

DEB_ARCH=`uname -m | sed -e 's/x86_64/amd64/' -e 's/i486/i386/' -e 's/i586/i386/' -e 's/i686/i386/'`
DEB_VERSION=$PACKAGE_VERSION

FI_VERSION=$PACKAGE_VERSION
FI_RELEASE=$RELEASE_DATE

# Final notices:
AC_MSG_NOTICE([generic: using CFLAGS: ${CFLAGS}])
AC_MSG_NOTICE([generic: using Debian architecture identifier: ${DEB_ARCH}])
AC_MSG_NOTICE([`lfit`: may use dynamic extensions: $dynamic_extensions])
AC_MSG_NOTICE([`fizip`: supported compression algorithms: $fizip_supported])

# Substitution for binary object archiver and linker:
AC_SUBST(AR, ${ac_prog_ar})
AC_SUBST(LD, ${ac_prog_ld})
AC_SUBST(INSTALL, ${ac_prog_install})

# Substitution for dynamic loading and additional optional supported methods:
AC_SUBST(DLEXT, ${DLEXT})
AC_SUBST(DLSWC, ${DLSWC})
AC_SUBST(DLDYN, ${DLDYN})
AC_SUBST(FIZIP_CFLAGS, ${FIZIP_CFLAGS})
AC_SUBST(FIZIP_LIBS, ${FIZIP_LIBS})

	
# Substitution for documentation and distribution management:
AC_SUBST(DEB_VERSION, ${DEB_VERSION})
AC_SUBST(DEB_ARCH, ${DEB_ARCH})
AC_SUBST(DPKG_DEB, ${DPKG_DEB})
AC_SUBST(HELP2MAN, ${HELP2MAN})
AC_SUBST(ROFF, ${ROFF})
AC_SUBST(GZIP, ${ac_prog_gzip})
AC_SUBST(MAN, ${ac_prog_man})
AC_SUBST(DVIPS, ${ac_prog_dvips})
AC_SUBST(PAPERSIZE, ${PAPERSIZE})
AC_SUBST(PAPER_WIDTH, ${PAPER_WIDTH})
AC_SUBST(PAPER_HEIGHT, ${PAPER_HEIGHT})

# Substitution for the config.h header:
AC_SUBST(FI_VERSION,${FI_VERSION})
AC_SUBST(FI_RELEASE,${FI_RELEASE})

ac_output_list="Makefile 
		config.h
		libastro/Makefile 
		libfits/Makefile
		libpsn/Makefile
		src/Makefile 
		src/index/Makefile 
		src/io/Makefile 
		src/link/Makefile 
		src/math/Makefile 
		src/math/dft/Makefile 
		src/math/elliptic/Makefile 
		src/math/elliptic/Makefile 
		src/math/expint/Makefile 
		src/math/fit/Makefile 
		src/math/intersec/Makefile 
		src/math/spline/Makefile 
		src/parallel/Makefile 
		src/psn/Makefile 
		doc/Makefile
		doc/examples/Makefile
		doc/examples/firandom/Makefile"

AC_OUTPUT([${ac_output_list}])

