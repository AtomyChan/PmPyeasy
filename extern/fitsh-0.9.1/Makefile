SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld
INSTALL=install

.PHONY: all install \
        libfits libpsn libastro src doc \
        info html pdf ps dvi \
        install install-strip \
        install-info install-html install-dvi install-pdf install-ps \
        uninstall \
        clean distclean mostlyclean maintainter-clean TAGS \
        deb deb-fitsh

DEB_VERSION = 0.9.1
DEB_ARCH = amd64
DPKG_DEB = false
GZIP	 = gzip

datarootdir = ${prefix}/share
prefix = /usr/local
exec_prefix = ${prefix}
mandir = ${datarootdir}/man
bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include

PACKAGE_NAME    = fitsh
PACKAGE_TARNAME = fitsh

infodir = ${datarootdir}/info
docdir  = ${datarootdir}/doc/${PACKAGE_TARNAME}
htmldir = ${docdir}
dvidir  = ${docdir}
pdfdir  = ${docdir}
psdir   = ${docdir}

DLEXT = dylib

HELP2MAN = true
ROFF	 = groff
MAN	 = man
DVIPS	 = dvips

ARCH_VERSION = $(DEB_VERSION)
DIST = fitsh-$(ARCH_VERSION)

all: src

libastro:
	$(MAKE) -C libastro

libfits: 
	$(MAKE) -C libfits

libpsn:
	$(MAKE) -C libpsn

src: libastro libfits libpsn
	$(MAKE) -C src

doc:
	$(MAKE) -C doc

clean:
	$(MAKE) -C libastro clean
	$(MAKE) -C libfits clean
	$(MAKE) -C libpsn clean
	$(MAKE) -C src clean
	$(MAKE) -C doc clean
	rm -f $(DIST).tar.gz

distclean: clean
	rm -f -r autom4te.cache
	rm -f ./config.status ./config.log
	rm -f ./Makefile
	rm -f ./config.h
	rm -f ./libastro/Makefile 
	rm -f ./libfits/Makefile
	rm -f ./libpsn/Makefile
	rm -f ./src/Makefile 
	rm -f ./src/index/Makefile 
	rm -f ./src/io/Makefile 
	rm -f ./src/link/Makefile 
	rm -f ./src/math/Makefile 
	rm -f ./src/math/dft/Makefile 
	rm -f ./src/math/elliptic/Makefile 
	rm -f ./src/math/elliptic/Makefile 
	rm -f ./src/math/expint/Makefile 
	rm -f ./src/math/fit/Makefile 
	rm -f ./src/math/intersec/Makefile 
	rm -f ./src/math/spline/Makefile 
	rm -f ./src/parallel/Makefile 
	rm -f ./src/psn/Makefile 
	rm -f ./doc/Makefile
	rm -f ./doc/examples/Makefile
	rm -f ./doc/examples/firandom/Makefile

install: src 
	$(INSTALL) -d ${bindir} 
	$(INSTALL) -m 0755 ./src/fiarith ${bindir}
	$(INSTALL) -m 0755 ./src/ficalib ${bindir}
	$(INSTALL) -m 0755 ./src/ficombine ${bindir}
	$(INSTALL) -m 0755 ./src/ficonv ${bindir}
	$(INSTALL) -m 0755 ./src/fiheader ${bindir}
	$(INSTALL) -m 0755 ./src/fiign ${bindir}
	$(INSTALL) -m 0755 ./src/fiinfo ${bindir}
	$(INSTALL) -m 0755 ./src/fiphot ${bindir}
	$(INSTALL) -m 0755 ./src/firandom ${bindir}
	$(INSTALL) -m 0755 ./src/fistar ${bindir}
	$(INSTALL) -m 0755 ./src/fitrans ${bindir}
	$(INSTALL) -m 0755 ./src/grcollect ${bindir}
	$(INSTALL) -m 0755 ./src/grmatch ${bindir}
	$(INSTALL) -m 0755 ./src/grtrans ${bindir}
	$(INSTALL) -m 0755 ./src/grselect ${bindir}
	$(INSTALL) -m 0755 ./src/grtfilter ${bindir}
	$(INSTALL) -m 0755 ./src/lfit ${bindir}
	$(INSTALL) -d ${libdir}/lfit
	$(INSTALL) -m 0644 ./src/linear.$(DLEXT) ${libdir}/lfit

uninstall:
	rm -f ${bindir}/fiarith
	rm -f ${bindir}/ficalib
	rm -f ${bindir}/ficombine
	rm -f ${bindir}/ficonv
	rm -f ${bindir}/fiheader
	rm -f ${bindir}/fiign
	rm -f ${bindir}/fiinfo
	rm -f ${bindir}/fiphot
	rm -f ${bindir}/firandom
	rm -f ${bindir}/fistar
	rm -f ${bindir}/fitrans
	rm -f ${bindir}/grcollect
	rm -f ${bindir}/grmatch
	rm -f ${bindir}/grtrans
	rm -f ${bindir}/grselect
	rm -f ${bindir}/grtfilter
	rm -f ${bindir}/lfit
	rm -f ${libdir}/lfit/linear.$(DLEXT)

HELP2MAN_LIST=fiarith ficalib ficombine ficonv fiheader fiign fiinfo fiphot firandom fistar fitrans grcollect grmatch grtrans lfit

help2man: src
	test -d ./doc/help2man || mkdir ./doc/help2man
	for b in $(HELP2MAN_LIST) ; do \
		$(HELP2MAN) -N -h --long-help ./src/$$b > ./doc/help2man/$$b.1 ; \
		$(ROFF) -Tdvi -mandoc ./doc/help2man/$$b.1 > ./doc/help2man/$$b.dvi 2>/dev/null ; \
		$(DVIPS) ./doc/help2man/$$b.dvi -q -o ./doc/help2man/$$b.ps ; \
		$(GZIP) -f ./doc/help2man/$$b.1 ; \
	done


dist:
	rm      -f -r $(DIST)
	mkdir   $(DIST)
	mkdir   $(DIST)/doc
	cp      -p doc/Makefile.in $(DIST)/doc
	mkdir	$(DIST)/doc/examples
	cp	-p doc/examples/Makefile.in $(DIST)/doc/examples
	mkdir	$(DIST)/doc/examples/firandom
	cp	-p doc/examples/firandom/Makefile.in $(DIST)/doc/examples/firandom
	mkdir   $(DIST)/include
	mkdir   $(DIST)/include/astro
	cp	-p include/astro/*.[ch] $(DIST)/include/astro
	mkdir   $(DIST)/include/fits
	cp	-p include/fits/*.[ch] $(DIST)/include/fits
	mkdir   $(DIST)/include/lfit
	cp	-p include/lfit/*.[ch] $(DIST)/include/lfit
	mkdir   $(DIST)/include/psn
	cp	-p include/psn/*.[ch] $(DIST)/include/psn
	mkdir   $(DIST)/libastro
	cp      -p libastro/*.[ch] libastro/Makefile.in $(DIST)/libastro
	mkdir   $(DIST)/libfits
	cp      -p libfits/*.[ch] libfits/Makefile.in $(DIST)/libfits
	mkdir   $(DIST)/libpsn
	cp      -p libpsn/*.[ch] libpsn/Makefile.in $(DIST)/libpsn
	mkdir   $(DIST)/src
	cp      -p src/*.[ch] src/Makefile.in $(DIST)/src
	mkdir   $(DIST)/src/index
	cp      -p src/index/*.[ch] src/index/Makefile.in $(DIST)/src/index
	mkdir   $(DIST)/src/io
	cp      -p src/io/*.[ch] src/io/Makefile.in $(DIST)/src/io
	mkdir   $(DIST)/src/link
	cp      -p src/link/*.[ch] src/link/Makefile.in $(DIST)/src/link
	mkdir   $(DIST)/src/math
	cp      -p src/math/*.[ch] src/math/Makefile.in $(DIST)/src/math
	mkdir   $(DIST)/src/math/dft
	cp      -p src/math/dft/*.[ch] src/math/dft/Makefile.in $(DIST)/src/math/dft
	mkdir   $(DIST)/src/math/elliptic
	cp      -p src/math/elliptic/*.[ch] src/math/elliptic/Makefile.in $(DIST)/src/math/elliptic
	mkdir   $(DIST)/src/math/expint
	cp      -p src/math/expint/*.[ch] src/math/expint/Makefile.in $(DIST)/src/math/expint
	mkdir   $(DIST)/src/math/fit
	cp      -p src/math/fit/*.[ch] src/math/fit/Makefile.in $(DIST)/src/math/fit
	mkdir   $(DIST)/src/math/intersec
	cp      -p src/math/intersec/*.[ch] src/math/intersec/Makefile.in $(DIST)/src/math/intersec
	mkdir   $(DIST)/src/math/spline
	cp      -p src/math/spline/*.[ch] src/math/spline/Makefile.in $(DIST)/src/math/spline
	mkdir   $(DIST)/src/parallel
	cp      -p src/parallel/*.[ch] src/parallel/Makefile.in $(DIST)/src/parallel
	mkdir   $(DIST)/src/psn
	cp      -p src/psn/*.[ch] src/psn/Makefile.in $(DIST)/src/psn
	cp      -p VERSION $(DIST)
	cp      -p ./configure configure.ac ./prepare $(DIST)
	cp      -p Makefile.in config.h.in $(DIST)
	tar     czf $(DIST).tar.gz $(DIST)
	rm      -f -r $(DIST)

deb: deb-fitsh

deb-fitsh: src help2man
	rm      -f -r deb/fitsh
	mkdir   -p deb/fitsh{,/DEBIAN}
	mkdir   -p deb/fitsh{,/usr/bin,/usr/share/man/man1,/usr/share/info,/usr/share/doc/fitsh}
	cat deb/fitsh.list | while read file target ; do \
		cp $$file deb/fitsh/$$target ; \
	done
	cp      README  deb/fitsh/usr/share/doc/fitsh/
	gzip    deb/fitsh/usr/share/doc/fitsh/README
	cp      AUTHORS  deb/fitsh/usr/share/doc/fitsh/
	gzip    deb/fitsh/usr/share/doc/fitsh/AUTHORS
	cat     deb/fitsh.control | \
	sed     -e "s/__ARCH__/$(DEB_ARCH)/" \
	        -e "s/__VERSION__/$(DEB_VERSION)/" \
	        -e "s/__INSTALLED_SIZE__/`du -k -s --exclude DEBIAN deb/fitsh | awk '{ print $$1; }'`/" > deb/fitsh/DEBIAN/control
	for f in `find deb/fitsh` ; do \
	        o=`echo $$f | sed 's/deb\/fitsh\///'` ; \
	        if [ -f "$$f" -a `dirname "$$f"` != deb/fitsh/DEBIAN ] ; then \
	                m=`md5sum $$f | awk '{ print $$1; }'` ; \
	                echo $$m $$o ; \
	        fi ; \
	done > deb/fitsh/DEBIAN/md5sums
	$(DPKG_DEB) --build deb/fitsh fitsh_$(DEB_VERSION)_$(DEB_ARCH).deb
	rm      -f -r deb/fitsh
