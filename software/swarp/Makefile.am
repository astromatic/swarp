# Main Makefile for SWarp
# Copyright (C) 2003-2007 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
dist_pkgdata_DATA	= doc xsl/swarp.xsl
EXTRA_DIST		= AUTHORS BUGS ChangeLog HISTORY INSTALL README \
			  THANKS TODO acx_pthread.m4 acx_prog_cc_optim.m4 \
			  acx_mpi.m4 test/test.fits acx_urbi_resolve_dir.m4
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -march=opteron" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec                                                               

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -axWP -ip -unroll" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec                                                               

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -axWP -ip -unroll" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

icc:
	$(MAKE) CC="icc" CFLAGS="-O3 -axWP -ip -unroll -pthread" LIBS="-lm -lpthread"

debug:
	$(MAKE) CFLAGS="-O3 -funroll-loops -fomit-frame-pointer -Wall -m32 -g" LIBS="-L/usr/lib/ -lpthread -lm"
