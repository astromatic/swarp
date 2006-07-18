# Main Makefile for SWarp
# Copyright (C) 2003-2005 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
EXTRA_DIST		= doc AUTHORS BUGS ChangeLog HISTORY INSTALL README \
			  THANKS TODO acx_pthread.m4 acx_prog_cc_optim.m4 \
			  acx_mpi.m4 test/test.fits
RPM_ROOTDIR		= /usr/src/RPM
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

rpm-opteron:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
#	ICCFLAGS="-O3 -xW" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec
	GCCFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -march=opteron" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec                                                               

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
#	ICCFLAGS="-O3 -xK" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec
	GCCFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec $(PACKAGE_NAME)-mp.spec

debug:
	$(MAKE) CFLAGS="-O3 -funroll-loops -fomit-frame-pointer -Wall -m32 -g" LIBS="-L/usr/lib/ -lpthread -lm"
