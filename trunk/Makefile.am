# Main Makefile for SWarp
# Copyright (C) 2003-2010 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
dist_pkgdata_DATA	= doc xsl/swarp.xsl
EXTRA_DIST		= AUTHORS BUGS ChangeLog COPYRIGHT HISTORY INSTALL \
			  LICENSE README THANKS TODO \
			  acx_pthread.m4 acx_prog_cc_optim.m4 \
			  acx_urbi_resolve_dir.m4 test/test.fits
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
dist-hook:
	rm -rf `find $(distdir) -name .svn`

rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

rpm-best:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_BEST="1" rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

