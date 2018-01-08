#
#				configure.ac
#
# Main Makefile. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SWarp
#
#	Copyright:		(C) 2002-2018 IAP/CNRS/UPMC
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
#	Last modified:		08/01/2018
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACLOCAL_AMFLAGS		= -I m4
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
dist_pkgdata_DATA	= xsl/swarp.xsl
EXTRA_DIST		= doc AUTHORS BUGS ChangeLog COPYRIGHT HISTORY INSTALL \
			  LICENSE README.md THANKS \
			  acx_pthread.m4 acx_prog_cc_optim.m4 \
			  acx_urbi_resolve_dir.m4 test/test.fits
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
DEB_SRCDIR		= ~/DEBIAN/SOURCES

dist-hook:
	rm -rf `find $(distdir) -name .svn`

rpm:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

rpm-best:	dist
	cp -f $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_BEST="1" rpmbuild -ba --clean --nodeps $(PACKAGE_TARNAME).spec

deb-best:	dist
	mkdir -p $(DEB_SRCDIR)
	tar -xvf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz -C $(DEB_SRCDIR) 
	cd $(DEB_SRCDIR)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)
	debuild  --set-envvar=PACKAGE_RELEASE=$(PACKAGE_RELEASE) --set-envvar=USE_BEST=1 --preserve-env --preserve-envvar=PATH -j4 -b -us -uc

