SUBDIRS = . surflib registration tools 
#doc

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = $(m4_files) \
	HACKING \
	test/tester.cc test/tester_meshes.hpp

m4_files = \
        m4/mni_REQUIRE_LIB.m4 \
        m4/mni_REQUIRE_MNILIBS.m4 \
        m4/mni_REQUIRE_OPENINVENTOR.m4 \
        m4/mni_cxx_have_koenig_lookup.m4 \
        m4/smr_OPTIONAL_LIB.m4 \
        m4/smr_REQUIRED_LIB.m4 \
        m4/smr_WITH_BUILD_PATH.m4 \
	m4/smr_CGAL_MAKEFILE.m4
