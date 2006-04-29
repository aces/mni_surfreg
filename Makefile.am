SUBDIRS = . surflib registration tools docs

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = $(m4_files)					       \
	epm-header.in          	                                       \
	HACKING							       \
	test/tester.cc test/tester_meshes.hpp			       \
	test/poly_20.obj test/poly_20_flipx.sm			       \
	test/poly_20_identity.sm test/poly_5120.obj		       \
	test/poly_5120_alpha.vv test/poly_5120_chull.vv		       \
	test/poly_5120_dt00.vv test/poly_5120_dt05.vv		       \
	test/poly_5120_dt15.vv test/poly_5120_sdt00.vv		       \
	test/poly_5120_sdt05.vv test/poly_5120_sdt15.vv

m4_files = \
        m4/mni_REQUIRE_LIB.m4 \
        m4/mni_REQUIRE_MNILIBS.m4 \
        m4/mni_REQUIRE_OPENINVENTOR.m4 \
        m4/mni_cxx_have_koenig_lookup.m4 \
        m4/smr_OPTIONAL_LIB.m4 \
        m4/smr_REQUIRED_LIB.m4 \
        m4/smr_WITH_BUILD_PATH.m4 \
	m4/smr_CGAL_MAKEFILE.m4
