# list your tests here
TESTS = tests/a5.test tests/m11.test

# the rest stays the same when you just add a test
bin_PROGRAMS = tcmain tcfrontend
dist_bin_SCRIPTS = bin/tcenum
tcmain_SOURCES = src/tcmain.F
tcfrontend_SOURCES = src/tcfrontend.c

TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/build-aux/tap-driver.sh

EXTRA_DIST = $(TESTS) tests/test.sh $(TESTS:.test=.res) $(TESTS:.test=)
