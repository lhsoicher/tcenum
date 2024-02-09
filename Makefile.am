bin_PROGRAMS = tcmain tcfrontend
dist_bin_SCRIPTS = bin/tcenum
tcmain_SOURCES = src/tcmain.F
tcfrontend_SOURCES = src/tcfrontend.c

TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/build-aux/tap-driver.sh

TESTS = tests/a5.test tests/m11.test
EXTRA_DIST = $(TESTS) tests/test.sh tests/a5 tests/a5.res tests/m11 tests/m11.res
