
CLEANFILES= src/*.lo src/*.o src/starscreamTest initial_conditions.dat

include_HEADERS=starscream.h
lib_LTLIBRARIES=libstarscream.la
libstarscream_la_SOURCES=starscream_init.c starscream_io.c starscream_pf.c \
                         starscream_structure.c starscream_vel.c 

check_PROGRAMS=starscreamTest

starscreamTest:
	$(CC) $(AM_CFLAGS) -o starscreamTest starscream.c \
        -L@top_builddir@/src/.libs -lstarscream $(GSL_LIBS) $(FFTW_LIBS)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):@top_builddir@/src/.libs:$(GSL_LIBDIR):$(FFTW_LIBDIR) \
        time ./starscreamTest
