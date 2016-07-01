CC=g++
EXE=
RM=rm -f

ifdef CLANG
	CC=clang
endif

ifdef MSWIN
	EXE=.exe
endif

ifdef XCOMPILE
	CC=x86_64-w64-mingw32-g++
	EXE=.exe
endif

all: get_high$(EXE) mergetle$(EXE) obs_tes2$(EXE) obs_test$(EXE) out_comp$(EXE) \
	test_sat$(EXE) test2$(EXE) sat_id$(EXE) sat_id2$(EXE) test_out$(EXE)

CFLAGS=-Wextra -Wall -O3 -pedantic -Wno-unused-parameter

clean:
	$(RM) *.o
	$(RM) get_high$(EXE)
	$(RM) mergetle$(EXE)
	$(RM) obs_tes2$(EXE)
	$(RM) obs_test$(EXE)
	$(RM) out_comp$(EXE)
	$(RM) libsatell.a
	$(RM) sat_id$(EXE)
	$(RM) sat_id2$(EXE)
	$(RM) test2$(EXE)
	$(RM) test_out$(EXE)
	$(RM) test_sat$(EXE)

install:
	cp libsatell.a /usr/local/lib
	cp norad.h     /usr/local/include

uninstall:
	rm /usr/local/lib
	rm /usr/local/include

OBJS= sgp.o sgp4.o sgp8.o sdp4.o sdp8.o deep.o basics.o get_el.o common.o tle_out.o

get_high$(EXE):	 get_high.o get_el.o
	$(CC) $(CFLAGS) -o get_high$(EXE) get_high.o get_el.o

mergetle$(EXE):	 mergetle.o
	$(CC) $(CFLAGS) -o mergetle$(EXE) mergetle.o

obs_tes2$(EXE):	 obs_tes2.o observe.o libsatell.a
	$(CC) $(CFLAGS) -o obs_tes2$(EXE) obs_tes2.o observe.o libsatell.a -lm

obs_test$(EXE):	 obs_test.o observe.o libsatell.a
	$(CC) $(CFLAGS) -o obs_test$(EXE) obs_test.o observe.o libsatell.a -lm

out_comp$(EXE):	 out_comp.o
	$(CC) $(CFLAGS) -o out_comp$(EXE) out_comp.o -lm

libsatell.a: $(OBJS)
	rm -f libsatell.a
	ar rv libsatell.a $(OBJS)

sat_id$(EXE):	 	sat_id.o	observe.o libsatell.a
	$(CC) $(CFLAGS) -o sat_id$(EXE) sat_id.o observe.o libsatell.a -lm

sat_id2$(EXE):	 	sat_id2.o sat_id.cpp observe.o  ../find_orb/cgi_func.cpp libsatell.a
	$(CC) $(CFLAGS) -o sat_id2$(EXE) -DON_LINE_VERSION sat_id2.o sat_id.cpp observe.o ../find_orb/cgi_func.cpp libsatell.a -lm

test2$(EXE):	 	test2.o sgp.o libsatell.a
	$(CC) $(CFLAGS) -o test2$(EXE) test2.o sgp.o libsatell.a -lm

test_out$(EXE):	 test_out.o tle_out.o get_el.o sgp4.o common.o
	$(CC) $(CFLAGS) -o test_out$(EXE) test_out.o tle_out.o get_el.o sgp4.o common.o -lm

test_sat$(EXE):	 test_sat.o libsatell.a
	$(CC) $(CFLAGS) -o test_sat$(EXE) test_sat.o libsatell.a -lm

.cpp.o:
	$(CC) $(CFLAGS) -c $<
