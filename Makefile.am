AUTOMAKE_OPTIONS = foreign
SUBDIRS = htslib src

install:  
	mkdir -p lib && cp src/libseqlib.a htslib/libhts.a lib

clean:
	rm -f lib/*.a src/*.o src/*.a
	(cd htslib && make clean)
