bin_PROGRAMS = anapipe anarpt

anapipe_LDADD = $(LDFLAGS)
anapipe_SOURCES = main.cpp genjob.cpp genpipe.cpp job.cpp options.cpp pipeline.cpp task.cpp

anarpt_LDADD = $(LDFLAGS)
anarpt_SOURCES = anarpt.cpp

clean:
	rm -rf .deps Makefile.in Makefile *.o ${bin_PROGRAMS}
