FPIC      = "-fPIC"

#MODE      =PRODUCTION

# ------------------------------------   comment to avoid optimization; -ggdb -O0 will be used

OPTIONS = -DUSE_QAGS #-DDEBUG

ifeq (${MODE}, PRODUCTION)
MMODE    = PRODUCTION
LIBSN    = libsn.so
LIBSNS   = libsn.a
CFLAGS_1 = -std=c11 -O3 -march=native $(OPTIONS) $(FPIC) -fno-stack-protector 
CFLAGS_2 = -Wall -W -I$(myINCLUDE)
LIBIMF = -lIMF

else

MMODE    = DEBUG
LIBSN    = libsn_dbg.so
LIBSNS   = libsn_dbg.a
CFLAGS_1 = -std=c11 -ggdb3 -Wall -W -I$(myINCLUDE) $(OPTIONS) $(FPIC) #-fsanitize=float-divide-by-zero -fsanitize=bounds -fsanitize=signed-integer-overflow -fsanitize=vla-bound -fsanitize=undefined -fsanitize=integer-divide-by-zero
CFLAGS_2 = -Wall -W -I$(myINCLUDE)
LIBIMF   = -lIMF_dbg

endif

CFLAGS = $(CFLAGS_1) $(CFLAGS_2)

# ------------------------------------ define dependence files for each target
OBJS   = libsn_vars.o libsn.o libsn_io.o libsn_options.o
INCL   = libsn_vars.h Makefile

# ------------------------------------ useful stuffs
LNKCMD = ld -shared $(FPIC)
LNKCMDS= ar rcs
LIBS   = -lgsl -lgslcblas -lm -L$(myLIB) $(LIBIMF)
CC     = gcc

# ------------------------------------ define rules

libsn: $(OBJS)
	@$(LNKCMD) $(OBJS) $(LIBS) -o $(LIBSN)
	@echo "========================  ::  "$(LIBSN)" built"
	@echo "                              with flags " $(CFLAGS)
	@echo "                              with opts  " $(OPTIONS)
	@$(LNKCMDS) $(LIBSNS) $(OBJS)
	@echo "========================  ::  "$(LIBSNS)" built"
	@echo "                              with flags " $(CFLAGS)
	@echo "                              with opts  " $(OPTIONS)
	@echo "indexing the library..."
	@ranlib $(LIBSNS)
	@rm -f libsn_options.c

$(OBJS): $(INCL) libsn_options.c

%.o : %.c
	@echo "........................[$(CC)]" $< "->" $@
	@$(CC) $(OPTIONS) $(CFLAGS) -c $< -o $@

libsn_options.c:
	@echo "char *options_string = \""$(MMODE) $(OPTIONS) $(CFLAGS_1) "\";" > libsn_options.c

clean:
	@echo "cleaning objects files: "
	@echo $(OBJS) 
	@rm -f $(OBJS) 

clean_all:
	@echo "cleaning objects files and libraries: "
	@echo $(OBJS) $(LIBSN) $(LIBSNS)
	@rm -f $(OBJS) $(LIBSN) $(LIBSNS)

