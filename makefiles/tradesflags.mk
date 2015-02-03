
CFLAGS = $(CFLAGS1) $(OPTFLAG) -I$(INCPATH) $(VIBFLAG) $(DEF)
LDFLAGS = -I$(INCPATH) $(LDFLAGS1) $(OPTFLAG) -L$(LIBPATH) $(VIBFLAG) $(DEF)


## To clean out the directory without removing make
##

## Implicit actions
##
## if need a .o, compile the .c
##

.c.o :
	$(CC) $(CFLAGS) $<
	mv $(@F) $(SRCDIR)
