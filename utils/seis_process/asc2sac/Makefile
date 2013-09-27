
CC=gcc
CFLAGS= -DSAC_NULL
PROG=asc2sac
SRCS=asc2sac.c  sacclose.c  sacdatawrite.c  sachdrwrite.c  sacopen.c
OBJS=${SRCS:%.c=%.o}

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *~ core

