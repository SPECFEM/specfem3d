FC      = mpif77
CC	= mpicc
FFLAGS  = -O
CFLAGS	= $(FFLAGS)
OBJS    = param.o matrix.o compwave.o msub.o sacconv.o sacbw.o dsmti.o


all:
	make dsmti

dsmti: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o dsmti

.f.o:
	$(FC) -c $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f $(OBJS) dsmti *~
