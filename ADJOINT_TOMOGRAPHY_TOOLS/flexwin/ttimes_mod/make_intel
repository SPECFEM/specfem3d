OPT=-g 
FC=ifort
FFLAGS=$(OPT)

CODE=${HOME}/code

BINDIR=${CODE}/bin
SHAREDIR=${CODE}/share

LIBDIR=${CODE}/lib

LIBSRC=libtau.f libsun.f emiasp91.f
LIBSRC_SOCAL=libtau.f libsun.f emsocal.f
LIB=libtau.a
LIB_SOCAL=libtau_socal.a
LIBOBJ = ${LIBSRC:%.f=${LIB}(%.o)}
LIBOBJ_SOCAL = ${LIBSRC_SOCAL:%.f=${LIB_SOCAL}(%.o)}

INSTALL = install

PROGS= remodl remodl_socal setbrn setbrn_socal 

default: ${PROGS} iaspei socal

SUFFIXES= .f .o
.SUFFIXES: ${SUFFIXES}

.f.o:
	${FC} ${FFLAGS} -c -o $@ $<

${LIB}:${LIBOBJ}
${LIB_SOCAL}:${LIBOBJ_SOCAL}

remodl:remodl.o limits.inc ${LIB}
	${FC} ${FFLAGS} -o $@ remodl.o ${LIB} 

setbrn:setbrn.o ${LIB}
	${FC} ${FFLAGS} -o $@ setbrn.o  ${LIB}

remodl_socal:remodl.o limits.inc ${LIB_SOCAL}
	${FC} ${FFLAGS} -o $@ remodl.o ${LIB_SOCAL} 

setbrn_socal:setbrn.o ${LIB_SOCAL}
	${FC} ${FFLAGS} -o $@ setbrn.o  ${LIB_SOCAL}


ttimes:ttimes.o ttlim.inc ${LIB}
	${FC} ${FFLAGS} -o $@ $@.o ${LIB} 

findit:findit.o 
	${FC} ${FFLAGS} -o $@ $@.o 

pick_index:pick_index.o 
	${FC} ${FFLAGS} -o $@ $@.o -L${LIBDIR} -lsac 

iaspei: remodl setbrn
	./remodl
	./setbrn

socal: remodl_socal setbrn_socal
	./remodl_socal
	./setbrn_socal

clean:
	rm -f *.o ${LIB} ${LIB_SOCAL} ${PROGS} *.hed *.tbl *.lis

install:
	for i in ${PROGS}; do ${INSTALL} -m 0755 -s $$i $(BINDIR); done
	for i in SEED_prep.sh data_prep.sh iasp_prep.sh run_rdseed ; do ${INSTALL} -m 0755 $$i $(BINDIR); done
	${INSTALL} -c -m 0644 iasp91.hed iasp91.tbl resp_wwlpbn $(SHAREDIR)
	${INSTALL} -c -m 0644 $(LIB) $(LIBDIR)
	ranlib  $(LIBDIR)/${LIB}

	
