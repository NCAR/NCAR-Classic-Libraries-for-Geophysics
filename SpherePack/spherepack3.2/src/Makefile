
include ../make.inc

SRC=advec.f alf.f divec.f dives.f divgc.f divgs.f gaqd.f geo2math.f \
    gradec.f grades.f gradgc.f gradgs.f helmsph.f hrfft.f idivec.f \
    idives.f idivgc.f idivgs.f idvtec.f idvtes.f idvtgc.f idvtgs.f \
    igradec.f igrades.f igradgc.f igradgs.f ihgeod.f isfvpec.f isfvpes.f \
    isfvpgc.f isfvpgs.f islapec.f islapes.f islapgc.f islapgs.f ivlapec.f \
    ivlapes.f ivlapgc.f ivlapgs.f ivrtec.f ivrtes.f ivrtgc.f ivrtgs.f \
    lfim.f lfin.f lfp.f lfpt.f sfvpec.f sfvpes.f sfvpgc.f sfvpgs.f shaec.f \
    shaes.f shagc.f shags.f shallow.f shigc.f shigs.f shpe.f shpg.f shsec.f \
    shses.f shsgc.f shsgs.f slapec.f slapes.f slapgc.f slapgs.f sphcom.f \
    sshifte.f trssph.f trvsph.f vhaec.f vhaes.f vhagc.f vhags.f vhsec.f \
    vhses.f vhsgc.f vhsgs.f visequ.f visgau.f visgeo.f vlapec.f vlapes.f \
    vlapgc.f vlapgs.f vrtec.f vrtes.f vrtgc.f vrtgs.f vshifte.f vsurf.f \
    vtsec.f vtses.f vtsgc.f vtsgs.f

OBJ=$(subst .f,.o,$(SRC))
OBJS=$(addprefix ../objs/,$(OBJ))

$(LIB) : $(OBJS)
	$(AR) -rv $@ $? 

../objs/%.o : %.f
	$(F90) -c $< -o $@

clean:
	rm -f $(LIB) $(OBJS) 
