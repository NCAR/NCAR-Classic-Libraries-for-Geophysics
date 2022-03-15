
include ../make.inc

SRC=testrvsph.f tdiv.f testrssph.f testsshifte.f testvshifte.f testvtsgs.f \
    tgaqd.f tgrad.f tidvt.f tsha.f tshpe.f tshpg.f tslap.f tvha.f tvlap.f \
    tvrt.f tvts.f

EXES=$(subst .f,.exe, $(SRC))

all : $(EXES)

$(EXES) : $(SRC)

%.exe : %.f
	rm -f $@
	$(F90) $< -o $@ -L../lib -l spherepack
	./$@


clean :
	rm -f $(EXES)
