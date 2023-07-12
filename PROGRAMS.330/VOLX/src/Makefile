#####
# Makefile for Volume X - Computer Programs in Seismology
#
# To move this to another machine, the following variables must be
# defined:
#
#       #Machine Type  - Inform -fative Comment
#       FCMP    - name of FORTRAN Compiler with appropriate flags
#       CCMP    - name of C Compiler with appropriage flags
#       MCHDEP  - Name of object file for machine dependent initialization
#       MCHCMD  - Name of object file for getting number of command
#                       line arguments
#       CALPLTF - Name of CALPLOT graphics library or absolute path,
#               e.g., CALPLTF=-lcalpltf if it in /lib /usr/.lib or
#                       /usr/local/lib
#               or CALPLTF=../../lib/libcalpltf.a
#       DEST    - Destination for the final execuable
#####
#       If one of the following models is appropriate, just uncomment
#       it by removing the first # sign
#
#       Otherwise, use the existing machine models as guides
#####

#####
# LINUX gcc/gfortran -frecord-marker=4 Compilers
#####
FCMP=gfortran -O -frecord-marker=4  
#FCMP=gfortran -pg
CCMP=gcc 
MCHDEP=mchdep.o
MCHCMD=mnmarg.o mgtarg.o ssytem.o
CALPLTF=../../lib/libcalpltf.a
RM= rm -f
MV= mv

DEST=../../bin

  
.c.o:
	$(CCMP) $(CFLAG) -c $<

.f.o:
	$(FCMP) $(FFLAG) -c $<

all:	alls  allt clean

mnmarg.o:
	cp ../../SUBS/mnmarg.f .
	$(FCMP) -c mnmarg.f

mgtarg.o:
	cp ../../SUBS/mgtarg.f .
	$(FCMP) -c mgtarg.f

mchdep.o:
	cp ../../SUBS/mchdep.f .
	$(FCMP) -c mchdep.f

lgstr.o:
	cp ../../SUBS/lgstr.f .
	$(FCMP) -c lgstr.f

igetmod.o:
	cp ../../SUBS/igetmod.f .
	$(FCMP) -c igetmod.f

iputmod.o:
	cp ../../SUBS/iputmod.f .
	$(FCMP) -c iputmod.f

sacsubf.o:
	cp ../../SUBS/sacsubf.f .
	$(FCMP) -c sacsubf.f

grphsubf.o:
	cp ../../SUBS/grphsubf.f .
	$(FCMP) -c grphsubf.f


#####
#	travel time  inversion routines
#####

#allt: ttinv96  ttinvdr96   ttinvpv96  ttinvvp ttinvpr96
allt:   ttinvdr96   ttinvpv96  ttinvvp ttinvpr96 ttinv96

ttinv96: ttinv96.o modlt.o tmenu.o ssamat.o stamat.o tamat.o igetmod.o iputmod.o \
		query.o lgstr.o shwctl.o shwwtv.o $(MCHDEP) $(MCHCMD)  
	$(FCMP) $(FFLAG) -o ttinv96 ttinv96.o  modlt.o tmenu.o ssamat.o stamat.o \
		query.o tamat.o \
		igetmod.o iputmod.o lgstr.o shwctl.o shwwtv.o $(MCHDEP) $(MCHCMD) 
	$(MV) ttinv96 $(DEST)


ttinvdr96: ttinvdr96.o lgstr.o sacsubf.o igetmod.o iputmod.o  \
		$(MCHDEP) $(MCHCMD)
	$(FCMP) $(FFLAG) -o ttinvdr96 ttinvdr96.o lgstr.o sacsubf.o \
		igetmod.o iputmod.o  $(MCHDEP) $(MCHCMD)
	$(MV) ttinvdr96 $(DEST)


ttinvpv96:	ttinvpv96.o igetmod.o grphsubf.o lgstr.o $(MCHDEP) $(MCHCMD)
	$(FCMP) $(FFLAG) -o ttinvpv96 ttinvpv96.o igetmod.o grphsubf.o lgstr.o \
		$(MCHDEP) $(MCHCMD) $(CALPLTF)
	$(MV) ttinvpv96 $(DEST)

ttinvpr96: ttinvpr96.o sacsubf.o  lgstr.o igetmod.o iputmod.o  setmod.o $(MCHDEP)
	$(FCMP) $(FFLAG) -o ttinvpr96 ttinvpr96.o sacsubf.o  lgstr.o \
		igetmod.o iputmod.o  setmod.o $(MCHDEP)
	$(MV) ttinvpr96 $(DEST)
	$(RM) ttinvpr96.o

ttinvvp: Makefile
	cp ttinvvp $(DEST)
	touch Makefile

#####
#	shallow joint inversion
#####

alls: shallow96 shalpr96

shalpr96: shalpr96.o sacsubf.o  lgstr.o igetmod.o iputmod.o setdsp.o setmod.o $(MCHDEP)
	$(FCMP) $(FFLAG) -o shalpr96 shalpr96.o sacsubf.o  lgstr.o \
		igetmod.o iputmod.o setdsp.o setmod.o $(MCHDEP)
	$(MV) shalpr96 $(DEST)
	$(RM) shalpr96.o

shallow96: shallow96.o modlt.o shmenu.o ssamat.o stamat.o shamat.o igetmod.o iputmod.o \
		query.o lgstr.o shwctl.o shwwtv.o $(MCHDEP) $(MCHCMD)  
	$(FCMP) $(FFLAG) -o shallow96 shallow96.o  modlt.o shmenu.o ssamat.o stamat.o \
		query.o shamat.o \
		igetmod.o iputmod.o lgstr.o shwctl.o shwwtv.o $(MCHDEP) $(MCHCMD) 
	$(MV) shallow96 $(DEST)


#####
#	clean up a work files
#####

clean:
	$(RM) *.o
	$(RM) grphsubf.f
	$(RM) igetmod.f
	$(RM) iputmod.f
	$(RM) lgstr.f
	$(RM) mnmarg.f
	$(RM) mgtarg.f
	$(RM) mchdep.f
	$(RM) sacsubf.f

