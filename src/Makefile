EXECUTABLE = bamrdplot
SAMTOOLS = samtools-0.1.18
ALGLIB = alglib
LOBJS = ${ALGLIB}/libalglib.a ${SAMTOOLS}/libbam.a
CC = g++
CFLAGS = -g -O2 -Wall
#CFLAGS = -g -Wall
INCLUDE = -I${SAMTOOLS} -I${ALGLIB}
LIBS =  -L${SAMTOOLS} -lbam -L${ALGLIB} -lalglib -lz
BACKUPDIR = ~/projects/bamrdplot
BACKUPLOG = ${BACKUPDIR}/bamrdplot.log
STAMP = $(shell  date +'%Y.%m.%d-%H.%M.%S')
BACKUPFOLDER = ${BACKUPDIR}/${STAMP}	

MAINCXX = plotcnvmain.cpp plotcnv.cpp wufunctions.cpp samfunctions.cpp
MAINHDR = wu2.h pstream.h $(MAINCXX:.cpp=.h)
MAINOBJ = $(MAINCXX:.cpp=.o)

MYSOURCE = ${MAINCXX} ${MAINHDR}

all : ${EXECUTABLE}

bamrdplot : $(LOBJS) ${MAINOBJ}
	${CC} ${CFLAGS} ${MAINOBJ} ${LIBS} -o $@

$(MAINOBJ): %.o: %.cpp %.h wu2.h
	${CC} ${INCLUDE} ${CFLAGS} -c $< -o $@

${ALGLIB}/libalglib.a : 
	cd ${ALGLIB} && make

${SAMTOOLS}/libbam.a :
	cd ${SAMTOOLS} && make libbam.a

clean : 
	cd ${SAMTOOLS} && make clean
	cd ${ALGLIB} && make clean
	rm -fr *.o ${EXECUTABLE}

backup :	
	mkdir -p ${BACKUPFOLDER}
	cp Makefile ${MYSOURCE} ${BACKUPFOLDER}
	cp -R ${SAMTOOLS} ${BACKUPFOLDER}
	cp -R ${ALGLIB} ${BACKUPFOLDER}
	if [ -e `echo *.R  | cut -d" " -f1` ]; then cp *.R  ${BACKUPFOLDER}; fi	
	if [ -e `echo *.sh | cut -d" " -f1` ]; then cp *.sh ${BACKUPFOLDER}; fi	
	if [ -e `echo *.pl | cut -d" " -f1` ]; then cp *.pl ${BACKUPFOLDER}; fi	
	if [ -e `echo *.gp | cut -d" " -f1` ]; then cp *.gp ${BACKUPFOLDER}; fi	
	tar -czf ${BACKUPDIR}/${STAMP}.tgz -C ${BACKUPDIR} ${STAMP}
	rm -fR ${BACKUPFOLDER}
	echo ${BACKUPLOG}
	echo -e "${STAMP}.tgz\t${comment}" >> ${BACKUPLOG}
