FFORT = gfortran
FFLAGS = -O3
INC = -I/home/dboateng/applications/netcdf_f/include
LIB = -L/home/dboateng/applications/netcdf_f/lib -lnetcdff -lnetcdf
OBJ = times.o iotra.o ioinp.o inter.o libcdfio.o libcdfplus.o 
OBJS_select = special.o select.o times.o iotra.o ioinp.o inter.o libcdfio.o libcdfplus.o
EXE = caltra create_startf density changet datelist difference extract getmima gettidiff getvars list2lsl lsl2list lsl2rdf mergetra newtime reformat timeres trainfo lidar trace select

all: $(EXE)

.PHONY: clean
clean:
	rm -f *.o *.mod $(EXE)

$(OBJ): %.o: %.f
	$(FFORT) $(FFLAGS) $(INC) -c $<



caltra: caltra.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

create_startf: create_startf.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

density: density.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

changet: changet.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

datelist: datelist.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

difference: difference.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

extract: extract.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

getmima: getmima.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

gettidiff: gettidiff.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

getvars: getvars.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

list2lsl: list2lsl.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

lsl2list: lsl2list.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

lsl2rdf: lsl2rdf.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

mergetra: mergetra.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

newtime: newtime.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

reformat: reformat.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

timeres: timeres.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

trainfo: trainfo.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

lidar: lidar.f $(OBJ)
	$(FFORT) $(FFLAGS) $(INC) -o $@ $< $(OBJ) $(LIB)

trace:	$(OBJ)
	${FFORT} -o trace trace.f90 calvar.f $(OBJ) ${INC} $(LIB)

select: $(OBJS_select)
	$(FFORT) -o select $(OBJS_select) ${INC} $(LIB)




