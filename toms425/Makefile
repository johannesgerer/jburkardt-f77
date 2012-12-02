F77=nagfor
F77OPTS=-gline -u -dcfuns -nan -save -C -C=undefined
all: Res

src.o: src.f
	$(F77) $(F77OPTS) -c src.f

driver.o: driver.f
	$(F77) $(F77OPTS) -c driver.f

DRIVERS= driver
RESULTS= Res

Objs1= driver.o src.o
driver: $(Objs1)
	$(F77) $(F77OPTS) -o driver $(Objs1) $(SRCLIBS)
Res: driver 
	./driver >Res

diffres:Res res
	echo "Differences in results from driver"
	$(DIFF) Res res

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
