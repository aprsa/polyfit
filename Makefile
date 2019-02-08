all: clean polyfit

polyfit:  polyfit.c
	gcc -g -Wall -O2 -I/usr/include polyfit.c -lgsl -lgslcblas -lm -o polyfit

clean:
	rm -f polyfit

install:
	cp polyfit /usr/local/bin

dist:
	cd ..; tar cvzf polyfit-0.12.tar.gz polyfit-0.12/polyfit.c polyfit-0.12/Makefile polyfit-0.12/lc.dat polyfit-0.12/README; cd polyfit-0.12
