all:
	cd src && make
	cp src/bamrdplot .

clean:
	cd src && make clean
	rm bamrdplot 
