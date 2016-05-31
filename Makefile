all: npstat


npstat: 
	gcc -o npstat NPStat-v1.c -lgsl -lgslcblas -lm

clean: 
	rm npstat