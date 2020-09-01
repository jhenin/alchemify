# 

all: alchemify

alchemify: alchemify.c libalchemify.c libalchemify.h
	gcc -c -ansi alchemify.c libalchemify.c
	gcc -o alchemify alchemify.o libalchemify.o

clean:
	rm -f *.o alchemify


