#

all: libalchemify.so

alchemify_tcl_wrap.c : alchemify_tcl.c libalchemify.c alchemify_tcl.i
	swig -tcl alchemify_tcl.i

libalchemify.so : alchemify_tcl_wrap.c alchemify_tcl.c libalchemify.c
	gcc -ansi -c -fpic alchemify_tcl.c alchemify_tcl_wrap.c libalchemify.c -I/usr/local/include
	gcc -ansi -shared  alchemify_tcl.o alchemify_tcl_wrap.o libalchemify.o -o libalchemify.so

