CC = gcc
INCLUDES=-I/opt/X11/include
LDFLAGS=-lX11 -lm -lrt
CFLAGS = -Wall -O3 -funroll-loops
RM = rm -f
EXECUTABLE = galsim

$(EXECUTABLE): galsim.o graphics.o
	$(CC) -o $(EXECUTABLE) galsim.o graphics.o $(LDFLAGS)

graphics.o: graphics.c graphics.h
	$(CC) $(CFLAGS) $(INCLUDES) -c graphics.c

galsim.o: galsim.c
	$(CC) $(CFLAGS) $(INCLUDES) -c galsim.c

clean: 
	$(RM) $(EXECUTABLE) *.o
