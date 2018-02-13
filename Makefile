CC = gcc
INCLUDES=-I/opt/X11/include
LDFLAGS=-lX11 -lm -avx
CFLAGS = -Wall -O3 -ffast-math
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
