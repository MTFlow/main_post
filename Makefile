# Makefile to generate executable

CC = g++
CFLAGS = -O3 -g -Wall 
DEPS = main.h func.h 
OBJ	= main.o

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

main: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ 

.PHONY: clean

clean:
	$(RM) main *.o *~ 
