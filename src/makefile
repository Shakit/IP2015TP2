# définition des actions
.PHONY: clean

# définition des variables
CC = gcc
CCFLAGS =  
LDFLAGS = -lglpk -lm
SOURCES = main.c
OBJECTS = $(SOURCES:.c =.o)

# regles
all: $(OBJECTS)
	$(CC) -o ../bin/solveLotSizing $^ $(CCFLAGS) $(LDFLAGS)

%.o: %.c
	$(CC) -c $< $(CCFLAGS) $(LDFLAGS)

# action clean
clean:
	-rm *.o ../bin/solveLotSizing ../res/*.lp *~
