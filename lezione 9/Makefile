# Compilatore e flag
CC = g++
CFLAGS = -Wall -O3 --std=c++11

# Nome dell'eseguibile
EXEC = genetic_alg

# File sorgente e oggetti
SRCS = main.cpp genetic_alg.cpp random.cpp
OBJS = $(SRCS:.cpp=.o)

# Regole principali
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) -g -o $(EXEC)

main.o: main.cpp genetic_alg.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

genetic_alg.o: genetic_alg.cpp genetic_alg.h random.h
	$(CC) -c genetic_alg.cpp -o genetic_alg.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

# Pulizia
clean:
	rm -f *.o $(EXEC) seed.out
