# Compilatore e opzioni di compilazione
CC = mpic++
CFLAGS = -Wall -O3 -std=c++11

# Nome dell'eseguibile finale
TARGET = tsp

# Oggetti da compilare
OBJECTS = main.o tsp.o random.o

# Regola generale di compilazione
all: $(TARGET)

# Regola per l'eseguibile finale
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET)

# Regola per compilare il file sorgente del programma principale
main.o: main.cpp tsp.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

# Regola per compilare il file sorgente del TSP
tsp.o: tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)

# Regola per compilare il file sorgente del generatore casuale
random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

# Pulizia degli oggetti e dell'eseguibile
clean:
	rm -f $(OBJECTS) $(TARGET) seed.out

# Regola per rigenerare il seed.out se necessario
seed.out:
	# Inserisci qui i comandi per generare seed.out, se necessario

.PHONY: all clean


