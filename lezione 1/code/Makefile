CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe

esercizio1 : esercizio1.o random.o
	$(CC) random.o esercizio1.o -o esercizio1

esercizio1.2 : esercizio1.2.o random.o
	$(CC) random.o esercizio1.2.o -o esercizio1.2

esercizio1.3 : esercizio1.3.o random.o
	$(CC) random.o esercizio1.3.o -o esercizio1.3

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
