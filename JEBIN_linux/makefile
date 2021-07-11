CC = g++
CFLAGS = -lm -pthread -Ofast -march=native -Wall -funroll-loops -Wno-unused-result -lgsl -lm -lgslcblas
LFLAGS = -lgsl -lm -lgslcblas
INCLUDES = -I/usr/local/include -I/home/gywu/LINE/env/eigen_install/include/eigen3
LIBS = -L/usr/local/lib



JEBIN_linux : linelib.o main.o
	$(CC) -g $(CFLAGS) -o JEBIN_linux linelib.o main.o $(INCLUDES) $(LIBS) $(LFLAGS)

linelib.o : linelib.cpp
	$(CC) -g $(CFLAGS) -c linelib.cpp $(INCLUDES) $(LIBS) $(LFLAGS)

main.o : main.cpp linelib.o
	$(CC) -g $(CFLAGS) -c main.cpp $(INCLUDES) $(LIBS) $(LFLAGS)

clean :
	rm -rf *.o JEBIN_linux
