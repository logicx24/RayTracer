CC = g++
RM = /bin/rm -f 

all: main
main: raytracer.o
	$(CC) -o raytracer raytracer.o lodepng.o
raytracer.o: lodepng raytracer.cpp
	$(CC) -c raytracer.cpp -o raytracer.o
lodepng:
	$(CC) -c lodepng.cpp -o lodepng.o
clean:
	$(RM) *.o raytracer