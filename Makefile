CC = g++
RM = /bin/rm -f 

all: main
main: raytracer.o lodepng.o
	$(CC) -o raytracer raytracer.o lodepng.o
raytracer.o:
	$(CC) -c raytracer.cpp -o raytracer.o
lodepng.o:
	$(CC) -c lodepng.cpp -o lodepng.o
clean:
	$(RM) *.o raytracer lodepng