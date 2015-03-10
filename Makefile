CC = g++
RM = /bin/rm -f 

all: clean main 
main: raytracer.o lodepng.o
	$(CC) -o raytracer raytracer.o lodepng.o
raytracer.o:
	$(CC) -c raytracer.cpp -o raytracer.o
lodepng.o:
	$(CC) -c lodepng.cpp -o lodepng.o
file:
	./raytracer > ~/Desktop/out.txt
clean:
	$(RM) *.o raytracer lodepng img.png