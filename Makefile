src = app/main.cpp bitmap/bitmap.cpp util/util.cpp

all:
	mpicxx -Wall -Werror -o main $(src) -Ibitmap -Iapp -Iutil

exec:
	mpiexec -n 12 ./main

clean:
	rm -rf main