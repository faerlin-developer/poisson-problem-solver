src = app/main.cpp io/bitmap.cpp io/csv.cpp

all:
	mpicxx -Wall -Werror -o main $(src) -Iapp -Iio

exec:
	mpiexec -n 12 ./main

clean:
	rm -rf main