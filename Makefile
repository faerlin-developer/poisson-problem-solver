main = src/main.cpp

all:
	mpicxx -Wall -Werror -o main $(main) -Iinclude

exec:
	mpiexec -n 12 ./main


