default:
	mpicc -Wall -g -std=c99 main.c -o main

profile:
	mpicc -Wall -g -pg -std=c99 main.c -o main

run:
	mpirun -np 2 main s5e50.in 1

debug:
	mpirun -np 2 xterm -e ddd ./main
