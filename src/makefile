EXEC = godunov_method
OBJS = main.o godunov.o sim_data.o
SRCS = .
INCS = ./inc
LIBS = ./inc/lib

CC = g++
CFLAGS = -Wall -O3 -c -std=c++0x
LFLAGS = -Wall -O3 -o $(EXEC)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

main.o: $(SRCS)/main.cpp $(INCS)/sim_data.hpp $(INCS)/godunov.hpp
	$(CC) $(CFLAGS) $(SRCS)/main.cpp

sim_data.o: $(SRCS)/sim_data.cpp $(INCS)/sim_data.hpp
	$(CC) $(CFLAGS) $(SRCS)/sim_data.cpp

godunov.o: $(SRCS)/godunov.cpp $(INCS)/godunov.hpp
	$(CC) $(CFLAGS) $(SRCS)/godunov.cpp

clean:
	@rm *.o