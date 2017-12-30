OBJECTS = buildproblemconj.o main.o
EIGEN_DIR=./dependencies/Eigen
#compiler to use, currently requires users to update to their desired one 

CC = g++-7

NAME = FEM_Laplace


all: $(NAME) $(OBJECTS)

$(NAME):$(OBJECTS)
	$(CC) -o $(NAME) $(OBJECTS)

buildproblemconj.o: buildproblemconj.cpp buildproblemconj.hpp
	$(CC) -c -I$(EIGEN_DIR) buildproblemconj.cpp 

main.o: buildproblemconj.o main.cpp
	$(CC) -c -I$(EIGEN_DIR) main.cpp
.PHONY : clean

clean:
	rm -f $(NAME) $(OBJECTS)

