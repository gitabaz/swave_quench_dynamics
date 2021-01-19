CXX=clang++-11
#CXX=g++
SRC_DIR = src
OBJ_DIR = obj
TEST_DIR = tests

INC_DIRS =
LIB_DIRS =
LIBSAK = -lsundials_arkode -lsundials_nvecopenmp -lm
#LIBSCV = -lsundials_cvode -lsundials_nvecopenmp -lm

.PHONY: tests aboost

a.out: $(SRC_DIR)/main.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -O2 -o $@ $^ $(INC_DIRS) $(LIB_DIRS) $(LIBSAK)

aboost: aboost.out

aboost.out: $(SRC_DIR)/boost/main.cpp $(SRC_DIR)/boost/swave.cpp #$(SRC_DIR)/boost/swave.hpp
	$(CXX) -fopenmp -O2 $^ -I/home/olabian/opt/boost_1_75_0 -o $@

tests: $(TEST_DIR)/tests.out

$(TEST_DIR)/tests.out: $(TEST_DIR)/swave_tests.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -o $@ $^ $(INC_DIRS) $(LIB_DIRS) $(LIBSCV)

$(OBJ_DIR):
	mkdir $(OBJ_DIR)
