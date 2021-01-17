SRC_DIR = src
OBJ_DIR = obj
TEST_DIR = tests

INC_DIRS =
LIB_DIRS =
#LIBS = -lsundials_cvode -lsundials_nvecopenmp -lm
LIBS = -lsundials_arkode -lsundials_nvecopenmp -lm

.PHONY: tests aprof

a.out: $(SRC_DIR)/main.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -O2 -o $@ $^ $(INC_DIRS) $(LIB_DIRS) $(LIBS)

aprof: $(SRC_DIR)/main.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -g -o aprof.out $^ $(INC_DIRS) $(LIB_DIRS) $(LIBS)

tests: $(TEST_DIR)/tests

$(TEST_DIR)/tests: $(TEST_DIR)/swave_tests.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -o $@ $^ $(INC_DIRS) $(LIB_DIRS) $(LIBS)

$(OBJ_DIR)/swave.o: $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h | $(OBJ_DIR)
	gcc -c $(SRC_DIR)/swave.c -o $@

$(OBJ_DIR):
	mkdir $(OBJ_DIR)
