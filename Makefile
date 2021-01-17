SRC_DIR = src
OBJ_DIR = obj
TEST_DIR = tests

INC_DIRS =
LIB_DIRS =
LIBSCV = -lsundials_cvode -lsundials_nvecopenmp -lm
LIBSAK = -lsundials_arkode -lsundials_nvecopenmp -lm

.PHONY: tests aprof

a.out: $(SRC_DIR)/main.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -O2 -o $@ $^ $(INC_DIRS) $(LIB_DIRS) $(LIBSAK)

tests: $(TEST_DIR)/tests.out

$(TEST_DIR)/tests.out: $(TEST_DIR)/swave_tests.c $(SRC_DIR)/swave.c $(SRC_DIR)/swave.h
	gcc -o $@ $^ $(INC_DIRS) $(LIB_DIRS) $(LIBSCV)

$(OBJ_DIR):
	mkdir $(OBJ_DIR)
