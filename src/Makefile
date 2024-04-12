GCC=g++
STD_FLAGS=-Wall -Werror -Wextra 
# -std=c11
LIB=-lcheck
GCOV_FLAGS=-fprofile-arcs -ftest-coverage
HTML=lcov -t test -o rep.info -c -d ./
CFLAGS=--std=c++17 -lstdc++ -lm
TESTFLAGS=-lgtest
TEST_FILE = tests/tests.cc

TEST_FOLDER=tests/
FUNCTIONS_FOLDER=functions/

LIB_NAME = s21_matrix_oop.a

UNIT_TEST=-o unit_test
GCOV_TEST=-o gcov_test

ifeq ($(shell uname), Linux)
 LIB +=-lpthread -lrt -lsubunit -lm
endif


all: ${LIB_NAME} test

${LIB_NAME}:
	${GCC} -std=c++17 ${STD_FLAGS} -c ${wildcard ${FUNCTIONS_FOLDER}*.cc}
	ar rc ${LIB_NAME} *.o
	ranlib ${LIB_NAME}
	rm *.o

clean: 
	rm -rf report a.out unit_test gcov_test *.gcda *.gcno *.info ${LIB_NAME} *.o


test: clean ${LIB_NAME}
	$(GCC) -g ${wildcard ${TEST_FOLDER}*.cc} ${LIB_NAME} $(CFLAGS) $(TESTFLAGS) -o test
	./test
	

valgrind_test: ${LIB_NAME}
	${GCC} ${STD_FLAGS} ${UNIT_TEST} ${wildcard ${TEST_FOLDER}*.cc} ${LIB_NAME}
	CK_FORK=no valgrind --trace-children=yes --track-fds=yes --leak-check=full --track-origins=yes --log-file=valgrind-out.txt --show-leak-kinds=all ./unit_test
	rm ${LIB_NAME}
	rm unit_test
	
gcov_report:
	${GCC} -lgtest -lstdc++ --coverage ${wildcard ${TEST_FOLDER}*.cc} ${wildcard ${FUNCTIONS_FOLDER}*.cc} ${GCOV_TEST} 
	./gcov_test
	lcov -t "stest" -o s21_test.info -c -d .
	genhtml -o report s21_test.info
	open ./report/index.html
	rm gcov_test


check_clang_format:
	mv ../materials/linters/.clang-format .
	clang-format -n tests/*.cc
	clang-format -n tests/*.h
	clang-format -n s21_matrix_oop.h
	mv .clang-format ../materials/linters/

fix_clang_format:
	mv ../materials/linters/.clang-format .
	clang-format -i functions/*.cc
	clang-format -i tests/*.cc
	clang-format -i tests/*.h
	clang-format -i s21_matrix_oop.h
	mv .clang-format ../materials/linters/
