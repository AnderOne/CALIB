INCPATH = ./inc
SRCPATH = ./src
OBJPATH = ./obj
LIBPATH = ./lib
BINPATH = ./bin

CALIB_INCLIST = $(wildcard $(INCPATH)/*.hpp  $(INCPATH)/*/*.hpp)
CALIB_SRCLIST = $(wildcard $(SRCPATH)/calib/*.cpp)
CALIB_OBJLIST = $(patsubst $(SRCPATH)/%.cpp, $(OBJPATH)/%.o, $(CALIB_SRCLIST))
CALIB = $(LIBPATH)/libCALIB.a

TEST_INCLIST = $(wildcard $(CALIB_INCLIST) $(SRCPATH)/test/*.hpp)
TEST_SRCLIST = $(wildcard $(SRCPATH)/test/*.cpp)
TEST_OBJLIST = $(patsubst $(SRCPATH)/%.cpp, $(OBJPATH)/%.o, $(TEST_SRCLIST))
TEST = $(patsubst $(OBJPATH)/test/%.o, $(BINPATH)/test/%, $(TEST_OBJLIST))

DEMO_INCLIST = $(wildcard $(CALIB_INCLIST) $(SRCPATH)/demo/*.hpp)
DEMO_SRCLIST = $(wildcard $(SRCPATH)/demo/*.cpp)
DEMO_OBJLIST = $(patsubst $(SRCPATH)/%.cpp, $(OBJPATH)/%.o, $(DEMO_SRCLIST))
DEMO = $(patsubst $(OBJPATH)/demo/%.o, $(BINPATH)/demo/%, $(DEMO_OBJLIST))

OBJLIST = $(CALIB_OBJLIST) $(TEST_OBJLIST) $(DEMO_OBJLIST)

INCLUDE = -I $(INCPATH)
LINKOPT = -L $(LIBPATH) -lCALIB -lblitz -lfftw3\
          -lboost_unit_test_framework -static
CC = g++ -std=c++14 -O3
AR = ar rcs

$(OBJPATH)/calib/%.o: $(SRCPATH)/calib/%.cpp $(CALIB_INCLIST)
	$(CC) -c $< $(INCLUDE) -o $@

$(OBJPATH)/test/%.o: $(SRCPATH)/test/%.cpp $(TEST_INCLIST)
	$(CC) -c $< $(INCLUDE) -o $@

$(OBJPATH)/demo/%.o: $(SRCPATH)/demo/%.cpp $(DEMO_INCLIST)
	$(CC) -c $< $(INCLUDE) -o $@

$(BINPATH)/test/%: $(OBJPATH)/test/%.o $(CALIB)
	$(CC) $< $(LINKOPT) -o $@

$(BINPATH)/demo/%: $(OBJPATH)/demo/%.o $(CALIB)
	$(CC) $< $(LINKOPT) -o $@

$(CALIB): $(CALIB_OBJLIST)
	$(AR) $(CALIB) $(CALIB_OBJLIST)

all: $(CALIB) $(TEST) $(DEMO) $(OBJLIST)

lib: $(CALIB)

delete: clean
	$(RM) $(CALIB) $(TEST) $(DEMO)
clean:
	$(RM) $(OBJLIST)

.PHONY: all lib clean delete
