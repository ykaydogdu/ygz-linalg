CC = g++
CFLAGS = -std=c++11
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include
TESTDIR = ./test
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# link object files
libygzlinalg.a: $(OBJS)
	ar rcs $@ $^

# compile the source files into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -I $(INCDIR) -c $< -o $@

test: libygzlinalg.a test.o
	$(CC) $(CFLAGS) -I $(INCDIR) -L. -o $(TESTDIR)/test $(TESTDIR)/test.o -lygzlinalg
	$(TESTDIR)/test

test.o:
	$(CC) $(CFLAGS) -I $(INCDIR) -c $(TESTDIR)/test.cpp -o $(TESTDIR)/test.o