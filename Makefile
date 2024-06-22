CC = g++
CFLAGS = -std=c++11
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# link object files
linalg.a: $(OBJS)
	ar rcs $@ $^

# compile the source files into object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -I $(INCDIR) -c $< -o $@