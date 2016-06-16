CXX = g++
TARGET = clustering_nonparametric_bayes
CXXFLAGS = -Wall -O3 --std=c++11 -I/usr/local/include/eigen3/
LDFLAGS = -lm 
SRCS := $(wildcard *.cc)
OBJS := $(SRCS:.cc=.o)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) 

clean:
	rm -f $(TARGET) $(OBJS)
