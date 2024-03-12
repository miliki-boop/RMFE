CXX = g++
CXXFLAGS = -g -O2 -std=c++11 -pthread -march=native -DNTL_EXCEPTIONS
LDFLAGS = -lntl -lgf2x -lgmp -lm

TARGET = RMFE_GR
SRCS = util.cpp test_RMFE_GR_com.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET)