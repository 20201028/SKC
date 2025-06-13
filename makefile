CC = g++
CXXFLAGS += -std=c++11

main:./*.cpp
	$(CXX) $(CXXFLAGS) -O2 -o main ./*.cpp

# clean:
# 	rm -rf main