CXX = g++
CXXFLAGS = -Wall

.PHONY : all
all : main

main : main.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY : clean
clean :
	rm main
