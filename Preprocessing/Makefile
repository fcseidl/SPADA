CXX = g++
CXXFLAGS = -std=c++1z -Wconversion -Wall -Werror -Wextra -pedantic

preprocess: driver.cpp Preprocessing.cpp Preprocessing.hpp
	$(CXX) $(CXXFLAGS) -c Preprocessing.cpp
	$(CXX) $(CXXFLAGS) -c driver.cpp
	$(CXX) $(CXXFLAGS) -o preprocess.exe driver.o Preprocessing.o
	
test: preprocess
	./preprocess.exe -ssf test_names.txt test_authors.txt ssf_names.txt ssf_authors.txt
	./preprocess.exe -corners ssf_names.txt
	./preprocess.exe -corners ssf_authors.txt

clean:
	rm *.exe
	rm ssf*
