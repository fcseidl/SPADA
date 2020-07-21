CXX = g++
CXXFLAGS = -std=c++1z -Wconversion -Wall -Werror -Wextra -pedantic

preprocess: driver.cpp csvtools.cpp csvtools.hpp
	$(CXX) $(CXXFLAGS) -c csvtools.cpp
	$(CXX) $(CXXFLAGS) -c driver.cpp
	$(CXX) $(CXXFLAGS) -o preprocess.exe driver.o csvtools.o
	
test: preprocess
	./preprocess.exe -ssf test_names.txt test_authors.txt ssf_names.txt ssf_authors.txt
	./preprocess.exe -corners ssf_names.txt
	./preprocess.exe -corners ssf_authors.txt

clean:
	rm *.o
	rm *.exe
	rm ssf*
