CXX = g++
CXXFLAGS = -std=c++1z -Wconversion -Wall -Werror -Wextra -pedantic

preprocess: driver.cpp csvtools.cpp csvtools.hpp topics.cpp topics.hpp
	$(CXX) $(CXXFLAGS) -c csvtools.cpp
	$(CXX) $(CXXFLAGS) -c driver.cpp
	$(CXX) $(CXXFLAGS) -c topics.cpp
	$(CXX) $(CXXFLAGS) -o preprocess.exe driver.o csvtools.o topics.o
	rm *.o
	
test: preprocess
	./preprocess.exe -ssf test_names.csv test_authors.csv ssf_names.csv ssf_authors.csv
	./preprocess.exe -corners ssf_names.csv
	./preprocess.exe -corners ssf_authors.csv
	./preprocess.exe -join test_names.csv test_authors.csv test_names_authors_joint.csv
	./preprocess.exe -corners test_names_authors_joint.csv -ignore 1
	./preprocess.exe -wordbag abstract.txt bags_abstract.csv 5
	./preprocess.exe -corners bags_abstract.csv

clean:
	rm *.exe
	rm ssf*
	rm *joint*
	rm bags*
