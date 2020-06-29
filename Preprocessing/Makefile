process: Preprocessing.cpp
	g++ -std=c++11 Preprocessing.cpp -o process.exe
	
test: process
	./process.exe test_names.txt test_authors.txt ssf_names.txt ssf_authors.txt
	cat ssf_names.txt
	cat ssf_authors.txt

clean:
	rm ssf*
	rm *.exe
