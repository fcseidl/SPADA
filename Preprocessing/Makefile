preprocess: Preprocessing.cpp
	g++ -std=c++11 Preprocessing.cpp -o preprocess.exe
	
test: preprocess
	./preprocess.exe -ssf test_names.txt test_authors.txt ssf_names.txt ssf_authors.txt
	./preprocess.exe -corners ssf_names.txt
	./preprocess.exe -corners ssf_authors.txt

clean:
	rm ssf*
	rm *.exe
