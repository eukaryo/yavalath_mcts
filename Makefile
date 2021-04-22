yavalath_cli: Source.cpp
	g++ -std=c++1z -static-libstdc++ -O2 -flto -march=native -DNDEBUG Source.cpp -o yavalath_cli
clean:
	rm -f *.o yavalath_cli