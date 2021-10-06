#include <omp.h>
#include <iostream>
#include <random>
#include <chrono>
#include <fstream>

#define REAL

#ifdef REAL
typedef double T;
#else
typedef int T;
#endif


void generate_random(T* Data, unsigned int n) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
#ifdef REAL
	std::uniform_real_distribution<T> dis(T(0), T(1));
#else
	std::uniform_int_distribution<T> dis(0);
#endif

	//Data = new T[n];
	for (int i = 0; i < n; ++i) {
		Data[i] = dis(gen);
	}
}


int main(int argc, char** argv) {
	int N = 100;
	int M = 16;
	int Mode = 0;

	if (argc >= 2) {
		N = std::atoi(argv[1]);
	}
	if (argc >= 3) {
		M = std::atoi(argv[2]);
	}
	if (argc >= 4) {
		Mode = std::atoi(argv[3]);
	}

	T** Data = new T*[M];
	for (int i = 0; i < M; ++i) {
		Data[i] = new T[N];
		generate_random(Data[i], N);
	}

	if (N <= 0 || M < 0 || Data == nullptr) throw std::overflow_error("error");

	//std::ofstream file("in.txt");
	if (Mode == 0 || Mode == 1){
		std::ofstream file("in.txt");
		for (int j = 0; j < M; ++j) {
			for (int i = 0; i < N; ++i) {
				file << Data[j][i] << " ";
			}
			file << std::endl;
		}
		file.close();
	}
	if (Mode != 0){
		std::ofstream file("in.b", std::ios::binary);
		for (int j = 0; j < M; ++j) {
			for (int i = 0; i < N; ++i) {
				T x = Data[j][i];
				file.write(reinterpret_cast<char*>(&x), sizeof(x));
			}
		}
		file.close();
	}
	//file.close();
	for (int i = 0; i < M; ++i) {
		delete[] Data[i];
	}
	delete[] Data;
	return 0;
}