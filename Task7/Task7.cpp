#include <omp.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

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
	std::uniform_real_distribution<T> dis(T(-1), T(1));
#else
	std::uniform_int_distribution<T> dis(0, 10);
#endif

	//Data = new T[n];
	for (int i = 0; i < n; ++i) {
		Data[i] = dis(gen);
	}
}

T dot_product(const T* Vec1, const T* Vec2, const int n) {
	T Val = T(0);

	for (int i = 0; i < n; ++i) {
		Val += Vec1[i] * Vec2[i];
	}
	return Val;
}

/*
T dot_product_omp(const T* Vec1, const T* Vec2, const int n) {
	T Val = T(0);

#pragma omp parallel for shared(Vec1, Vec2) reduction(+:Val)
	for (int i = 0; i < n; ++i) {
		//std::cout << "hello from " << omp_get_thread_num() << std::endl;
		T tmp = Vec1[i] * Vec2[i];
		Val += tmp;
	}
	return Val;
}
*/

T dot_product_omp(const T* Vec1, const T* Vec2, const int n, int cores = 4) {
	if (n <= cores) cores = 1;

	T* A = new T[cores];

#pragma omp parallel for schedule(static,1)
	for (int c = 0; c < cores; ++c) {
		//std::cout << "hello from " << omp_get_thread_num() << std::endl;
		int b = c * (n / cores);
		int e = (c + 1) * (n / cores);
		if (c == cores - 1) e = n;
		T tmp = T(0);
		for (int i = b; i < e; ++i) {
			tmp += Vec1[i] * Vec2[i];
		}
		A[c] = tmp;
	}

	T Val = T(0);

	for (int i = 0; i < cores; ++i) {
		Val += A[i];
	}

	delete[] A;
	return Val;
}

void load_vec_txt(std::ifstream &file, T* vec, const int N) {
	T v = T(0);
	int i = 0;
	while (i < N){
		for (; i < N && (file >> v); ++i) {
			vec[i] = v;
		}
		if (file.eof()) { file.clear(); file.seekg(0, std::ios::beg); /*std::cout << "go to begining";*/ }
	}
}

void load_vec_binary(std::ifstream& file, T* vec, const int N) {
	//T v = T(0);
	int i = 0;
	int ToRead = N;
	if (file.eof()) { file.clear(); file.seekg(0, std::ios::beg); /*std::cout << "go to begining";*/ }
	while (i < N) {
		file.read(reinterpret_cast<char*>(&vec[i]), ToRead *sizeof(T));
		if (file) {
			//All read sucessfuly
			i += N;
			ToRead -= N;
		}
		else {
			int read = file.gcount();
			i += read;
			ToRead -= read;
		}
		if (file.eof()) { file.clear(); file.seekg(0, std::ios::beg); /*std::cout << "go to begining";*/ }
	}
}

#define eps T(0.0001)
//#define REPEATS 10

int main(int argc, char** argv) {
	int N = 100;
	int vector_pairs = 16;
	int cores = omp_get_num_procs();
	bool silent = false;
	if (argc >= 2) {
		N = std::atoi(argv[1]);
	}
	if (argc >= 3) {
		cores = std::atoi(argv[2]);
	}
	if (argc >= 4) {
		silent = true;
	}
	omp_set_num_threads(cores);
	int vec_read = 0;
	int vec_computed = 0;

	//std::ifstream file("in.txt");
	std::ifstream file("in.b");
	if (!file) { file.close();  std::cout << "no file"; return 1; }

	if (!silent) {
		std::cout << "vector length: " << N << std::endl;
		std::cout << "number of omp threads: " << cores << std::endl;
	}
	T* Vec1 = new T[N];
	T* Vec2 = new T[N];
	T* VecR1 = new T[N];
	T* VecR2 = new T[N];
	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	auto diff1 = diff;
	auto startL = start;
	auto endL = start;
	auto diffL = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	T DP0, DP;

	load_vec_binary(file, VecR1, N);
	load_vec_binary(file, VecR2, N);

	omp_set_nested(1);

	for (int i = 0; i < vector_pairs; ++i) {
	#pragma omp parallel num_threads(2)
		{
		#pragma omp master //shared(Vec1,Vec2,VecR1,VecR2)
		{
			std::swap(Vec1, VecR1);
			std::swap(Vec2, VecR2);
		}
		#pragma omp sections
		{
			#pragma omp section
			{
				startL = std::chrono::high_resolution_clock::now();
				//std::cout << "load file from " << omp_get_thread_num() << std::endl;
				load_vec_binary(file, VecR1, N);
				load_vec_binary(file, VecR2, N);
				endL = std::chrono::high_resolution_clock::now();
				diffL += std::chrono::duration_cast<std::chrono::microseconds>(endL - startL).count();
			}
			#pragma omp section
			{
				//std::cout << "work from " << omp_get_thread_num() << std::endl;
				start = std::chrono::high_resolution_clock::now();
				DP0 = dot_product(Vec1, Vec2, N);
				end = std::chrono::high_resolution_clock::now();
				diff += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

				start = std::chrono::high_resolution_clock::now();
				DP = dot_product_omp(Vec1, Vec2, N);
				end = std::chrono::high_resolution_clock::now();
				diff1 += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
			}
		}
		#pragma omp barrier
		}
	}
	diffL /= vector_pairs;
	diff /= vector_pairs;
	diff1 /= vector_pairs;

	if (!silent) {
		std::clog << "time(us): \t\t" << diff << std::endl;
		std::clog << "time(us) omp: \t" << diff1 << std::endl;
		std::clog << "time(us) load: \t" << diffL << std::endl;
		//std::clog << "time(us) omp1: \t" << diff2 << std::endl;

		if (std::abs(DP0 - DP) <= eps)// && std::abs(DP0 - DP1) <= eps)
			std::cout << "Dot product found OK: " << DP << std::endl;
		else
			std::cout << "Error: " << DP /* << " , " << DP1*/ << "; Should be: " << DP0 << std::endl;
	}
	else {
		if (std::abs(DP0 - DP) <= eps)// && std::abs(DP0 - DP1) <= eps)
		{
			std::cout << N << " " << cores << " ";
			std::cout << diff << " " << diff1 << " " << diffL << std::endl;
		}
		else
			std::cout << 0 << " " << 0 << std::endl;
	}
	delete[] Vec1;
	delete[] Vec2;
	delete[] VecR1;
	delete[] VecR2;
	file.close();
	return 0;
}