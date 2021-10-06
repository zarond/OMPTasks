#include <omp.h>
#include <iostream>
#include <random>
#include <chrono>
//#include <ctime>
#include <limits>
#include <cstddef>

#define REAL

#ifdef REAL
typedef double T;
#else
typedef int T;
#endif

#define MINLIMIT -std::numeric_limits<T>::max()

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

T work(const T* Mat, const int n, const int m) {
	T min = Mat[0];
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			T v = Mat[i * m + j];
			if (min > v)
				min = v;
		}
	}
	return min;
}

T work_omp(const T* Mat, const int n, const int m) {
	T min = Mat[0];
#pragma omp parallel for schedule(static,1) firstprivate(min) lastprivate(min)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			T v = Mat[i * m + j];
			if (min > v)
				min = v;
		}
	}
	return min;
}

void array_delay(int delaylength, double a[1]) {

	int i;
	a[0] = 1.0;
	for (i = 0; i < delaylength; i++)
		a[0] += i;
	if (a[0] < 0)
		printf("%f \n", a[0]);

}

double* light(const long long n) {
	double a[1];
	{
		array_delay(n, a);
	}
	return a;
}

double* light_omp(const long long n) {
	double a[1];
#pragma omp parallel private(a) //numthreads()
	{
		array_delay(n, a);
	}
	return a;
}


#define eps T(0.00001)
#define REPEATS 1000

int main(int argc, char** argv) {
	int N = 2000;
	//int Ns[] = { 16, 5, 10, 20, 30, 40, 50, 60 };
	//int Ms[] = { 1, 10, 100, 200, 500, 1000, 2000, 5000 };
	int M = 2000;
	int cores = omp_get_num_procs();
	bool silent = false;
	/*
	if (argc >= 2) {
		N = std::atoi(argv[1]);
		M = N;
	}
	if (argc >= 3) {
		cores = std::atoi(argv[2]);
	}
	if (argc >= 4) {
		silent = true;
	}*/
	/*
	omp_set_num_threads(cores);
	if (!silent) {
		std::cout << "matrix N, M: " << N << ", " << M << std::endl;
		std::cout << "number of omp threads: " << cores << std::endl;
	}
	*/
	T* Mat = new T[N * M * 16];
	generate_random(Mat, N * M);

	if (N <= 0 || M <= 0 || cores <= 0 || Mat == nullptr) throw std::overflow_error("error");
	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	T DP0, DP1, DP2;

	std::cout << "----------------------" << "bench work 1" << std::endl;
	/*
	for (int j = 0; j < 8; ++j) {
		int n0 = N;// Ns[j];
		int m = Ms[j];

		for (int k = 1; k <= 16; ++k) {
			omp_set_num_threads(k);
			
			int n = n0 * k;

			start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < REPEATS; ++i)
				DP0 = work(Mat, n0, m);
			end = std::chrono::high_resolution_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

			start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < REPEATS; ++i)
				DP1 = work_omp(Mat, n0, m);
			end = std::chrono::high_resolution_clock::now();
			auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

			start = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < REPEATS; ++i)
				DP2 = work_omp(Mat, n, m);
			end = std::chrono::high_resolution_clock::now();
			auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

			std::cout << n0 << " " << m << " " << k << " ";
			std::cout << diff << " " << diff1 << " " << diff2 <<  std::endl;
		}

	}*/
	std::cout << "----------------------" << "bench work 2" << std::endl;

	int n0 = 1;// Ns[j];
	int m = 1;// Ms[j];

	for (int k = 1; k <= 32; ++k) {
		omp_set_num_threads(k);

		int n = n0 * k;

		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < REPEATS; ++i)
			DP0 = work(Mat, n0, m);
		end = std::chrono::high_resolution_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / REPEATS;

		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < REPEATS; ++i)
			DP2 = work_omp(Mat, n, m);
		end = std::chrono::high_resolution_clock::now();
		auto diff1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / REPEATS;

		std::cout << n0 << " " << m << " " << k << " ";
		std::cout << diff << " " << diff1 << std::endl;
	}

	std::cout << "----------------------" << "light bench" << std::endl;

	for (int z = 12; z <= 20; z+=2)
	for (int k = 1; k <= 32; ++k) {
		omp_set_num_threads(k);

		long long n = 1; n <<= z;

		//std::clock_t c_start = std::clock();
		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < REPEATS; ++i)
			DP0 = *light(n);
		//std::clock_t c_end = std::clock();
		end = std::chrono::high_resolution_clock::now();
		auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / REPEATS;
		//auto c_diff = (1000000000.0 / REPEATS) * (c_end - c_start) / (CLOCKS_PER_SEC);// (c_end - c_start) / REPEATS;

		//c_start = std::clock();
		start = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < REPEATS; ++i)
			DP2 = *light_omp(n);
		//c_end = std::clock();
		end = std::chrono::high_resolution_clock::now();
		auto diff1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / REPEATS;
		//auto c_diff1 = (1000000000.0 / REPEATS) * (c_end - c_start) / (CLOCKS_PER_SEC);

		std::cout << n << " " << k << " ";
		std::cout << diff << " " << diff1 << std::endl;
		//std::cout << " " << c_diff << " " << c_diff1 << std::endl;
	}


	/*
	if (!silent) {
		std::clog << "time(us): \t\t" << diff << std::endl;
		std::clog << "time(us) omp: \t\t" << diff1 << std::endl;
		//std::clog << "time(us) omp threadprivate: \t" << diff3 << std::endl;

	}
	else {
		std::cout << N << " " << cores << " ";
		std::cout << diff << " " << diff1 << std::endl;
	}
	*/
	delete[] Mat;
	return 0;
}