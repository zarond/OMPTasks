#include <omp.h>
#include <iostream>
#include <random>
#include <chrono>
#include <time.h>
#include< limits >

typedef double T;

T function(T x) {
	//T Val = x * x + std::cos(x);
	T Val = std::exp(x);
	return Val;
}

T integral(T(*func)(T), const T a, const T b, const int n) {
	T h = (b - a) / n;
	T Val = T(0);

	for (int i = 0; i < n; ++i) {
		T x = a + i * h;
		Val += func(x) * h;
	}
	return Val;
}

T integral_omp_reduction(T(*func)(T), const T a, const T b, const int n) {
	T h = (b - a) / n;
	T Val = T(0);

#pragma omp parallel for reduction(+:Val)
	for (int i = 0; i < n; ++i) {
		T x = a + i * h;
		Val += func(x) * h;
	}
	return Val;
}

T integral_omp_divide(T(*func)(T), const T a, const T b, const int n, int cores = 4) {
	if (n <= cores) cores = 1;
	T h = (b - a) / n;
	T Val = T(0);

	T* A = new T[cores];

#pragma omp parallel for schedule(static,1)
	for (int c = 0; c < cores; ++c) {
		int b = c * (n / cores);
		int e = (c + 1) * (n / cores);
		if (c == cores - 1) e = n;
		T tmp = T(0);
		for (int i = b; i < e; ++i) {
			T x = a + i * h;
			tmp += func(x) * h;
		}
		A[c] = tmp;
	}

	for (int i = 0; i < cores; ++i) {
		Val += A[i];
	}

	delete[] A;
	return Val;
}

#define eps T(0.0001)
#define REPEATS 10
#define MIN_EPS std::numeric_limits<T>::min()

int main(int argc, char** argv) {
	int N = 1000000;
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
	if (!silent) {
		std::cout << "N: " << N << std::endl;
		std::cout << "number of omp threads: " << cores << std::endl;
	}

	T a = T(0);
	T b = T(10);

	if (N <= 0 || cores <= 0) throw std::overflow_error("error");

	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	T I0, I, I1;

	//T sum = 0;
	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i){
		I0 = integral(function, a + MIN_EPS, b, N);
		//sum += I0 - i;
	}
	end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		I = integral_omp_reduction(function, a + MIN_EPS, b, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		I1 = integral_omp_divide(function, a + MIN_EPS, b, N, cores);
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	if (!silent) {
		std::clog << "time(us): \t\t" << diff << std::endl;
		std::clog << "time(us) omp reduction: " << diff1 << std::endl;
		std::clog << "time(us) omp divide: \t" << diff2 << std::endl;

		if (std::abs(I0 - I) <= eps )//&& std::abs(I0 - (std::exp(b) - std::exp(a))) <= eps)
			std::cout << "integral found OK: " << I << /*"ignore:" /<< sum <<*/ std::endl;
		else
			std::cout << "Error: " << I << "; Should be: " << I0 << /*"ignore:" << sum <<*/ std::endl;
	}
	else {
		if (std::abs(I0 - I) <= eps) {
			std::cout << N << " " << cores << " ";
			std::cout << diff << " " << diff1 << " " << diff2 << std::endl;
		}
		else
			std::cout << 0 << " " << 0 << std::endl;
	}
	return 0;
}