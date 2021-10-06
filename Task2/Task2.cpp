#include <omp.h>
#include <iostream>
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
	std::uniform_int_distribution<T> dis(0,10);
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

T dot_product_omp_reduction(const T* Vec1, const T* Vec2, const int n) {
	T Val = T(0);

#pragma omp parallel for reduction(+:Val)
	for (int i = 0; i < n; ++i) {
		T tmp = Vec1[i] * Vec2[i];
		Val += tmp;
	}
	return Val;
}


T dot_product_omp_divide(const T* Vec1, const T* Vec2, const int n, int cores = 4) {
	if (n <= cores) cores = 1;

	T* A = new T[cores];

#pragma omp parallel for schedule(static,1)
	for (int c = 0; c < cores; ++c) {
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

T dot_product_omp_atomic(const T* Vec1, const T* Vec2, const int n) {
	T Val = T(0);

#pragma omp parallel for
	for (int i = 0; i < n; ++i) {
		T tmp = Vec1[i] * Vec2[i];
		#pragma omp atomic
		Val += tmp;
	}
	return Val;
}

T dot_product_omp_critical(const T* Vec1, const T* Vec2, const int n) {
	T Val = T(0);

#pragma omp parallel for
	for (int i = 0; i < n; ++i) {
		T tmp = Vec1[i] * Vec2[i];
#pragma omp critical
		Val += tmp;
	}
	return Val;
}

T dot_product_omp_lock(const T* Vec1, const T* Vec2, const int n) {
	T Val = T(0);
	omp_lock_t lock;
	omp_init_lock(&lock);
#pragma omp parallel for
	for (int i = 0; i < n; ++i) {
		T tmp = Vec1[i] * Vec2[i];
		omp_set_lock(&lock);
		Val += tmp;
		omp_unset_lock(&lock);
	}
	omp_destroy_lock(&lock);
	return Val;
}

#define eps T(0.00001)
#define REPEATS 10

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
		std::cout << "vector length: " << N << std::endl;
		std::cout << "number of omp threads: " << cores << std::endl;
	}
	T* Vec1 = new T[N];
	T* Vec2 = new T[N];
	generate_random(Vec1, N);
	generate_random(Vec2, N);
	if (N <= 0 || cores <= 0 || Vec1 == nullptr || Vec2 == nullptr) throw std::overflow_error("error");
	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	T DP0, DP, DP1, DP2, DP3, DP4;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP0 = dot_product(Vec1, Vec2, N);	
	end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP = dot_product_omp_reduction(Vec1, Vec2, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP1 = dot_product_omp_divide(Vec1, Vec2, N, cores);	
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP2 = dot_product_omp_atomic(Vec1, Vec2, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff3 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP3 = dot_product_omp_critical(Vec1, Vec2, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff4 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP4 = dot_product_omp_lock(Vec1, Vec2, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff5 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	if (!silent) {
		std::clog << "time(us): \t\t\t" << diff << std::endl;
		std::clog << "time(us) omp reduction: \t" << diff1 << std::endl;
		std::clog << "time(us) omp divide: \t\t" << diff2 << std::endl;
		std::clog << "time(us) omp atomic: \t\t" << diff3 << std::endl;
		std::clog << "time(us) omp critical: \t\t" << diff4 << std::endl;
		std::clog << "time(us) omp lock: \t\t" << diff5 << std::endl;

		if (std::abs(DP0 - DP) <= eps && std::abs(DP0 - DP1) <= eps && std::abs(DP0 - DP2) <= eps && std::abs(DP0 - DP3) <= eps && std::abs(DP0 - DP4) <= eps)
			std::cout << "Dot product found OK: " << DP << std::endl;
		else
			std::cout << "Error: " << DP << " , " << DP1 << " , " << DP2 << " , " << DP3 << " , " << DP4 << "; Should be: " << DP0 << std::endl;
	}
	else {
		if (std::abs(DP0 - DP) <= eps && std::abs(DP0 - DP1) <= eps && std::abs(DP0 - DP2) <= eps && std::abs(DP0 - DP3) <= eps && std::abs(DP0 - DP4) <= eps) {
			std::cout << N << " " << cores << " ";
			std::cout << diff << " " << diff1 << " " << diff2 << " " << diff3 << " " << diff4 << " " << diff5 << std::endl;
		}
		else
			std::cout << 0 << " " << 0 << std::endl;
	}
	delete[] Vec1;
	delete[] Vec2;
	return 0;
}