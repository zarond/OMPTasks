#include <omp.h>
#include <iostream>
#include <random>
#include <chrono>
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

T maxmin(const T* Mat, const int n, const int m) {
	T max = MINLIMIT;
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * m];
		for (int j = 0; j < m; ++j) {
			T v = Mat[i * m + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	return max;
}

T maxmin_omp(const T* Mat, const int n, const int m) {
	T max = MINLIMIT;
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		//std::cout << "hello from " << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
		T min = Mat[i * m];
		for (int j = 0; j < m; ++j) {
			T v = Mat[i * m + j];
			if (min > v)
				min = v;
		}
#pragma omp critical
		if (max < min) max = min;
	}
	return max;
}

// Experimental function. Wins at low N, because critical section used only NUM_THREADS times, not N. At N>5000 gives no bost over maxmin_omp().
// Unusable in the production because it requires static variables to run.
T maxmin_omp_threadprivate(const T* Mat, const int n, const int m, int cores) {
	static T max = MINLIMIT;
	#pragma omp threadprivate(max)

	#pragma omp parallel for copyin(max) schedule(static)
	for (int i = 0; i < n; ++i) {
		//std::cout << "hello from " << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
		T min = Mat[i * m];
		for (int j = 0; j < m; ++j) {
			T v = Mat[i * m + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	T Amax = MINLIMIT;
	#pragma omp parallel for shared(Amax) schedule(static,1)
	for (int i = 0; i < cores; ++i){
		//std::cout << "hello from " << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
		#pragma omp critical
		if (Amax < max) Amax = max;
	}
	return Amax;
}

T min_omp(const T* Data, const int n, int cores = 4) {
	if (n <= cores) cores = 1;

	T* A = new T[cores];

	//#pragma omp parallel for schedule(static,1) firstprivate(Data,n,A,cores) default(none)
#pragma omp parallel for schedule(static,1)
	for (int c = 0; c < cores; ++c) {
		//std::cout << "hello from " << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
		int b = c * (n / cores);
		int e = (c + 1) * (n / cores);
		if (c == cores - 1) e = n;
		T tmp = Data[b];
		for (int i = b; i < e; ++i) {
			if (tmp > Data[i]) tmp = Data[i];
		}
		A[c] = tmp;
	}

	T Val = A[0];

	for (int i = 0; i < cores; ++i) {
		if (Val > A[i]) Val = A[i];
	}

	delete[] A;
	return Val;
}

T maxmin_omp_nested(const T* Mat, const int n, const int m, const int cores = 4) {
	T max = MINLIMIT;
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		//std::cout << "hello from " << omp_get_thread_num() << "/" << omp_get_num_threads() << std::endl;
		T min = min_omp( &Mat[i * m], n, cores);
		#pragma omp critical
		if (max < min) max = min;
	}
	return max;
}

#define eps T(0.00001)
#define REPEATS 10

int main(int argc, char** argv) {
	int N = 5000;
	int M = N;
	int cores = omp_get_num_procs();
	bool silent = false;
	if (argc >= 2) {
		N = std::atoi(argv[1]);
		M = N;
	}
	//if (argc >= 3) {
	//	M = std::atoi(argv[2]);
	//}
	if (argc >= 3) {
		cores = std::atoi(argv[2]);
	}
	if (argc >= 4) {
		silent = true;
	}
	omp_set_num_threads(cores);
	if (!silent) {
		std::cout << "matrix N, M: " << N << ", " << M << std::endl;
		std::cout << "number of omp threads: " << cores << std::endl;
	}
	T* Mat = new T[N * M];
	generate_random(Mat, N * M);

	if (N <= 0 || M <= 0 || cores <= 0 || Mat == nullptr) throw std::overflow_error("error");
	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	T DP0, DP, DP1;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i) 
		DP0 = maxmin(Mat, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP = maxmin_omp(Mat, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	omp_set_nested(1);
	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP1= maxmin_omp_nested(Mat, N, M, cores);
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;
	omp_set_nested(0);
	
	/*
	T DP2;
	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		DP2 = maxmin_omp_threadprivate(Mat, N, M, cores);
	end = std::chrono::high_resolution_clock::now();
	auto diff3 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;
	*/
	/*
	start = std::chrono::high_resolution_clock::now();
	T DP1 = dot_product_omp1(Vec1, Vec2, N, cores);
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	*/
	if (!silent) {
		std::clog << "time(us): \t\t" << diff << std::endl;
		std::clog << "time(us) omp: \t\t" << diff1 << std::endl;
		std::clog << "time(us) omp nested: \t" << diff2 << std::endl;
		//std::clog << "time(us) omp threadprivate: \t" << diff3 << std::endl;

		if (std::abs(DP0 - DP) <= eps && std::abs(DP0 - DP1) <= eps)//&& std::abs(DP0 - DP2) <= eps)
			std::cout << "maxmin found OK: " << DP << std::endl;
		else
			std::cout << "Error: " << DP << "; Should be: " << DP0 << std::endl;
	}
	else {
		if (std::abs(DP0 - DP) <= eps && std::abs(DP0 - DP1) <= eps) {
			std::cout << N << " " << cores << " ";
			std::cout << diff << " " << diff1 << " " << diff2 << std::endl;
		}
		else
			std::cout << 0 << " " << 0 << std::endl;
	}
	delete[] Mat;
	return 0;
}