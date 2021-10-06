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
	std::uniform_int_distribution<T> dis(0);
#endif

	//Data = new T[n];
	for (int i = 0; i < n; ++i) {
		Data[i] = dis(gen);
	}
}

/*
* reduction over min operator not supported in visual studio, because version of openmp 2.0 
T min(const T* Data, const int n) {
	T Val = Data[0];

#pragma omp parallel for schedule(static) firstprivate(Data,n) default(none) reduction(min:Val)
	for (int i = 0; i < n; ++i) {
		if (Val > Data[i]) Val = Data[i];
	}
}
*/

T min(const T* Data, const int n) {
	T Val = Data[0];

	for (int i = 0; i < n; ++i) {
		if (Val > Data[i]) Val = Data[i];
	}
	return Val;
}

T max(const T* Data, const int n) {
	T Val = Data[0];

	for (int i = 0; i < n; ++i) {
		if (Val < Data[i]) Val = Data[i];
	}
	return Val;
}

T min_omp_critical(const T* Data, const int n, int cores = 4) {
	if (n <= cores) cores = 1;
	
	T min = Data[0];

//#pragma omp parallel for schedule(static,1) firstprivate(Data,n,A,cores) default(none)
#pragma omp parallel for shared(min) schedule(static,1)
	for (int c = 0; c < cores; ++c) {
		//std::cout << "hello from " << omp_get_thread_num()<<"/" << omp_get_num_threads() << std::endl;
		int b = c * (n / cores);
		int e = (c + 1) * (n / cores);
		if (c == cores - 1) e = n;
		T tmp = Data[b];
		for (int i = b; i < e; ++i) {
			if (tmp > Data[i]) tmp = Data[i];
		}
		#pragma omp critical
		if (min > tmp) min = tmp;
	}
	return min;
}

T max_omp_critical(const T* Data, const int n, int cores = 4) {
	if (n <= cores) cores = 1;

	T max = Data[0];

//#pragma omp parallel for schedule(static,1) firstprivate(Data,n,A,cores) default(none)
#pragma omp parallel for shared(max) schedule(static,1)
	for (int c = 0; c < cores; ++c) {
		//std::cout << "hello from " << omp_get_thread_num() << std::endl;
		int b = c * (n / cores);
		int e = (c + 1) * (n / cores);
		if (c == cores - 1) e = n;
		T tmp = Data[b];
		for (int i = b; i < e; ++i) {
			if (tmp < Data[i]) tmp = Data[i];
		}
		#pragma omp critical
		if (max < tmp) max = tmp;
	}
	return max;
}

T min_omp_array(const T* Data, const int n, int cores = 4) {
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

T max_omp_array(const T* Data, const int n, int cores = 4) {
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
			if (tmp < Data[i]) tmp = Data[i];
		}
		A[c] = tmp;
	}

	T Val = A[0];

	for (int i = 0; i < cores; ++i) {
		if (Val < A[i]) Val = A[i];
	}

	delete[] A;
	return Val;
}

#define REPEATS 100

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

	T* Data = new T[N];
	generate_random(Data, N);
	if (N <= 0 || cores <= 0 || Data == nullptr) throw std::overflow_error("error");

	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	//auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	//auto diff1 = diff;
	T Extr0_1, Extr0_2, Extr_1, Extr_2;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i) {
		Extr0_1 = min(Data, N);
		Extr0_2 = max(Data, N);
	}
	end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i) {
		Extr_1 = min_omp_critical(Data, N, cores);
		Extr_2 = max_omp_critical(Data, N, cores);
	}
	end = std::chrono::high_resolution_clock::now();
	auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i) {
		Extr_1 = min_omp_array(Data, N, cores);
		Extr_2 = max_omp_array(Data, N, cores);
	}
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	/*
	start = std::chrono::high_resolution_clock::now();
	T Extr_3 = min_omp1(Data, N, cores);
	T Extr_4 = max_omp1(Data, N, cores);
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	*/
	if (!silent) {
		std::clog << "time(us): \t\t" << diff << std::endl;
		std::clog << "time(us) omp critical: \t" << diff1 << std::endl;
		std::clog << "time(us) omp array: \t" << diff2 << std::endl;
	
		if (Extr0_1 == Extr_1 && Extr0_2 == Extr_2)// && Extr0_1 == Extr_3 && Extr0_2 == Extr_4)
			std::cout << "Extremums found OK: " << Extr_1 << " , " << Extr_2 << std::endl;
		else
			std::cout << "Error: " << Extr_1 << " , " << Extr_2 /* << Extr_3 << " , " << Extr_4 */<< "; Should be: " << Extr0_1 << " , " << Extr0_2 << std::endl;
	}
	else {
		if (Extr0_1 == Extr_1 && Extr0_2 == Extr_2){// && Extr0_1 == Extr_3 && Extr0_2 == Extr_4) {
			std::cout << N << " " << cores << " ";
			std::cout << diff << " " << diff1 << " "  << diff2 << std::endl;
		}
		else 
			std::cout << 0 << " " << 0 << std::endl;
	}
	delete[] Data;
	return 0;
}