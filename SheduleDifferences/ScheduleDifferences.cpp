#include <omp.h>
#include <iostream>
#include <random>
#include <chrono>


/* Iterative Function to calculate (x^y)%p in O(log y) */
int power(long long x, unsigned int y, int p)
{
	int res = 1;     // Initialize result

	x = x % p; // Update x if it is more than or
				// equal to p

	if (x == 0) return 0; // In case x is divisible by p;

	while (y > 0)
	{
		// If y is odd, multiply x with result
		if (y & 1)
			res = (res * x) % p;

		// y must be even now
		y = y >> 1; // y = y/2
		x = (x * x) % p;
	}
	return res;
}

void workload_uneven(int* Data, unsigned int n) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> dis(0,1000);

	for (int i = 0; i < n; ++i) {
		int x = i + 10;
		long long val = power(long long (x), 65537,  x / 2);
		if (val % 20 == 0) {
			Data[i] = val;
			for (int j = 0; j < 10; ++j)
				Data[i] += dis(gen); 
		}
		else Data[i] = int(val);
	}
}

void workload_uneven_omp(int* Data, unsigned int n) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> dis(0, 1000);

#pragma omp parallel for schedule(runtime) firstprivate(gen,dis)
	for (int i = 0; i < n; ++i) {
		int x = i + 10;
		long long val = power(long long(x), 65537, x / 2);
		if (val % 20 == 0) {
			Data[i] = val;
			for (int j = 0; j < 10; ++j)
				Data[i] += dis(gen);
		}
		else Data[i] = int(val);
	}
}

void Collatz(int* Data, const int n) {
	for (int i = 0; i < n; ++i) {
		int steps = 0;
		long long val = i+1;
		while (val > 1) {
			val = (val % 2 == 0) ? val/2: 3*val+1;
			++steps;
		}
		Data[i] = steps;
	}
}

void Collatz_omp(int* Data, const int n) {
	#pragma omp parallel for schedule(runtime)
	for (int i = 0; i < n; ++i) {
		int steps = 0;
		long long val = i + 1;
		while (val > 1) {
			val = (val % 2 == 0) ? val / 2 : 3 * val + 1;
			++steps;
		}
		Data[i] = steps;
	}
}

#define eps T(0.00001)
#define REPEATS 10

int main(int argc, char** argv) {
	int N = 10000;
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
	int* Data = new int[N];
	int* DataInt = new int[N];

	//std::cout << getenv("OMP_SCHEDULE") << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	auto end = start;
	
	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		workload_uneven(Data, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		workload_uneven_omp(Data, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		Collatz(DataInt, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < REPEATS; ++i)
		Collatz_omp(DataInt, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff3 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / REPEATS;

	if (!silent) {
		std::clog << "schedule mode: " << getenv("OMP_SCHEDULE") << std::endl;
		std::clog << "time(us) uneven:\t" << diff << std::endl;
		std::clog << "time(us) uneven omp:\t" << diff1 << std::endl;
		std::clog << "time(us) Collatz:\t" << diff2 << std::endl;
		std::clog << "time(us) Collatz omp:\t" << diff3 << std::endl;
	}
	else {
		std::cout << N << " " << cores << " ";
		std::cout << diff << " " << diff1 << " " << diff2 << " " << diff3 << std::endl;
	}
	delete[] Data;
	delete[] DataInt;
	return 0;
}