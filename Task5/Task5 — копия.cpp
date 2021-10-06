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

#define MINLIMIT std::numeric_limits<T>::min()

void generate_random(T* Data, unsigned int n) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
#ifdef REAL
	std::uniform_real_distribution<T> dis(T(-1), T(1));
#else
	std::uniform_int_distribution<T> dis(0, 10);
#endif

	for (int i = 0; i < n; ++i) {
		Data[i] = dis(gen);
	}
}

void generate_offsets(int* Data, unsigned int n, unsigned int l) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> dis(1, n / 4);
	for (int i = 0; i < l; ++i) {
		Data[i] = dis(gen);
	}
}

void upper_triangle(T* Data, unsigned int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			Data[i * n + j] = T(0);
		}
	}
}

void lower_triangle(T* Data, unsigned int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = i; j < n; ++j) {
			Data[i * n + j] = T(0);
		}
	}
}

void band_matrix(T* Data, unsigned int n, int k1 = 1, int k2 = 3) {
	k1 = std::max(0, k1); k2 = std::max(0, k2);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i - k1; ++j) {
			Data[i * n + j] = T(0);
		}
		for (int j = i + k2 + 1; j < n; ++j) {
			Data[i * n + j] = T(0);
		}
	}
}

//???
/*
void block_matrix(T* Data, unsigned int n) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> dis(0, n/4);
	int c = 0, offset = 0, offsettotal = 0, i0 = dis(gen);
	for (int i = i0; i < n; ++i) {
		if (c == 0) { c = dis(gen); offset = dis(gen); offsettotal += offset;}
		for (int j = 0; j < n && j < offsettotal; ++j) {
			Data[i * n + j] = T(0);
		}
		--c;
	}
}*/

//???
/*
void block_diag_matrix(T* Data, unsigned int n) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<int> dis(0, n / 5);
	int c = 0, offset = 0, i0 = dis(gen);
	for (int i = i0; i < n; ++i) {
		if (c == 0) { c = dis(gen); offset = i; }
		for (int j = 0; j < n && j < offset; ++j) {
			Data[i * n + j] = T(0);
		}
		for (int j = 0; j < n && j < offset; ++j) {
			Data[j * n + i] = T(0);
		}
		--c;
	}
}
*/
void block_diag_matrix(T* Data, unsigned int n, unsigned  int l, int* offsets) {
	int c = offsets[0], offset = 0;
	int counter = 1;
	for (int i = 0; i < n; ++i) {
		if (c == 0 && counter < l) { c = offsets[counter++]; offset = i; }
		for (int j = 0; j < n && j < offset; ++j) {
			Data[i * n + j] = T(0);
		}
		for (int j = 0; j < n && j < offset; ++j) {
			Data[j * n + i] = T(0);
		}
		--c;
	}
}

#pragma optimize( "", off )
/* unoptimized code section */
T maxmin(const T* Mat, const int n, const int m) {
	T max = MINLIMIT;
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n];
		for (int j = 0; j < m; ++j) {
			T v = Mat[i * m + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	return max;
}

T maxmin_upper(const T* Mat, const int n) {
	T max = MINLIMIT;
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n];
		for (int j = i; j < n; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	return max;
}

T maxmin_lower(const T* Mat, const int n) {
	T max = MINLIMIT;
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n + n - 1];
		for (int j = 0; j <= i; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	return max;
}

T maxmin_band(const T* Mat, const int n, int k1 = 1, int k2 = 3) {
	k1 = std::max(0, k1); k2 = std::max(0, k2);
	T max = MINLIMIT;
	for (int i = 0; i < n; ++i) {
		T min = T(0);// Mat[i * n];
		int x1 = std::max(0, i - k1 + 1); int x2 = std::min(n, i + k2 + 1);
		for (int j = x1; j < x2; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	return max;
}

//???
/*
T maxmin_block(const T* Mat, const int n) {
	T max = MINLIMIT;
	T offset = 0;
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n];
		for (int j = offset; j < n && Mat[i * n + j] == T(0); ++j, ++offset) {}
		for (int j = offset; j < n; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
	}
	return max;
} */

//???
T maxmin_block_diag(const T* Mat, const int n, unsigned  int l, int* offsets) {
	T max = MINLIMIT;
	int c = offsets[0], offset = 0;// , i0 = offsets[0];
	int counter = 1;
	int x1 = 0, x2 = c;
	for (int i = 0; i < n; ++i) {
		if (c == 0 && counter < l) { c = offsets[counter++]; offset = i; x1 = x2; x2 += c; }
		T min = T(0);//Mat[i * n];
		for (int j = x1; j < x2 && j < n ; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
		if (max < min) max = min;
		--c;
	}
	return max;
}
#pragma optimize("", on)

T maxmin_omp(const T* Mat, const int n, const int m) {
	T max = MINLIMIT;
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n];
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

T maxmin_omp_upper(const T* Mat, const int n) {
	T max = MINLIMIT;
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n];
		for (int j = i; j < n; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
#pragma omp critical
		if (max < min) max = min;
	}
	return max;
}

T maxmin_omp_lower(const T* Mat, const int n) {
	T max = MINLIMIT;
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		T min = Mat[i * n + n - 1];
		for (int j = 0; j <= i; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
#pragma omp critical
		if (max < min) max = min;
	}
	return max;
}

T maxmin_omp_band(const T* Mat, const int n, int k1 = 1, int k2 = 3) {
	k1 = std::max(0, k1); k2 = std::max(0, k2);
	T max = MINLIMIT;
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		T min = T(0);// Mat[i * n];
		int x1 = std::max(0, i - k1 + 1); int x2 = std::min(n, i + k2 + 1);
		for (int j = x1; j < x2; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
#pragma omp critical
		if (max < min) max = min;
	}
	return max;
}

T maxmin_omp_block_diag(const T* Mat, const int n, unsigned  int l, int* offsets) { //???
	T max = MINLIMIT;
	T* lims = new T[2 * n];
	int c = offsets[0], offset = 0;// , i0 = offsets[0];
	int counter = 1;
	int x1 = 0, x2 = c;
	for (int i = 0; i < n; ++i, --c) {
		if (c == 0 && counter < l) { c = offsets[counter++]; offset = i; x1 = x2; x2 += c; }
		lims[2 * i] = x1;
		lims[2 * i + 1] = x2;
	}
#pragma omp parallel for shared(max)
	for (int i = 0; i < n; ++i) {
		int l1 = lims[2 * i], l2 = lims[2 * i + 1];
		T min = T(0);//Mat[i * n];
		for (int j = l1; j < l2 && j < n; ++j) {
			T v = Mat[i * n + j];
			if (min > v)
				min = v;
		}
#pragma omp critical
		if (max < min) max = min;
	}
	delete[] lims;
	return max;
}

#define eps T(0.00001)

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
	T* MatU = new T[N * M];
	generate_random(Mat, N * M);
	//std::memcpy(MatU, Mat, sizeof(T) * N * M);

	if (N <= 0 || M <= 0 || cores <= 0 || Mat == nullptr) throw std::overflow_error("error");

	auto start = std::chrono::high_resolution_clock::now();
	T DP0 = maxmin(Mat, N, M);
	auto end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP = maxmin_omp(Mat, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	// upper triangle
	std::memcpy(MatU, Mat, sizeof(T) * N * M);
	upper_triangle(MatU, N);

	start = std::chrono::high_resolution_clock::now();
	T DP0_upper = maxmin(MatU, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff_upper = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP_upper_omp = maxmin_omp(MatU, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff_upper_omp = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP_upper_omp_s = maxmin_omp_upper(MatU, N);
	end = std::chrono::high_resolution_clock::now();
	auto diff_upper_omp_s = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	//band
	std::memcpy(MatU, Mat, sizeof(T) * N * M);
	int k1 = 2, k2 = 6;
	band_matrix(MatU, N, k1, k2);

	start = std::chrono::high_resolution_clock::now();
	T DP0_band = maxmin(MatU, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff_band = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP_band_omp = maxmin_omp(MatU, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff_band_omp = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP_band_omp_s = maxmin_omp_band(MatU, N, k1, k2);
	end = std::chrono::high_resolution_clock::now();
	auto diff_band_omp_s = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	//block diagonal
	std::memcpy(MatU, Mat, sizeof(T) * N * M);
	int l = 20;
	int* offsets = new int[l];
	generate_offsets(offsets, N, l);
	block_diag_matrix(MatU, N, l, offsets);

	start = std::chrono::high_resolution_clock::now();
	T DP0_blockdiag = maxmin(MatU, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff_blockdiag = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP_blockdiag_omp = maxmin_omp(MatU, N, M);
	end = std::chrono::high_resolution_clock::now();
	auto diff_blockdiag_omp = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	start = std::chrono::high_resolution_clock::now();
	T DP_blockdiag_omp_s = maxmin_omp_block_diag(MatU, N, l, offsets);
	end = std::chrono::high_resolution_clock::now();
	auto diff_blockdiag_omp_s = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

	if (!silent) {
		std::clog << "time(nano):\t\t\t" << diff << std::endl;
		std::clog << "time(nano) omp:\t\t\t" << diff1 << std::endl;
		std::clog << "time(nano) triangle:\t\t" << diff_upper << std::endl;
		std::clog << "time(nano) triangle omp:\t" << diff_upper_omp << std::endl;
		std::clog << "time(nano) triangle omp spc.:\t" << diff_upper_omp_s << std::endl;
		std::clog << "time(nano) band: \t\t" << diff_band << std::endl;
		std::clog << "time(nano) band omp:\t\t" << diff_band_omp << std::endl;
		std::clog << "time(nano) band omp spc.:\t" << diff_band_omp_s << std::endl;
		std::clog << "time(nano) blockdiag:\t\t" << diff_blockdiag << std::endl;
		std::clog << "time(nano) blockdiag omp:\t" << diff_blockdiag_omp << std::endl;
		std::clog << "time(nano) blockdiag omp spc.:\t" << diff_blockdiag_omp_s << std::endl;


		if (std::abs(DP0 - DP) <= eps 
			&& std::abs(DP0_upper - DP_upper_omp) <= eps
			&& std::abs(DP0_upper - DP_upper_omp_s) <= eps
			&& std::abs(DP0_band - DP_band_omp) <= eps
			&& std::abs(DP0_band - DP_band_omp_s) <= eps
			&& std::abs(DP0_blockdiag - DP_blockdiag_omp) <= eps
			&& std::abs(DP0_blockdiag - DP_blockdiag_omp_s) <= eps)
			std::cout << "minmax found OK: " << DP << std::endl;
		else
			std::cout << "Error: " /* << DP << "; Should be: " << DP0*/ << std::endl;
	}
	else {
		if (std::abs(DP0 - DP) <= eps
			&& std::abs(DP0_upper - DP_upper_omp) <= eps
			&& std::abs(DP0_upper - DP_upper_omp_s) <= eps
			&& std::abs(DP0_band - DP_band_omp) <= eps
			&& std::abs(DP0_band - DP_band_omp_s) <= eps
			&& std::abs(DP0_blockdiag - DP_blockdiag_omp) <= eps
			&& std::abs(DP0_blockdiag - DP_blockdiag_omp_s) <= eps)	{ 
				std::cout << N << " " << cores << " ";
				std::cout << diff << " " << diff1 << " " << diff_upper << " " << diff_upper_omp << " " << diff_upper_omp_s << " " 
					<< diff_band << " " << diff_band_omp << " " << diff_band_omp_s << " "
					<< diff_blockdiag << " " << diff_blockdiag_omp << " " << diff_blockdiag_omp_s << std::endl;
			}
		else
			std::cout << 0 << " " << 0 << std::endl;
	}
	delete[] Mat;
	delete[] MatU;
	return 0;
}