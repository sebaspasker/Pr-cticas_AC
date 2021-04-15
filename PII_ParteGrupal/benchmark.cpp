#include <stdlib.h>
#include <iostream>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include <synchapi.h>
#include <cmath>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <emmintrin.h>
#include <vector>

using namespace std;

// int, short, double, float order
struct timer_struct {
	clock_t random_matrix_creation;
	clock_t pure_C_time;
	clock_t pure_assem_time;
	clock_t C_row_col_time[4];
	clock_t SSE_time[2];
	clock_t SSE2_time;
	clock_t SSE3_time;
	clock_t fibb_C_time;
	clock_t fibb_assem_time;
};

timer_struct timer;

// ----------------------- UTILS ---------------------------------------
// ---------------------------------------------------------------------

bool equal_arr(int* arr1, int* arr2, unsigned int arr_size) {
	for (unsigned int i = 0; i < arr_size; i++) {
		if (arr1[i] != arr2[i]) {
			return false;
		}
	}

	return true;
}

int get_random_val() {
	return (int)rand() % 10 + 1;
}

void create_random_matrix(int** matrix, unsigned int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		matrix[i] = new int[size];
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			int* row = matrix[i];
			__asm {
				mov esi, row
				mov eax, j
				mov ebx, get_random_val
				mov [esi + eax * 4], ebx
			}
		}
	}
}

int* randomArrayNumbers(unsigned int size) {
	int* arr = new int[size];
	time_t t;
	srand((unsigned) time(&t));

	for (unsigned int i = 0; i < size; i++) {
		arr[i] = rand() % size + 1;
	};

	return arr;
}

// ----------------------- PRINT ---------------------------------------
// ---------------------------------------------------------------------

void print_arr(int* arr, unsigned int size) {
	printf_s("\n");

	for (unsigned int i = 0; i < size; i++) {
		printf_s("%d ", arr[i]);
	}

	printf_s("\n");
}

void print_matrix(int** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			printf_s("%d ", matrix[i][j]);
		}
		printf_s("\n");
	}
}

void print_double_matrix(double** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			printf_s("%.0f ", matrix[i][j]);
		}
		printf_s("\n");
	}
}

void print_short_matrix(short** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			printf_s("%hu ", matrix[i][j]);
		}
		printf_s("\n");
	}
}

void print_float_matrix(float** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			printf_s("%0.f ", matrix[i][j]);
		}
		printf_s("\n");
	}
}

void print_timer_mult(bool overvalue) {
	printf_s("********************************************************\n");
	printf_s("Matrix multiplication ejecution time:\n\n");
	printf_s("Matrix random values creation time\n");
	printf_s("INT: %.6f\n\n", (float)timer.random_matrix_creation / CLOCKS_PER_SEC);

	if (!overvalue) {
		printf_s("Pure C Algoritm:\n");
		printf_s("INT: %.6f\n\n", (float)timer.pure_C_time / CLOCKS_PER_SEC);

		printf_s("ROW-COL C Algorithm (without inverse matrix calcul)\n");
		printf_s("INT: %.6f\n", (float)timer.C_row_col_time[0] / CLOCKS_PER_SEC);
		printf_s("SHORT: %.6f\n", (float)timer.C_row_col_time[1] / CLOCKS_PER_SEC);
		printf_s("FLOAT: %.6f\n", (float)timer.C_row_col_time[2] / CLOCKS_PER_SEC);
		printf_s("DOUBLE: %.6f\n\n", (float)timer.C_row_col_time[3] / CLOCKS_PER_SEC);
	}

	printf_s("Pure Assembler X86 Algorithm:\n");
	printf_s("INT: %.6f\n\n", (float)timer.pure_assem_time / CLOCKS_PER_SEC);

	printf_s("Assembler X86 + SSE Algorithm:\n");
	printf_s("FLOAT: %.6f\n", (float)timer.SSE_time[0] / CLOCKS_PER_SEC);
	printf_s("DOUBLE: %.6f\n\n", (float)timer.SSE_time[1] / CLOCKS_PER_SEC);

	printf_s("Assembler X86 + SSE2 Algorithm:\n");
	printf_s("SHORT: %.6f\n\n", (float)timer.SSE2_time / CLOCKS_PER_SEC);

	printf_s("Assembler X86 + SSE3 Algorithm:\n");
	printf_s("FLOAT: %.6f\n\n", (float)timer.SSE3_time / CLOCKS_PER_SEC);
	printf_s("********************************************************\n");
}

// 1: INT, 2:SHORT, 3:FLOAT, 4:DOUBLE, 5:FIBONACCI
void print_individual_timer(bool int_, bool short_, bool float_, bool double_, bool fibb_, bool overvalue) {
	printf_s("********************************************************\n");
	if (int_ || short_ || float_ || double_) {
		printf_s("Matrix multiplication ejecution time:\n\n");
		printf_s("Matrix random values creation time\n");
		printf_s("INT: %.6f\n", (float)timer.random_matrix_creation / CLOCKS_PER_SEC);
	}

	if(!overvalue) {
		if (int_) {
			printf_s("Pure C algorithm\n");
			printf_s("INT: %.6f\n", (float)timer.pure_C_time / CLOCKS_PER_SEC);
		}

		if (int_ || short_ || float_ || double_) {
			printf_s("ROW-COL C Algorithm (without inverse matrix calcul)\n");
		}
		if(int_) printf_s("INT: %.6f\n", (float)timer.C_row_col_time[0] / CLOCKS_PER_SEC);
		if(short_) printf_s("SHORT: %.6f\n", (float)timer.C_row_col_time[1] / CLOCKS_PER_SEC);
		if(float_) printf_s("FLOAT: %.6f\n", (float)timer.C_row_col_time[2] / CLOCKS_PER_SEC);
		if(double_) printf_s("DOUBLE: %.6f\n", (float)timer.C_row_col_time[3] / CLOCKS_PER_SEC);
	}

	if (int_) {
		printf_s("Pure Assembler X86 Algorithm:\n");
		printf_s("INT: %.6f\n", (float)timer.pure_assem_time / CLOCKS_PER_SEC);
	}
	
	if(float_ || double_) printf_s("Assembler X86 + SSE Algorithm:\n");
	if(float_) printf_s("FLOAT: %.6f\n", (float)timer.SSE_time[0] / CLOCKS_PER_SEC);
	if(double_) printf_s("DOUBLE: %.6f\n", (float)timer.SSE_time[1] / CLOCKS_PER_SEC);

	if (short_) {
		printf_s("Assembler X86 + SSE2 Algorithm:\n");
		printf_s("SHORT: %.6f\n", (float)timer.SSE2_time / CLOCKS_PER_SEC);
	}

	if (float_) {
		printf_s("Assembler X86 + SSE3 Algorithm:\n");
		printf_s("FLOAT: %.6f\n", (float)timer.SSE3_time / CLOCKS_PER_SEC);
	}

	if (fibb_) {
		printf_s("Fibonacci ejecution time:\n");
		printf_s("C: %.6f\n", (float)timer.fibb_C_time / CLOCKS_PER_SEC);
		printf_s("Assembler X86: %.6f\n", (float)timer.fibb_assem_time / CLOCKS_PER_SEC);
	}
	printf_s("********************************************************\n");
}

void print_time(clock_t start, clock_t end) {
	printf_s("Time it took to run the function: %.6f", (float)(end - start) / CLOCKS_PER_SEC);
}

// ----------------------- DELETE --------------------------------------
// ---------------------------------------------------------------------

void delete_int_matrix(int** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) delete[] matrix[i];
	delete[] matrix;
}

void delete_double_matrix(double** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) delete[] matrix[i];
	delete[] matrix;
}

void delete_short_matrix(short** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) delete[] matrix[i];
	delete[] matrix;
}

void delete_float_matrix(float** matrix, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) delete[] matrix[i];
	delete[] matrix;
}

// ----------------------- FIBONACCI ---------------------
// --------------------------- C -------------------------

void fibonacci_C(int* arr_, unsigned int n) {
	int t1, t2, sum;
	t1 = 0;
	t2 = 1;
	arr_[0] = t1;
	arr_[1] = t2;
	for (unsigned int i = 2; i < n; i++) {
		sum = t1 + t2;
		arr_[i] = sum;
		t1 = t2;
		t2 = sum;
	}
}

// ----------------------- FIBONACCI ---------------------
// --------------------- ASSEMBLER X86 -------------------

void fibonacci_assemb(int* arr, unsigned int n) {
	int* sorted_arr = arr;
	__asm {
		pusha // save state

		mov esi, sorted_arr // empty arr
		mov eax, 0 // first number
		mov ebx, 1 // second number
		mov ecx, n // counter
		sub ecx, 2 // rest first two numbers

		// save first values
		mov [esi], eax
		add esi, 4
		mov [esi], ebx
		add esi, 4

		jmp _end_loop_1
_loop_1:
		mov eax, [esi - 4]
		mov ebx, [esi - 8]

		// a = arr[i] + arr[i+1]
		mov edx, ebx
		add edx, eax

		// save edx in arr
		mov [esi], edx

		// i++
		add esi, 4

		// n--
		sub ecx, 1
_end_loop_1:
		cmp ecx, 0
			ja _loop_1

		popa
	}
}

// -------------------------------- MATRIX MULTIPLICATION --------------------
// --------------------------------	    	C			  --------------------

void mult_matrix_C(int** arr, int** arr_mult, int size) {
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			arr_mult[i][j] = 0;
			for (unsigned int k = 0; k < (unsigned)size; k++) {
				arr_mult[i][j] = arr_mult[i][j] + arr[i][k] * arr[k][j];
			}
		}
	}
}

// -------------------------------- MATRIX MULTIPLICATION --------------------
// --------------------------------   ASSEMBLER X86       --------------------

void mult_matrix_assem(int** arr, int** arr_mult, int size) {
	int i = 0, j = 0, k = 0;
	int max_size = size;
	int value = 0;
	__asm {
		mov esi, arr
_loop_1:
		mov j, 0
_loop_2:
		mov value, 0
		mov k, 0
_loop_3:
		// arr[i][k]
		mov ebx, i
		mov edi, [esi + ebx*4]
		mov ebx, k
		mov eax, [edi + ebx*4]

		// arr[k][j]
		mov edi, [esi + ebx*4]
		mov ebx, j
		mov ecx, [edi + ebx*4]

		// arr[i][k] * arr[k][j]
		mul ecx
		add value, eax
		
		// k < size
		mov eax, k
		mov ecx, max_size
		add eax, 1 // k++
		cmp eax, ecx
		je _end_loop_3
		mov k, eax
		jmp _loop_3
_end_loop_3:
		// save edx in arr[i][j]
		mov ebx, i
		mov ecx, arr_mult
		mov edi, [ecx + ebx * 4]
		mov ebx, j
		mov edx, value
		mov [edi + ebx * 4], edx

		// j < size
		mov eax, j
		mov ecx, max_size
		add eax, 1 // j++
		cmp eax, ecx
		je _end_loop_2
		mov j, eax
		jmp _loop_2
_end_loop_2:
		// j < size
		mov eax, i
		mov ecx, max_size
		add eax, 1 // i++
		cmp eax, ecx
		je _end_loop_1
		mov i, eax
		jmp _loop_1
_end_loop_1:
	}
}

// ----------------------- ROW-COL MATRIX ------------------------------------------
// ----------------------- MULTIPLICATION ------------------------------------------
// -----------------------      INT       ------------------------------------------

int mult_sum_int(const int* row, const int* col, unsigned int size) {
	int res = 0;
	for (unsigned int i = 0; i < size; i++) {
		res += row[i] * col[i];
	}

	return res;
}

// ----------------------- ROW-COL MATRIX ------------------------------------------
// ----------------------- MULTIPLICATION ------------------------------------------
// -----------------------    DOUBLE     -------------------------------------------

double mult_sum_double(const double* row, const double* col, unsigned int size) {
	double res = 0.0;
	for (unsigned int i = 0; i < size; i++) {
		res += row[i] * col[i];
	};

	return res;
}

double mult_sum_double_SSE(const double *row_matrix, const double *col_matrix,unsigned int size) {
	double val = 0.0, res = 0.0;
	__declspec(align(16)) double tmp[2] = { 0.0, 0.0 };
	__m128d mres;

	if ((size / 2) != 0) {
		mres = _mm_load_sd(&val);
		for (unsigned int i = 0; i < size / 2; i++) {
			mres = _mm_add_pd(mres, _mm_mul_pd(_mm_loadu_pd(&row_matrix[2 * i]), _mm_loadu_pd(&col_matrix[2 * i])));
		}

		_mm_store_pd(tmp, mres);

		res = tmp[0] + tmp[1];
	}

	if ((size % 2) != 0) {
		for (unsigned int i = size - size % 2; i < size; i++) {
			res += row_matrix[i] * col_matrix[i];
		}
	}

	return res;
}

// ----------------------- ROW-COL MATRIX ------------------------------------------
// ----------------------- MULTIPLICATION ------------------------------------------
// -----------------------     SHORT      -------------------------------------------

short mult_sum_short(short* row, short* col, int size) {
	short sum = 0;
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		sum += row[i] * col[i];
	}

	return sum;
}

short mult_sum_short_SSE2(const short* row, const short* col, int size) {
	__m128i* mem_pointer_row = (__m128i*)row;
	__m128i* mem_pointer_col = (__m128i*)col;
	__m128i mem_res = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 0);

	for (unsigned int i = 0; i < (unsigned)size / 8; i++) {
		__m128i mem_tmp = _mm_mullo_epi16(_mm_loadu_si128(mem_pointer_row), 
			_mm_loadu_si128(mem_pointer_col));
		mem_res = _mm_add_epi16(mem_res, mem_tmp);
		mem_pointer_row++;
		mem_pointer_col++;
	}

	short res[8];
	__m128i* pointer_mem_res = (__m128i*)res;
	_mm_storeu_si128(pointer_mem_res, mem_res);
	short tot = 0;
	for (unsigned int i = 0; i < 8; i++) tot += res[i];
	return tot;
}

// ----------------------- ROW-COL MATRIX ------------------------------------------
// ----------------------- MULTIPLICATION ------------------------------------------
// -----------------------     FLOAT      -------------------------------------------


// Regular C row-col multiplication
float mult_sum_float(const float* row, const float* col, unsigned int size) {
	float tot = 0.0f;
	for (unsigned int i = 0; i < size; i++) {
		tot += row[i] * col[i];
	}

	return tot;
}

// SSE row-col multiplication
float mult_sum_float_SSE(const float* row_matrix, const float* col_matrix, unsigned int size) {
	float val = 0.0, res = 0.0;
	__declspec(align(16)) float tmp[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
	__m128 mres;

	if ((size / 2) != 0) {
		mres = _mm_load_ss(&val);
		for (unsigned int i = 0; i < size / 4; i++) {
			mres = _mm_add_ps(mres, _mm_mul_ps(_mm_loadu_ps(&row_matrix[4 * i]), _mm_loadu_ps(&col_matrix[4 * i])));
		}

		__m128 mv1 = _mm_movelh_ps(mres, mres); // a, b, a, b
		__m128 mv2 = _mm_movehl_ps(mres, mres); // c, d, c, d
		mres = _mm_add_ps(mv1, mv2);

		_mm_store_ps(tmp, mres);

		res = tmp[0] + tmp[1];
	}

	if ((size % 4) != 0) {
		for (unsigned int i = size - size % 4; i < size; i++) {
			res += row_matrix[i] * col_matrix[i];
		}
	}

	return res;
}

// SSE3 row-col multiplication
float mult_sum_float_SSE3(const float* row_matrix, const float* col_matrix, unsigned int size) {
	float val = 0.0f, res = 0.0f;

	if ((size / 4) != 0) {
		const float* pointer_row = row_matrix;
		const float* pointer_col = col_matrix;

		__asm {
			movss xmm0, xmmword ptr[val];
		}

		for (unsigned int i = 0; i < size / 4; i++) {
			__asm {
				mov eax, dword ptr[pointer_row]
				mov ebx, dword ptr[pointer_col]
				movups xmm1, [eax]
				movups xmm2, [ebx]
				mulps xmm1, xmm2
				addps xmm0, xmm1
			}

			pointer_row += 4;
			pointer_col += 4;
		}
		__asm {
			haddps xmm0, xmm0
			haddps xmm0, xmm0
			movss dword ptr[res], xmm0
		}
	}

	return res;
}

// ---------------------------------------------------------------------------------
// ------------------------ MATRIX MULTIPLICATION BENCHMARKS -----------------------

// I N T
void int_benchmark(int** arr_matrix, unsigned int size, bool print_arr, const bool overvalue) {
	clock_t start, end;
	int** matrix = new int*[size];
	int** mult_C = new int*[size];
	int** mult_C2 = new int*[size];
	int** mult_assem = new int*[size];
	int** matrix_inverse = new int*[size];

	for (unsigned int i = 0; i < (unsigned)size; i++) {
		matrix[i] = new int[size];
		mult_assem[i] = new int[size];
		mult_C[i] = new int[size];
		mult_C2[i] = new int[size];
		matrix_inverse[i] = new int[size];
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			matrix[i][j] = arr_matrix[i][j]; // Create security copy
			mult_C[i][j] = 0;
			mult_C2[i][j] = 0;
			mult_assem[i][j] = 0;
			matrix_inverse[i][j] = 0;
		}
	}

	// Create invers matrix 
	for (unsigned int i = 0; i < (unsigned)size; i++) {
		for (unsigned int j = 0; j < (unsigned)size; j++) {
			matrix_inverse[i][j] = matrix[j][i];
		}
	}

	if (!overvalue) {
		// Pure C matrix multiplication algorithm
		start = clock();
		mult_matrix_C(matrix, mult_C, size);
		end = clock();
		timer.pure_C_time = end - start;

		start = clock();
		// Row-Col C multiplication
		for (unsigned int i = 0; i < (unsigned)size; i++) {
			for (unsigned int j = 0; j < (unsigned)size; j++) {
				int* row = matrix[i];
				int* col = matrix_inverse[j];
				mult_C2[i][j] = mult_sum_int(row, col, size);
			}
		}
		end = clock();
		timer.C_row_col_time[0] = end - start;
	}

	start = clock();
	// Assembler X86 matrix multiplication
	mult_matrix_assem(matrix, mult_assem, size);
	end = clock();
	timer.pure_assem_time = end - start;


	if (print_arr) {
		printf_s("INT\n\n");
		print_matrix(matrix, size);
		printf_s("\n\n");
		print_matrix(mult_C, size);
		printf_s("\n\n");
		print_matrix(mult_C2, size);
		printf_s("\n\n");
		print_matrix(mult_assem, size);
		printf_s("\n\n");
	}

	delete_int_matrix(matrix, size);
	delete_int_matrix(mult_C, size);
	delete_int_matrix(mult_C2, size);
	delete_int_matrix(mult_assem, size);
	delete_int_matrix(matrix_inverse, size);
}

// S H O R T
void short_benchmark(int** standard_matrix, unsigned int size, bool print_arr, const bool overvalue) {
	clock_t start, end;
	short** matrix = new short* [size];
	short** inverse_matrix = new short* [size];
	short** mult_matrix = new short* [size];
	short** mult_matrix_SSE = new short* [size];

	for (unsigned int i = 0; i < size; i++) {
		matrix[i] = new short[size];
		inverse_matrix[i] = new short[size];
		mult_matrix_SSE[i] = new short[size];
		mult_matrix[i] = new short[size];
		for (unsigned int j = 0; j < size; j++) {
			matrix[i][j] = (short)standard_matrix[i][j];
			inverse_matrix[i][j] = (short)standard_matrix[j][i];
			mult_matrix_SSE[i][j] = (short)0.0;
			mult_matrix[i][j] = (short)0.0;
		}
	}
	
	short* row, *col; 
	if (!overvalue) {
		start = clock();
		// Simple C short calculation
		for (unsigned int i = 0; i < size; i++) {
			for (unsigned int j = 0; j < size; j++) {
				row = matrix[i];
				col = inverse_matrix[j];
				mult_matrix[i][j] = mult_sum_short(row, col, size);
			}
		}
		end = clock();
		timer.C_row_col_time[1] = end - start;
	}

	// SSE2 short calculation
	start = clock();
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			row = matrix[i];
			col = inverse_matrix[j];
			mult_matrix_SSE[i][j] = mult_sum_short_SSE2(row, col, size);
		}
	}
	end = clock();
	timer.SSE2_time = end - start;

	if (print_arr) {
		printf_s("SHORT\n\n");
		print_short_matrix(matrix, size);
		printf_s("\n\n");
		print_short_matrix(mult_matrix, size);
		printf_s("\n\n");
		print_short_matrix(mult_matrix_SSE, size);
		printf_s("\n\n");
	}

	delete_short_matrix(matrix, size);
	delete_short_matrix(inverse_matrix, size);
	delete_short_matrix(mult_matrix, size);
	delete_short_matrix(mult_matrix_SSE, size);
}

// F L O A T
void float_benchmark(int** standard_matrix, unsigned int size, bool print_arr, const bool overvalue) {
	clock_t start, end;
	float** matrix = new float*[size];
	float** inverse_matrix = new float*[size];
	float** mult_matrix = new float* [size];
	float** mult_matrix_SSE = new float* [size];
	float** mult_matrix_SSE3 = new float* [size];

	for (unsigned int i = 0; i < size; i++) {
		matrix[i] = new float[size];
		inverse_matrix[i] = new float[size];
		mult_matrix_SSE[i] = new float[size];
		mult_matrix_SSE3[i] = new float[size];
		mult_matrix[i] = new float[size];
		for (unsigned int j = 0; j < size; j++) {
			matrix[i][j] = (float)standard_matrix[i][j];
			inverse_matrix[i][j] = (float)standard_matrix[j][i];
			mult_matrix_SSE[i][j] = (float)0.0;
			mult_matrix_SSE3[i][j] = (float)0.0;
			mult_matrix[i][j] = (float)0.0;
		}
	}

	float* col, * row;
	// C matrix multiplication
	if (!overvalue) {
		start = clock();
		for (unsigned int i = 0; i < size; i++) {
			for (unsigned int j = 0; j < size; j++) {
				row = matrix[i];
				col = inverse_matrix[j];
				mult_matrix[i][j] = mult_sum_float(row, col, size);
			}
		}
		end = clock();
		timer.C_row_col_time[2] = end - start;
	}

	start = clock();
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			row = matrix[i];
			col = inverse_matrix[j];
			mult_matrix_SSE[i][j] = mult_sum_float_SSE(row, col, size);
		}
	}
	end = clock();
	timer.SSE_time[0] = end - start;

	start = clock();
	for (unsigned int i = 0; i < size; i++) { 
		for (unsigned int j = 0; j < size; j++) {
			row = matrix[i];
			col = inverse_matrix[j];
			mult_matrix_SSE3[i][j] = mult_sum_float_SSE3(row, col, size);
		}
	}
	end = clock();
	timer.SSE3_time = end - start;

	if (print_arr) {
		printf_s("FLOAT\n\n");
		print_float_matrix(matrix, size);
		printf_s("\n\n");
		print_float_matrix(mult_matrix, size);
		printf_s("\n\n");
		print_float_matrix(mult_matrix_SSE, size);
		printf_s("\n\n");
		print_float_matrix(mult_matrix_SSE3, size);
		printf_s("\n\n");
	}

	delete_float_matrix(matrix, size);
	delete_float_matrix(inverse_matrix, size);
	delete_float_matrix(mult_matrix, size);
	delete_float_matrix(mult_matrix_SSE, size);
	delete_float_matrix(mult_matrix_SSE3, size);
}

// D O U B L E
void double_benchmark(int** standard_matrix, unsigned int size, bool print_arr, const bool overvalue) {
	clock_t start, end;
	double** matrix = new double*[size];
	double** inverse_matrix = new double*[size];
	double** mult_matrix = new double* [size];
	double** mult_matrix_SSE = new double* [size];

	for (unsigned int i = 0; i < size; i++) {
		matrix[i] = new double[size];
		inverse_matrix[i] = new double[size];
		mult_matrix_SSE[i] = new double[size];
		mult_matrix[i] = new double[size];
		for (unsigned int j = 0; j < size; j++) {
			matrix[i][j] = (double)standard_matrix[i][j];
			inverse_matrix[i][j] = (double)standard_matrix[j][i];
			mult_matrix_SSE[i][j] = (double)0.0;
			mult_matrix[i][j] = (double)0.0;
		}
	}
	
	double* row, *col; 
	if (!overvalue) {
		start = clock();
		// Simple C calculation
		for (unsigned int i = 0; i < size; i++) {
			for (unsigned int j = 0; j < size; j++) {
				row = matrix[i];
				col = inverse_matrix[j];
				mult_matrix[i][j] = mult_sum_double(row, col, size);
			}
		}
		end = clock();
		timer.C_row_col_time[3] = end - start;
	}

	start = clock();
	// SSE double optimized calculation
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			row = matrix[i];
			col = inverse_matrix[j];
			mult_matrix_SSE[i][j] = mult_sum_double_SSE(row, col, size);
		}
	}
	end = clock();
	timer.SSE_time[1] = end - start;

	if (print_arr) {
		printf_s("DOUBLE\n\n");
		print_double_matrix(matrix, size);
		printf_s("\n\n");
		print_double_matrix(mult_matrix, size);
		printf_s("\n\n");
		print_double_matrix(mult_matrix_SSE, size);
		printf_s("\n\n");
	}

	delete_double_matrix(matrix, size);
	delete_double_matrix(inverse_matrix, size);
	delete_double_matrix(mult_matrix_SSE, size);
	delete_double_matrix(mult_matrix, size);
}

// ---------------------------------------------------------------------------------
// ------------------------  FIBONACCI BENCHMARK  ----------------------------------

void fibonacci_benchmark() {
	clock_t start, end;
	int size = (int)pow(10, 9);
	for (unsigned int i = 10; i < pow(10,9); i *= 10) {
		int size = i;

		// Fibonacci C algorithm -------------------------------

		printf_s("C sorted array size %d with fibonacci algorithm\n", size);
		int* arr_C = new int[size];
		start = clock();
		fibonacci_C(arr_C, size);
		end = clock();
		print_time(start, end);
		printf_s("\n\n");
		timer.fibb_C_time = end - start;

		// Fibonacci Assembler algorithm -----------------------

		printf_s("Assembler array size %d with fibonacci algorithm\n", size);
		int* arr_assemb = new int[size];
		start = clock();
		fibonacci_assemb(arr_assemb, size);
		end = clock();
		print_time(start, end);
		printf_s("\n\n");
		timer.fibb_assem_time = end - start;
		

		if (!equal_arr(arr_C, arr_assemb, size)) {
			printf("NOT EQUAL\n");
		}

		delete[] arr_C;
		delete[] arr_assemb;
		printf_s("-----------------------------------------\n\n");
	}
}

void print_options(int& option, char& print) {
	printf_s("Welcome, Options:\n\n");
	printf_s("1. INT Matrix multiplication benchmark\n");
	printf_s("2. SHORT Matrix multiplication benchmark\n");
	printf_s("3. FLOAT Matrix multiplication benchmark\n");
	printf_s("4. DOUBLE Matrix multiplication benchmark\n");
	printf_s("5. Fibbonaci benchmark\n");
	printf_s("6. All Matrix multiplication benchmark (default)\n");
	printf_s("\n");
	cin >> option;
	cin.ignore();
	printf_s("Print matrix: (y/n) default no:\n");
	print = getchar();
	printf_s("\n");
}

int main() {
	bool print_arr = false;
	int option;
	char print;
	print_options(option, print);
	if (print == 'y') print_arr = true;


	if(option != 5) {
		printf_s("Choose matrix max size multiple of 8 (max_size = size*8) (max: 4):\n");
		int size;
		cin >> size;
		cin.ignore();
		if (size > 6 || size < 0) size = 4;
		clock_t start, end;
		bool overvalue = false;
		for (unsigned int i = 1; i < (unsigned)size; i++) {
			unsigned int max_size = pow(8, i);
			int** arr_matrix = new int* [max_size];
			overvalue = i > 4;

			srand((unsigned)time(NULL));
			start = clock();
			create_random_matrix(arr_matrix, max_size);
			end = clock();
			timer.random_matrix_creation = end - start;

			printf_s("Matrix size: %d\n", max_size);
			if(option == 1 || option == 6) int_benchmark(arr_matrix, max_size, print_arr, overvalue);
			if(option == 2 || option == 6) short_benchmark(arr_matrix, max_size, print_arr, overvalue);
			if(option == 3 || option == 6) double_benchmark(arr_matrix, max_size, print_arr, overvalue);
			if(option == 4 || option == 6) float_benchmark(arr_matrix, max_size, print_arr, overvalue);

			if(option == 6) print_timer_mult(overvalue);
			if (option == 1) print_individual_timer(true, false, false, false, false, overvalue);
			if (option == 2) print_individual_timer(false, true, false, false, false, overvalue);
			if (option == 3) print_individual_timer(false, false, true, false, false, overvalue);
			if (option == 4) print_individual_timer(false, false, false, true, false, overvalue);

			delete_int_matrix(arr_matrix, i);
		}
	}
	else {
		fibonacci_benchmark();
	}
}