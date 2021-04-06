#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include <synchapi.h>

using namespace std;

bool equal_arr(int* arr1, int* arr2, unsigned int arr_size) {
	for (unsigned int i = 0; i < arr_size; i++) {
		if (arr1[i] != arr2[i]) {
			return false;
		}
	}

	return true;
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

void print_arr(int* arr, unsigned int size) {
	printf_s("\n");

	for (unsigned int i = 0; i < size; i++) {
		printf_s("%d ", arr[i]);
	}

	printf_s("\n");
}

void print_time(clock_t start, clock_t end) {
	printf_s("Time it took to run the function: %.6f", (float)(end - start) / CLOCKS_PER_SEC);
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
		mov esi, sorted_arr // empty arr
		mov eax, 0 // first number
		mov ebx, 1 // second number
		mov ecx, n // counter
		sub ecx, 2

		mov [esi], eax
		add esi, 4
		mov [esi], ebx
		add esi, 4

		jmp _end_loop_1
_loop_1:
		mov eax, [esi - 4]
		mov ebx, [esi - 8]
		mov edx, ebx
		add edx, eax
		mov [esi], edx
		add esi, 4
		sub ecx, 1
_end_loop_1:
		cmp ecx, 0
		ja _loop_1

	}
}

int main() {
	clock_t start, end;
	for(unsigned int i=10; i<1000000000; i*=10) {
		int size = i;
		// Fibonnacci Part 1
		// Fibonacci C algorithm
		printf_s("C sorted array size %d with fibonacci algorithm\n", size);
		int* arr_C = new int[size];
		start = clock();
		fibonacci_C(arr_C, size);
		end = clock();
		// print_arr(arr_C, size);
		print_time(start, end);
		printf_s("\n\n");

		// Fibonacci Assembler algorithm
		printf_s("Assembler array size %d with fibonacci algorithm\n", size);
		int* arr_assemb = new int[size];
		start = clock();
		fibonacci_assemb(arr_assemb, size);
		end = clock();
		// print_arr(arr_assemb, size);
		print_time(start, end);
		printf_s("\n\n");

		if (!equal_arr(arr_C, arr_assemb, size)) {
			printf("NOT EQUAL\n");
			return -1;
		}

		delete arr_C;
		delete arr_assemb;
		printf_s("-----------------------------------------\n\n");
	}

	// Fibonnacci Part 2
	for (unsigned int i = 10; i < 1000000000; i *= 10) {
		int size = i;
		// Fibonacci C algorithm
	}
}