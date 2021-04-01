#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace std;

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
	time_t time = (end - start) / CLOCKS_PER_SEC;
	printf("Time it took to run the function: %f\n", (double) time);
}

// ---------------------INSERTION SORT--------------------
// --------------------------- C -------------------------

void insertionSortC(int* arr, unsigned int n) {
	int* sorted_arr = arr;
	unsigned int j;
	int key;
	for (unsigned int i = 0; i < n; i++) {
		// Reference number
		key = sorted_arr[i];
		// Inferior numbers to reference
		j = i - 1;
		
		// All values inferior to key go to his left side
		while (j >= 0 && sorted_arr[j] > key) {
			sorted_arr[j + 1] = sorted_arr[j];
			j--;
		}

		arr[j + 1] = key;
	}
}

// --------------------- INSERTION SORT --------------------
// --------------------- ASSEMBLER X86 -----------------------

int* insertionSortAssembler(int* arr, int size_arr) {
	int* arr_sorted = arr;
	__asm {
		push ebp
		//mov ebp, esp
		push ebx
		push esi
		push edi

		mov esi, arr_sorted // array
		mov eax, 1 // i
		mov ebx, 0 // j
		mov ecx, size_arr // number of numbers
		
_loop_1:
			cmp eax, ecx
			//jnbe _end_loop_1
			jbe _end_loop_1

			push ecx // save number of array

			// ecx = array[i]
			mov ecx, [esi + eax * 4]

			// j = i - 1
			mov ebx, eax
			sub ebx, 1

_loop_2:
					// if j < 0 exit
					cmp ebx, 0
					//jnae _end_loop_2
					jae _end_loop_2

					// if arr[j] <= key exit
					cmp [esi + ebx * 4], ecx
					//jnb _end_loop_2
					jb _end_loop_2

					// array[j+1] = array[j]
					push [esi + ebx * 4]
					pop [esi + ebx * 4 + 4]

					// j--
					sub ebx, 1

					jmp _loop_2
_end_loop_2:

				// array[j+1] = key
				mov [esi + ebx*4 + 4], ecx

				// i++
				add eax, 1

				// restore items
				pop ecx

				jmp _loop_1
_end_loop_1:

			pop edi
			pop esi
			pop ebx
			pop ebp
	}

	return arr_sorted;
}

int main() {
	clock_t start, end;
	unsigned int size = 20;
	for (unsigned int size = 20; size <= 100; size *= 2) {
		int* arr = randomArrayNumbers(size);

		printf_s("Random numbers size %d array\n", size);
		print_arr(arr, size);
		printf_s("\n");

		// InsertSort C algorithm
		printf_s("C sorted array size %d with insertion sort algorithm\n", size);
		int* sortedArrC = new int[size];
		for (unsigned int i = 0; i < size; i++) {
			sortedArrC[i] = arr[i];
		}
		start = clock();
		insertionSortC(sortedArrC, size);
		end = clock();
		print_arr(sortedArrC, size);
		print_time(start, end);
		printf_s("\n\n");
		// delete sortedArrC;

		// InsertSort Assembler algorithm
		printf_s("Assembler sorted array size %d with insertion sort algorithm\n", size);
		int* sortedArrAssem = new int[size];
		for (unsigned int i = 0; i < size; i++) {
			sortedArrAssem[i] = arr[i];
		}
		start = clock();
		int* arr2 = insertionSortAssembler(sortedArrAssem, size);
		end = clock();
		print_arr(arr2, size);
		print_time(start, end);

		if (sortedArrC != arr2) {
			printf("NOT EQUAL");
			return -1;
		}
		printf_s("\n\n");
		// delete sortedArrAssem;\n;
		printf_s("-----------------------------------------\n\n");
	}
}