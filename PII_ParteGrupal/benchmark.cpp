#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace std;

int* randomArrayNumbers(unsigned int size) {
	int* arr = new int[size];
	time_t t;
	srand((unsigned) time(&t));

	for (unsigned int i = 0; i < size; i++) {
		arr[i] = rand() % 20 + 1;
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

// -------------------- QUICKSORT ------------------------
// --------------------     C     ------------------------

void swap(int* a, int* b) {
	int x = *a;
	*a = *b;
	*b = x;
}

int partition(int* arr, int low_num, int high_num) {
	int pivot = arr[high_num];
	int i = (low_num - 1);

	for (unsigned int j = low_num; j <= ((unsigned) high_num - 1); j++) {
		if (arr[j] <= pivot) {
			i++;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high_num]);

	return i + 1;
}

void quickSortC_rec(int* arr, int low_num, int high_num) {
	if (low_num < high_num) {
		int part = partition(arr, low_num, high_num);
		quickSortC_rec(arr, low_num, part - 1);
		quickSortC_rec(arr, part + 1, high_num);
	}
}


int* quickSortC(int* arr, unsigned int n) {
	int* arr_sorted = arr;
	quickSortC_rec(arr_sorted, 0, n - 1);

	return arr_sorted;
}


// -------------------- QUICKSORT ------------------------
// -------------------- ASSEMBLER ------------------------

int* quickSortAssembler(int* arr, int n) {
	int* sort_arr = arr;
	// l, h, pivot
	int* values_quick = new int[3];
	values_quick[0] = 0; values_quick[1] = 0; values_quick[2] = 0;
	// pivot, i, j
	int* values_part = new int[3];
	values_part[0] = 0; values_part[1] = 0; values_part[2] = 0; 

	__asm {
//		mov edi, values_quick
		// l
//		mov [edi], 0
		// h
//		mov eax, n
//		sub edx, 1
//		mov [edi + 4], edx
		mov eax, 0
		mov ebx, n
		sub ebx, 1
		
// quicksort(eax, ebx)
quicksort_:
		// mov edx, values_quick
		// mov [edx + 4*0], eax
		// mov [edx + 4*1], ebx
		cmp eax, ebx
		jbe end
		call partition_
		push eax
		push ebx
		push ecx
		mov ebx, ecx
		call quicksort_
		pop ecx
		pop ebx
		pop eax
		add ecx, 1
		push eax
		push ebx
		push ecx
		mov eax, ecx
		call quicksort_
		pop ecx
		pop ebx
		pop eax
		/*mov edx, values_quick
		mov eax, [edx + 4*0]
		mov ebx, ecx
		mov [edx + 4*0], ebx
		call quicksort_*/
		jmp end


		// partition(b, h) {
		//	pivot = arr[b]
		//	i = b
		//	j = h
		//	while(i < j) {
		//		do {
		//			i++
		//		} while(A[i] <= pivot)
		//		
		//		do {
		//			j--
		//		} while(A[j] > pivot)
		//		
		//		if (i < j) {
		//			swap(A[i], A[j])
		//		} 
		//
		//	swap(A[b], A[h])
		//		
		//	return i + 1
		//	
		// }
		//
		// partition(eax, ebx) -> ecx
partition_:
		mov edi, sort_arr
		mov edx, [edi + 4 * eax]
		mov edi, values_part
		mov [edi + 4 * 0], edx // pivot
		mov [edi + 4 * 1], eax // i
		mov [edi + 4 * 2], ebx // j
		
loop_1: // while(i<j) {
		mov edx, [edi + 4 * 1] 
		mov ecx, [edi + 4 * 2]
		cmp edx, ecx
		jbe end_loop_1

		mov ecx, [edi + 4 * 0] // pivot
	
loop_2:	// do {} while(A[i] <= pivot)
		add edx, 1
		mov edi, sort_arr
		mov esi, [edi + 4 * edx] // A[i]
		cmp esi, ecx 
		jb loop_2

		mov edi, values_part
		mov [edi + 4 * 1], edx // save i
		mov edx, [edi + 4 * 2]

loop_3:	// do {} while(A[j] > pivot)
		sub edx, 1
		mov edi, sort_arr
		mov esi, [edi + 4 * edx] // A[j]
		cmp esi, ecx
		ja loop_3

		mov edi, values_part
		mov [edi + 4 * 2], edx // save j

		// if(i < j)
		mov ecx, [edi + 4 * 1]
		cmp ecx, edx
		jbe end_conditional_1
		// swap
		mov edi, sort_arr
		mov esi, [edi + 4 * ecx]
		mov esp, [edi + 4 * edx]
		mov [edi + 4 * ecx], esp
		mov [edi + 4 * edx], esi
end_conditional_1:
		jmp loop_1
end_loop_1:
		mov edi, sort_arr
		mov edx, [edi + 4 * eax]
		mov ecx, [edi + 4 * ebx]
		mov [edi + 4 * eax], ecx
		mov [edi + 4 * ebx], edx

		mov ecx, [values_part + 4 * 1]
		add ecx, 1
		ret


end:
		ret
/*
loop_partition_1: // while(i < j) 1
		cmp edx, edi 
		jbe end_loop_partition_1

loop_partition_2: // do while 2
		sub edx, 1
		cmp [esi + 4 * edx], ecx
		jb loop_partition_2 // end do while 2

loop_partition_3: // do while 3
		add edi, 1
		cmp [esi + 4 * edi], ecx
		ja loop_partition_3 // end do while 3

		cmp edx, edi // if (i < j)
		jbe loop_partition_1
		// swap(A[i], A[j])
		mov edi, [esi + 4 * eax]
		mov ebp, [esi + 4 * ebx]
		mov [esi + 4 * eax], ebp
		mov [esi + 4 * ebx], edi

		jmp loop_partition_1 // end while 1 
end_loop_partition_1:

		// swap(A[b], A[h])
		mov edi, [esi + 4 * eax]
		mov ebp, [esi + 4 * ebx]
		mov [esi + 4 * eax], ebp
		mov [esi + 4 * ebx], edi
		jmp end_partition
		
		*/
	}
	return sort_arr;
}

int main() {
	unsigned int size = 20;
	int* arr = randomArrayNumbers(size);

	printf_s("Random numbers size %d array:\n", size);
	print_arr(arr, size);
	printf_s("\n");

	// QuickSort C algorithm
	int* sorted_arr_C = quickSortC(arr, size);

	printf_s("C sorted array size %d:\n", size);
	print_arr(sorted_arr_C, size);
	printf_s("\n");

	delete sorted_arr_C;

	// QuickSort Assembler algorithm
//	print_arr(arr, size);
//	int* sorted_arr_ass = quickSortAssembler(arr, size);
//	print_arr(sorted_arr_ass, size);
}