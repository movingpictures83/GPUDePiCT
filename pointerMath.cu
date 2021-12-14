#include "pointerMath.h"

int* ptrMath2D(int* array, int row, int column, int rowlength){
	return (array+row*rowlength+column);
}
