#include <stdio.h>
#include <stdlib.h>

#include "comp.h"

int main(void) {
	double* data = malloc(18 * sizeof(double));
	for (size_t i = 0; i < 18; i++)
		data[i] = i;

	for (size_t i = 0; i < 18; i++) {
		printf("%lf ", data[i]);
		if (i % 6 == 5)
			printf("\n");
	}
	printf("\n");

	free(data);

	return 0;
}

