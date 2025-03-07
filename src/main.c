#include <stdio.h>
#include <stdlib.h>
#include "comp.h"

#define PPMTOPERCENT 1e-3
#define PERCENTTOPPM 1e3	

static char line[256];

int main(void) {
	const double total = 100;
	const size_t observations = 1908573;
	const size_t channels = 4;

	int *x = malloc(observations * sizeof *x); 
	if (x == NULL) {
		fprintf(stderr, "Failed to allocate x\n");
		return 1;
	}
	int *y = malloc(observations * sizeof *y);
	if (y == NULL) {
		fprintf(stderr, "Failed to allocate y\n");
		free(x);
		return 1;
	}
	comp_t *comp = comp_alloc(observations, channels);
	if (comp == NULL) {
		fprintf(stderr, "Failed to allocate composition\n");
		free(x);
		free(y);
		return 1;
	}
	FILE *fin = fopen("share/merged.xyz", "r");	
	if (fin == NULL) {
		fprintf(stderr, "Failed to open share/merged.xyz\n");
		comp_free(&comp);
		free(x);
		free(y);
		return 1;
	}
	
	size_t c = 0;

	while (fgets(line, sizeof(line), fin)) {
		int xt = 0, yt = 0;
		double k = 0, t = 0, u = 0;
		sscanf(line, "%d %d %lf %lf %lf", &xt, &yt, &k, &t, &u);

		t *= PPMTOPERCENT;
		u *= PPMTOPERCENT;

		double i = total - k - t - u;
		x[c] = xt;
		y[c] = yt;

		comp->data[c * channels + 0] = k;
		comp->data[c * channels + 1] = t;
		comp->data[c * channels + 2] = u;
		comp->data[c * channels + 3] = i;

		c++;
	}

	fclose(fin);

	//
	//  The sequential binary partition looks
	//  like this:
	//
	//  Structure:            Interp:
	//     K    T    U    I     +     -
	//     1    1    1   -1    KTU vs I
	//    -1    1    1    0     TU vs K
	//     0    1   -1    0      T vs U
	//
	int sbp[12] = {1, 1, 1, -1, -1, 1, 1, 0, 0, 1, -1, 0};
	double contrast_matrix[12] = {0};

	comp_contrast_matrix(sbp, contrast_matrix, channels);

	comp_isometric_log_ratio(comp, contrast_matrix);

	FILE *fout = fopen("share/processed.xyz", "w");
	if (fout == NULL) {
		fprintf(stderr, "Failed to open share/processed.xyz\n");
		comp_free(&comp);
		free(x);
		free(y);
		return 1;
	}

	for (size_t i = 0; i < observations; i++) {
		double c0 = comp->data[i * channels + 0];
		double c1 = comp->data[i * channels + 1];
		double c2 = comp->data[i * channels + 2];

		fprintf(fout, "%d %d %.15lf %.15lf %.15lf\n", x[i], y[i], c0, c1, c2);
	}
	
	fclose(fout);
	comp_free(&comp);
	free(x);
	free(y);

	return 0;
}

