#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * data format:
 * observation 1: (channel 1, channel 2, ... , channel m)
 * observation 2: (channel 1, channel 2, ... , channel m)
 * observation 3: (channel 1, channel 2, ... , channel m)
 *       :
 * observation n: (channel 1, channel 2, ... , channel m)
 *
 * m -> # of channels
 * n -> # of observations
 *
 */

double sum(double *data, size_t size) {
	double s = 0;
	for ( ; size > 0; size--)
		s += data[size - 1];
	return s;
}

void scale(double* data, double a, size_t size) {
	for ( ; size > 0; size--)
		data[size - 1] *= a;
}

void closure(double *data, double total, size_t observations, size_t channels) {
	for ( ; observations > 0; observations--) {
		double k = total / sum(data, channels);
		scale(data, k, channels);
		if (observations > 1)
			data += channels;
	}
}

double geometric_mean(double *data, size_t size)  {
	double s = 0;
	for (size_t i = 0; i < size; i++)
		s += log(data[i]);
	return exp(s / size);
}

void perturbation(double *first, double *second, double total, size_t observations, size_t channels) {
	for (size_t i = 0; i < observations * channels; i++)
		first[i] *= second[i];
	return closure(first, total, observations, channels);
}

void power(double *data, double a, double total, size_t observations, size_t channels) {
	for (size_t i = 0; i < observations * channels; i++)
		data[i] = pow(data[i], a);
	return closure(data, total, observations, channels);
}

void centered_log_ratio(double *data, size_t observations, size_t channels) {
	for ( ; observations > 0; observations--) {
		double g = geometric_mean(data, channels);
		double log_g = log(g);
		for (size_t i = 0; i < channels; i++)
			data[i] = log(data[i]) - log_g;
		if (observations > 1)
			data += channels;
	}
}

void inverse_centered_log_ratio(double *data, double total, size_t observations, size_t channels) {
	for (size_t i = 0; i < observations * channels; i++)
		data[i] = exp(data[i]);
	closure(data, total, observations, channels);
}

void ait_contrast_matrix(int *sequential_binary_parition, double *contrast_matrix, size_t channels) {
	for (size_t i = 0; i < channels - 1; i++) {
		double r = 0, s = 0;
		for (size_t j = 0; j < channels; j++) {
			switch (sequential_binary_parition[j]) {
				case 0:  break;
				case 1:  r++; break;
				case -1: s++; break;
				default: assert(0 && "UNREACHABLE\n"); break;
			}
		}

		double psip = sqrt(s / (r * (r + s)));
		double psim =-sqrt(r / (s * (r + s)));

		for (size_t j = 0; j < channels; j++) {
			switch (sequential_binary_parition[j]) {
				case 0:  contrast_matrix[j] = 0;    break;
				case 1:  contrast_matrix[j] = psip; break;
				case -1: contrast_matrix[j] = psim; break;
				default: assert(0 && "UNREACHABLE\n"); break;
			}
		}

		if (i < channels - 2) {
			sequential_binary_parition += channels;
			contrast_matrix += channels;
		}
	}
}

double ddot(double *first, double *second, size_t size, size_t stride_f, size_t stride_s) {
	double s = 0;
	for ( ; size > 0; size--)
		s += first[(size - 1) * stride_f] * second[(size - 1) * stride_s];
	return s;
}

/*
 * aux size is channels - 1
 */
void isometric_log_ratio(double *data, double *contrast_matrix, double *aux, size_t observations, size_t channels) {
	centered_log_ratio(data, observations, channels);

	double *data_ptr = data;
	for (size_t i = 0; i < observations; i++) {
		for (size_t j = 0; j < channels - 1; j++)
			aux[j] = ddot(data_ptr, contrast_matrix + j * channels, channels, 1, 1);

		for (size_t j = 0; j < channels - 1; j++)
			(data_ptr - i)[j] = aux[j];

		if (i < observations - 1)
			data_ptr += channels;
	}
}

void swap(double *first, double *second) {
	double temp = *first;
	*first = *second;
	*second = temp;
}

void roll(double *data, size_t size, size_t steps) {
	for (size_t i = 0; i < size - steps % size; i++) {
		printf("%lf %lf\n", data[size - i - 1], data[size - i - 1 - steps % size]);
		swap(&data[size - i - 1], &data[size - i - 1 - steps % size]);
	}
	printf("\n");
}

void inverse_isometric_log_ratio(double *data, double *contrast_matrix, double *aux, double total, size_t observations, size_t channels) {
	roll(data, observations * channels, observations * (channels - 1));

	double *data_ptr = data + observations;
	for (size_t i = 0; i < observations; i++) {
		for (size_t j = 0; j < channels; j++)
			aux[j] = ddot(data_ptr, contrast_matrix + j, channels - 1, 1, channels);

		for (size_t j = 0; j < channels; j++)
			data[j] = aux[j];

		if (i < observations - 1) {
			data += channels;
			data_ptr += channels - 1;
		}
	}
	data -= (observations - 1) * channels; // reset the pointer

	closure(data, total, observations, channels);
}

int main(void) {

	size_t observations = 6, channels = 3;
	double* data = malloc(18 * sizeof(double));
	for (size_t i = 0; i < 18; i++)
		data[i] = (double)i;

	for (size_t i = 0; i < 18; i++) {
		printf("%lf ", data[i]);
		if (i % observations == 5) printf("\n");
	}
	printf("\n");

	roll(data, observations * channels, 1);
	// for (size_t i = 0; i < observations * (channels - 1); i++)
	// swap(&data[observations * channels - i - 1], &data[observations * (channels - 1) - i - 1]);

	for (size_t i = 0; i < 18; i++) {
		printf("%lf ", data[i]);
		if (i % observations == 5) printf("\n");
	}
	printf("\n");

	return 0;
}
