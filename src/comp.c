#include <assert.h>
#include <math.h>
#include <stddef.h>

#include "comp.h"

double comp_sum(double *data, size_t size, size_t stride) {
	double s = 0;
	for ( ; stride * size > 0; size -= stride)
		s += data[size - 1];
	return s;
}

void comp_scale(double* data, double a, size_t size, size_t stride) {
	for ( ; stride * size > 0; size -= stride)
		data[size - 1] *= a;
}

double comp_ddot(double *first, double *second, size_t size, size_t stride_f, size_t stride_s) {
	double s = 0;
	for ( ; size > 0; size--)
		s += first[(size - 1) * stride_f] * second[(size - 1) * stride_s];
	return s;
}

double comp_geometric_mean(double *data, size_t size, size_t stride)  {
	double s = 0;
	for (size_t i = 0; i < size * stride; i += stride)
		s += log(data[i * stride]);
	return exp(s / size);
}

void comp_swap(double *first, double *second) {
	double temp = *first;
	*first = *second;
	*second = temp;
}

void comp_roll(double *data, size_t size, size_t stride, size_t times) {
	for ( ; times > 0; times--)
		for (size_t i = size - 1; i > 0; i--)
			comp_swap(data + i * stride, data + (i - 1) * stride);
}

void comp_closure(comp_t* comp, double total) {
	size_t observations = comp->observations;
	size_t channels = comp->channels;
	double* data = comp->data;

	for (size_t i = 0; i < observations; i++) {
		double k = total / comp_sum(data, channels, 1);
		comp_scale(data, k, channels, 1);
		if (i < observations - 1)
			data += channels;
	}
}

void comp_perturbation(comp_t *first, comp_t *second, double total) {
	size_t data_size = first->observations * first->channels;
	double *data_first = first->data;
	double *data_second = second->data;

	for (size_t i = 0; i < data_size; i++)
		data_first[i] *= data_second[i];
	comp_closure(first, total);
}

void comp_power(comp_t *comp, double exponent, double total) {
	size_t data_size = comp->observations * comp->channels;
	double *data = comp->data;
	for (size_t i = 0; i < data_size; i++)
		data[i] = pow(data[i], exponent);
	comp_closure(comp, total);
}

void comp_centered_log_ratio(comp_t *comp) {
	size_t observations = comp->observations;
	size_t channels = comp->observations;
	double *data = comp->data;
	for ( ; observations > 0; observations--) {
		double geomean = comp_geometric_mean(data, channels, 1);
		double log_geomean = log(geomean);
		for (size_t i = 0; i < channels; i++)
			data[i] = log(data[i]) - log_geomean;
		if (observations > 1)
			data += channels;
	}
}

void comp_inverse_centered_log_ratio(comp_t *comp, double total) {
	size_t data_size = comp->observations * comp->channels;
	double *data = comp->data;
	for (size_t i = 0; i < data_size; i++)
		data[i] = exp(data[i]);
	comp_closure(comp, total);
}

void comp_contrast_matrix(int *sequential_binary_parition, double *contrast_matrix, size_t channels) {
	for (size_t i = 0; i < channels - 1; i++) {
		double r = 0, s = 0;

		for (size_t j = 0; j < channels; j++)
			switch (sequential_binary_parition[j]) {
				case 0:
					break;
				case 1:
					r++;
					break;
				case -1:
					s++; 
					break;
				default:
					assert(0 && "UNREACHABLE\n");
					break;
			}

		double psip = sqrt(s / (r * (r + s)));
		double psim =-sqrt(r / (s * (r + s)));

		for (size_t j = 0; j < channels; j++)
			switch (sequential_binary_parition[j]) {
				case 0:
					contrast_matrix[j] = 0;
					break;
				case 1:
					contrast_matrix[j] = psip;
					break;
				case -1:
					contrast_matrix[j] = psim;
					break;
				default:
					assert(0 && "UNREACHABLE\n");
					break;
			}

		if (i < channels - 2) {
			sequential_binary_parition += channels;
			contrast_matrix += channels;
		}
	}
}

void comp_isometric_log_ratio(comp_t *comp, double *contrast_matrix) {
	comp_centered_log_ratio(comp);

	double *data = comp->data;
	double *aux = comp->aux;
	size_t observations = comp->observations;
	size_t channels = comp->channels;

	for (size_t i = 0; i < observations; i++) {
		for (size_t j = 0; j < channels - 1; j++)
			aux[j] = comp_ddot(data, contrast_matrix + j * channels, channels, 1, 1);

		for (size_t j = 0; j < channels - 1; j++)
			(data - i)[j] = aux[j];

		if (i < observations - 1)
			data += channels;
	}
}

void comp_inverse_isometric_log_ratio(comp_t *comp, double *contrast_matrix, double total) {
	double *data = comp->data;
	double *aux = comp->aux;
	size_t observations = comp->observations;
	size_t channels = comp->channels;

	comp_roll(data, observations * channels, 1, observations);

	double* head = data + observations;
	for (size_t i = 0; i < observations; i++) {
		//
		// Cannot simply write into data[j + i * channels]
		// This would overwrite during the final iteration.
		// Instead write to aux[j] for simplicity
		//
		for (size_t j = 0; j < channels; j++)
			aux[j] = comp_ddot(head, contrast_matrix + j, channels - 1, 1, channels);
		for (size_t j = 0; j < channels; j++)
			data[j] = aux[j];

		if (i < observations - 1) {
			data += channels;
			head += channels - 1;
		}
	}

	comp_closure(comp, total);
}
