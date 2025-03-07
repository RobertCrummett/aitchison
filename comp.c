#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "comp.h"

double comp_sum(double *data, size_t size, size_t stride) {
	double s = 0;
	for (size_t i = 0; i < size; i++)
		s += data[i * stride];
	return s;
}

void comp_scale(double* data, double a, size_t size, size_t stride) {
	for (size_t i = 0; i < size; i++)
		data[i * stride] *= a;
}

double comp_ddot(double *first, double *second, size_t size, size_t stride_f, size_t stride_s) {
	double s = 0;
	for (size_t i = 0; i < size; i++)
		s += first[i * stride_f] * second[i * stride_s];
	return s;
}

double comp_geometric_mean(double *data, size_t size, size_t stride)  {
	double s = 0;
	for (size_t i = 0; i < size; i++)
		s += log(data[i * stride]);
	return exp(s / size);
}

void comp_swap(double *first, double *second) {
	double temp = *first;
	*first = *second;
	*second = temp;
}

void comp_reverse(double *data, size_t size, size_t stride) {
	for (size_t i = 0; i < size / 2; i++)
		comp_swap(data + i * stride, data + size - 1 - i * stride);
}

void comp_roll(double *data, size_t size, size_t stride, size_t times) {
	comp_reverse(data, times, stride);
	comp_reverse(data + times, size - times, stride);
	comp_reverse(data, size, stride);
}

void comp_closure(comp_t* comp, double total) {
	size_t observations = comp->observations;
	size_t channels = comp->channels;
	double* data = comp->data;

	for (size_t i = 0; i < observations; i++) {
		double k = total / comp_sum(&data[i * channels], channels, 1);
		comp_scale(&data[i * channels], k, channels, 1);
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
	size_t channels = comp->channels;
	double *data = comp->data;

	for (size_t i = 0; i < observations; i++) {
		double geomean = comp_geometric_mean(&data[i * channels], channels, 1);
		double log_geomean = log(geomean);

		for (size_t j = 0; j < channels; j++)
			data[i * channels + j] = log(data[i * channels + j]) - log_geomean;
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
			aux[j] = comp_ddot(&data[i * channels], &contrast_matrix[j * channels], channels, 1, 1);

		for (size_t j = 0; j < channels - 1; j++)
			data[i * channels + j] = aux[j];
	}
}

void comp_inverse_isometric_log_ratio(comp_t *comp, double *contrast_matrix, double total) {
	double *data = comp->data;
	double *aux = comp->aux;
	size_t observations = comp->observations;
	size_t channels = comp->channels;

	for (size_t i = 0; i < observations; i++) {
		/*
		 * Cannot simply write into data[j + i * channels]
		 * This would overwrite during the final iteration.
		 * Instead write to aux[j] for simplicity
		 */
		for (size_t j = 0; j < channels; j++)
			aux[j] = comp_ddot(&data[i * channels], &contrast_matrix[j], channels - 1, 1, channels);

		for (size_t j = 0; j < channels; j++)
			data[i * channels + j] = aux[j];
	}

	comp_inverse_centered_log_ratio(comp, total);
}

comp_t* comp_alloc(size_t observations, size_t channels) {
	comp_t *comp = malloc(sizeof *comp);
	if (comp == NULL)
		return NULL;

	double *data = malloc(observations * channels * sizeof *data);
	if (data == NULL) {
		free(comp);
		return NULL;
	}

	double *aux = malloc(channels * sizeof *aux);
	if (aux == NULL) {
		free(comp);
		free(data);
		return NULL;
	}

	comp->data = data;
	comp->aux = aux;
	comp->observations = observations;
	comp->channels = channels;

	return comp;
}

void comp_free(comp_t** comp) {
	if (comp == NULL || *comp == NULL)
		return;

	if ((*comp)->data != NULL)
		free((*comp)->data);
	
	if ((*comp)->aux != NULL)
		free((*comp)->aux);

	free(*comp);
	*comp = NULL;
	return;
}
