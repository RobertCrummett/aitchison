#ifndef COMP_H
#define COMP_H

typedef struct comp_t comp_t;

struct comp_t {
	double* data;
	double* aux;
	size_t observations;
	size_t channels;
};

extern void comp_closure(comp_t* comp, double total);
extern void comp_perturbation(comp_t *first, comp_t *second, double total);
extern void comp_power(comp_t *comp, double exponent, double total);

extern void comp_centered_log_ratio(comp_t *comp);
extern void comp_inverse_centered_log_ratio(comp_t *comp, double total);

extern void comp_contrast_matrix(int *sequential_binary_parition, double *contrast_matrix, size_t channels);
extern void comp_isometric_log_ratio(comp_t *comp, double *contrast_matrix);
extern void comp_inverse_isometric_log_ratio(comp_t *comp, double *contrast_matrix, double total);

extern comp_t* comp_alloc(size_t observations, size_t channels);
extern void comp_free(comp_t** comp);

#endif
