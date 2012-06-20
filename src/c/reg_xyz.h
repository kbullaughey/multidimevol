#ifndef __REG_XYZ_H__
#define __REG_XYZ_H__

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#define EVOLVABLE_PARAMS 5

extern "C" {
   int c_evolution_matrix (
      double *m_in, int *m_len, double *pm, int *pm_len, int *pm_rows_out,
      double *max_duration, double *time_resolution, double *tolerance, 
      double *max_step, int *error, int *evolvable_param_bitmask);
   int c_fitness(double *m_in, int *m_len, int *error, double *fit);
   int c_gradient_field(double *m_in, int *m_len, double *param_matrix, 
      int *pm_len, double *vector_field, int *vf_len, int *error);
}

int c_evolution_matrix (
   double *m_in, int *m_len, double *pm, int *pm_len, int *pm_rows_out,
   double *max_duration, double *time_resolution, double *tolerance, 
   double *max_step, int *error, int *evolvable_param_bitmask);
int c_fitness(double *m_in, int *m_len, int *error, double *fit);
int c_gradient_field(double *m_in, int *m_len, double *param_matrix, 
   int *pm_len, double *vector_field, int *vf_len, int *error);


#endif /* __REG_XYZ_H__ */
