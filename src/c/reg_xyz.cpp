#include "reg_xyz.h"
#include <stdio.h>

#define OUT_OF_BOUNDS -1000

typedef struct _gene_params {
   double a;
   double b;
   double n;
} gene_params_t;

typedef struct _model {
   gene_params_t Y;
   gene_params_t Z;
   double k;
   double c[2];
   double delta;
   double max_step;
   int ep_bitmask;
} model_t;

int dFdt(double f[], model_t &m);
int dFdt_ode (double t, const double evol_par[], double f[], void *mod);
void store_row(double *r, double t, double fit, double *y);
void model_setup(double *m_in, model_t &m);
double F (model_t &m);
double dFda_Y (model_t &m);
double dFdb_Y (model_t &m);
double dFda_Z (model_t &m);
double dFdb_Z (model_t &m);
double dFdk (model_t &m);

void store_row(double *r, double t, double fit, double *y) {
   if (r == 0 || y == 0) {
      fprintf(stderr, "C: store_row: null pointer\n");
      return;
   }
   r[0] = t;
   r[1] = fit;
   for (int k = 0; k < EVOLVABLE_PARAMS; k++) {
      r[k + 2] = y[k];
   }
   return;
}

double F (model_t &m) {
   double p1, p2, p3;
   if (m.k < m.Y.b/m.Y.a * (1-exp(-4*m.Y.a))) {
      p1 = m.Z.b*m.Z.n*pow(4*m.Y.a+log(1-m.Y.a/m.Y.b*m.k),2)/(8*pow(m.Y.a,2)*m.Z.a);
      p3 = (m.Z.b*m.c[1]*m.delta*(
        2 - 8*m.Z.a - 2*exp(-m.Z.a*(4+log(1-m.Y.a/m.Y.b*m.k)/m.Y.a)) - 
          (2*m.Z.a*log(1-m.Y.a/m.Y.b*m.k))/m.Y.a+pow(m.Z.a,2)*pow(4+log(1-m.Y.a/m.Y.b*m.k)/m.Y.a,2)
        )) / (8*pow(m.Z.a,3));
   } else {
      p1 = 0;
      p3 = 0;
   }
   if (m.k < m.Y.b/m.Y.a * (1-exp(-m.Y.a/2))) {
      p2 = (m.Z.b*m.Z.n*(m.Y.a-2*log(1-m.Y.a/m.Y.b*m.k))*pow(m.Y.a+2*log(1-m.Y.a/m.Y.b*m.k),3)) / 
        (4*pow(m.Y.a,4)*m.Z.a);
   } else {
      p2 = 0;
   }

   double fitness = - (m.Y.b*(m.c[0]+8*m.c[1])*m.Y.n) / (4*m.Y.a) - m.c[1]*p1 - m.c[0]*p2 + p3;
   double penalty1 = exp((m.Y.a + m.Y.b + m.Z.a + m.Z.b + m.k)/25);
   double penalty2 = 1/(m.Y.a * m.Y.b * m.Z.a * m.Z.b * m.k);
   
   return fitness - penalty1 - penalty2;
}

double dFda_Y (model_t &m) {
  double p1, p2, p3;
  double tg = log(1-(m.Y.a*m.k)/m.Y.b);
  if (m.k < m.Y.b/m.Y.a * (1-exp(-4*m.Y.a))) {
    p1 = -(m.Z.b*m.Z.n*(4*m.Y.a+tg)*(-m.Y.a*m.k+(-m.Y.b+m.Y.a*m.k)*tg)) / (4*pow(m.Y.a,3)*m.Z.a*(-m.Y.b+m.Y.a*m.k));
    p3 = -(1/(4*pow(m.Y.a,3)*pow(m.Z.a,2)*(-m.Y.b+m.Y.a*m.k))) * 
      exp(-m.Z.a*(4+tg/m.Y.a)) * m.Z.b * m.c[1] * m.delta * 
      (m.Y.a*(1+exp(m.Z.a*(4+tg/m.Y.a))*(-1+4*m.Z.a))+exp(m.Z.a*(4+tg/m.Y.a))*m.Z.a*tg) * 
      (-m.Y.a*m.k+(-m.Y.b+m.Y.a*m.k)*tg);
  } else {
    p1 = 0;
    p3 = 0;
  }

  if (m.k < m.Y.b/m.Y.a * (1-exp(-m.Y.a/2))) {
    p2 = -(m.Z.b*m.Z.n*(m.Y.a-4*tg)*pow(m.Y.a+2*tg,2)*(-m.Y.a*m.k+(-m.Y.b+m.Y.a*m.k)*tg)) / 
      (pow(m.Y.a,5)*m.Z.a*(-m.Y.b+m.Y.a*m.k));
  } else {
    p2 = 0;
  }

  double ret = (m.Y.b*m.c[0]*m.Y.n)/(4*pow(m.Y.a,2)) + (2*m.Y.b*m.c[1]*m.Y.n)/pow(m.Y.a,2) - m.c[1]*p1 - m.c[0]*p2 + p3;
  return ret;
}

double dFdb_Y (model_t &m) {
  double p1, p2, p3;
  double tg = log(1-(m.Y.a*m.k)/m.Y.b);
  if (m.k < m.Y.b/m.Y.a * (1-exp(-4*m.Y.a))) {
    p1 = (4*m.Y.a*m.Z.b*m.k*m.Z.n+m.Z.b*m.k*m.Z.n*tg) / (4*m.Y.a*m.Z.a*pow(m.Y.b,2)-4*pow(m.Y.a,2)*m.Z.a*m.Y.b*m.k);
    p3 = (m.Z.b*m.c[1]*m.delta*m.k*(-1+exp(-m.Z.a*(4+tg/m.Y.a))+4*m.Z.a+(m.Z.a*tg)/m.Y.a)) / (4*pow(m.Z.a,2)*m.Y.b*(m.Y.b-m.Y.a*m.k));
  } else {
    p1 = 0;
    p3 = 0;
  }

  if (m.k < m.Y.b/m.Y.a * (1-exp(-m.Y.a/2))) {
    p2 = -(m.Z.b*m.k*m.Z.n*(m.Y.a-4*tg)*pow(m.Y.a+2*tg,2)) / (pow(m.Y.a,3)*m.Z.a*m.Y.b*(-m.Y.b+m.Y.a*m.k));
  } else {
    p2 = 0;
  }

  double ret = -(m.c[0]*m.Y.n)/(4*m.Y.a) - (2*m.c[1]*m.Y.n)/m.Y.a - m.c[1]*p1 - m.c[0]*p2 + p3;
  return ret;
}

double dFda_Z (model_t &m) {
  double p1, p2, p3;
  double tg = log(1-(m.Y.a*m.k)/m.Y.b);
  if (m.k < m.Y.b/m.Y.a * (1-exp(-4*m.Y.a))) {
    p1 = -(m.Z.b*m.Z.n*pow(4*m.Y.a+tg,2)) / (8*pow(m.Y.a,2)*pow(m.Z.a,2));
    p3 = -(1/(8*pow(m.Y.a,2)*pow(m.Z.a,4))) * 
      exp(-m.Z.a*(4+tg/m.Y.a)) * m.Z.b * m.c[1] * m.delta * 
      (2*pow(m.Y.a,2)*(-3-4*m.Z.a+exp(m.Z.a*(4+tg/m.Y.a))*(3-8*m.Z.a+8*pow(m.Z.a,2))) + 
        2*m.Y.a*m.Z.a*(-1+exp(m.Z.a*(4+tg/m.Y.a))*(-2+4*m.Z.a))*tg + 
        exp(m.Z.a*(4+tg/m.Y.a))*pow(m.Z.a,2)*pow(tg,2));
  } else {
    p1 = 0;
    p3 = 0;
  }

  if (m.k < m.Y.b/m.Y.a * (1-exp(-m.Y.a/2))) {
    p2 = -(m.Z.b*m.Z.n*(m.Y.a-2*tg)*pow(m.Y.a+2*tg,3)) / (4*pow(m.Y.a,4)*pow(m.Z.a,2));
  } else {
    p2 = 0;
  }

  double ret = -m.c[1]*p1 - m.c[0]*p2 + p3;
  return ret;
}

double dFdb_Z (model_t &m) {
  double p1, p2, p3;
  double tg = log(1-(m.Y.a*m.k)/m.Y.b);
  if (m.k < m.Y.b/m.Y.a * (1-exp(-4*m.Y.a))) {
    p1 = (m.Z.n*pow(4*m.Y.a+tg,2)) / (8*pow(m.Y.a,2)*m.Z.a);
    p3 = (m.c[1]*m.delta*(2-2*exp(-m.Z.a*(4+tg/m.Y.a))-8*m.Z.a-(2*m.Z.a*tg)/m.Y.a+pow(m.Z.a,2)*pow(4+tg/m.Y.a,2)))/(8*pow(m.Z.a,3));
  } else {
    p1 = 0;
    p3 = 0;
  }

  if (m.k < m.Y.b/m.Y.a * (1-exp(-m.Y.a/2))) {
    p2 = (m.Z.n*(m.Y.a-2*tg)*pow(m.Y.a+2*tg,3)) / (4*pow(m.Y.a,4)*m.Z.a);
  } else {
    p2 = 0;
  }

  double ret = -m.c[1]*p1 - m.c[0]*p2 + p3;
  return ret;
}

double dFdk (model_t &m) {
  double p1, p2, p3;
  double tg = log(1-(m.Y.a*m.k)/m.Y.b);
  if (m.k < m.Y.b/m.Y.a * (1-exp(-4*m.Y.a))) {
    p1 = (m.Z.b*m.Z.n*(4*m.Y.a+tg)) / (4*m.Y.a*m.Z.a*(-m.Y.b+m.Y.a*m.k));
    p3 = (m.Z.b*m.c[1]*m.delta*(1-exp(-m.Z.a*(4+tg/m.Y.a))-m.Z.a*(4+tg/m.Y.a))) / (4*pow(m.Z.a,2)*(m.Y.b-m.Y.a*m.k));
  } else {
    p1 = 0;
    p3 = 0;
  }

  if (m.k < m.Y.b/m.Y.a * (1-exp(-m.Y.a/2))) {
    p2 = (m.Z.b*m.Z.n*(m.Y.a-4*tg)*pow(m.Y.a+2*tg,2)) / (pow(m.Y.a,3)*m.Z.a*(-m.Y.b+m.Y.a*m.k));
  } else {
    p2 = 0;
  }

  double ret = -m.c[1]*p1 - m.c[0]*p2 + p3;
  return ret;
}

int dFdt_ode (double t, const double evol_par[], double f[], void *mod) {
   model_t *m = (model_t *)mod;

   /* use the current parameters with the ones that have evolved */
   m->Y.a = evol_par[0];
   m->Y.b = evol_par[1];
   m->Z.a = evol_par[2];
   m->Z.b = evol_par[3];
   m->k = evol_par[4];

   if (dFdt(f, *m) != 0) {
      fprintf(stderr, "C: dFdt_ode: dFdt failed\n");
      return 2;
   }

   if (m->max_step > 0) {
      /* I compute the fold-change in any dimension: c = (x2-x1)/x1. 
       * I only allow a maximum fold change, t = m->max_step.
       * So if t/c < 1 then I rescale the move vector by s = t/c
       */
      double s = 1;
      double s_try, c;
      double t = m->max_step;
      for (int k = 0; k < 5; k++) {
         c = fabs(f[k])/evol_par[k];
         s_try = t/c;
         if (s_try < s) s = s_try;
      }
      if (s < 1) {
        for (int k = 0; k < 5; k++) 
           f[k] = s * f[k];
      }
   }

   /* loop through and zero out non-evolvable parameters */
   for (int k = 0; k < 5; k++) {
      if (((1 << k) & m->ep_bitmask) == 0) {
         f[k] = 0;
      }
   }

   return GSL_SUCCESS;
}

int dFdt(double f[], model_t &m) {
   double penalty1 = exp((m.Y.a + m.Y.b + m.Z.a + m.Z.b + m.k)/25)/25;
   double penalty2 = 1/(m.Y.a * m.Y.b * m.Z.a * m.Z.b * m.k);
   
   /* compute the new gradient */
   f[0] = dFda_Y(m) - penalty1 + penalty2/m.Y.a;
   f[1] = dFdb_Y(m) - penalty1 + penalty2/m.Y.b;
   f[2] = dFda_Z(m) - penalty1 + penalty2/m.Z.a;
   f[3] = dFdb_Z(m) - penalty1 + penalty2/m.Z.b;
   f[4] = dFdk(m) - penalty1 + penalty2/m.k;

   return GSL_SUCCESS;
}

int c_gradient_field(double *m_in, int *m_len, double *param_matrix, 
      int *pm_len, double *vector_field, int *vf_len, int *error) {
   if (m_in == 0 || m_len == 0 || param_matrix == 0 || pm_len == 0 || 
         vector_field == 0 || vf_len == 0 || error == 0) {
      fprintf(stderr, "C: c_gradient_field: one or more pointers are null\n");
      return 1;
   }
   *error = 0;

   if (*pm_len != *vf_len || (*pm_len) % EVOLVABLE_PARAMS != 0) {
      fprintf(stderr, "C: c_gradient_field: length mismatch. pm_len=%d, vf_len=%d\n",
         *pm_len, *vf_len);
      *error = 1;
      return 0;
   }

   if (*m_len != 10) {
      *error = 2;
      fprintf(stderr, "C: c_gradient_field: m_len: %d\n", *m_len);
      return 0;
   }

   /* set up model parameters */
   model_t m;
   model_setup(m_in, m);

   int rows = (*pm_len) / EVOLVABLE_PARAMS;
   for (int i = 0; i < rows; i++) {
      m.Y.a = param_matrix[i*EVOLVABLE_PARAMS + 0];
      m.Y.b = param_matrix[i*EVOLVABLE_PARAMS + 1];
      m.Z.a = param_matrix[i*EVOLVABLE_PARAMS + 2];
      m.Z.b = param_matrix[i*EVOLVABLE_PARAMS + 3];
      m.k = param_matrix[i*EVOLVABLE_PARAMS + 4];
      if (dFdt(vector_field + i*EVOLVABLE_PARAMS, m) != 0) {
         fprintf(stderr, "C: gradient_field: dFdt failed\n");
         return 3;
      }
   }

   return 0;
}

/* fill in the model parameters, assumes m_in is correct length */
void model_setup(double *m_in, model_t &m) {
   m.Y.a = m_in[0];
   m.Y.b = m_in[1];
   m.Y.n = m_in[2];
   m.Z.a = m_in[3];
   m.Z.b = m_in[4];
   m.Z.n = m_in[5];
   m.k = m_in[6];
   m.c[0] = m_in[7];
   m.c[1] = m_in[8];
   m.delta = m_in[9];
   return;
}

int c_fitness(double *m_in, int *m_len, int *error, double *fit) {
   if (m_in == 0 || m_len == 0 || error == 0 || fit == 0) {
      fprintf(stderr, "C: c_fitness: one or more pointers are null\n");
      return 1;
   }
   *error = 0;

   if (*m_len != 10) {
      *error = 1;
      fprintf(stderr, "C: c_fitness: m_len: %d\n", *m_len);
      return 0;
   }

   /* set up model parameters */
   model_t m;
   model_setup(m_in, m);

   *fit = F(m);

   return 0;
}
     
/* this is the entry point from R. It computes the parameter evolution matrix with each column giving a path over time of a parameter */
int c_evolution_matrix(
      double *m_in, int *m_len,
      double *pm, int *pm_len, 
      int *pm_rows_out,
      double *max_duration, double *time_resolution, double *tolerance, double *max_step, 
      int *error, int *evolvable_param_bitmask) {

   if (m_in == 0 || m_len == 0 || pm == 0 || pm_len == 0 || pm_rows_out == 0 ||
         max_duration == 0 || time_resolution == 0 || tolerance == 0 || max_step == 0 || error == 0) {
      fprintf(stderr, "C: c_evolution_matrix: one or more pointers are null\n");
      fprintf(stderr, "C: c_evolution_matrix: %p, %p, %p, %p, %p, %p, %p, %p, %p, %p\n",
         m_in, m_len, pm, pm_len, pm_rows_out, max_duration, time_resolution, tolerance, max_step, error); 
      return 1;
   }

   *error = 0;

   if (*m_len != 10) {
      *error = 1;
      fprintf(stderr, "C: c_evolution_matrix: m_len: %d\n", *m_len);
      return 0;
   }
   /* make sure there's enough space to store the matrix */
   if ((*pm_len) < (2+EVOLVABLE_PARAMS) * (int)ceil((*max_duration) / (*time_resolution) + 2.0)) {
      fprintf(stderr, "C: c_evolution_matrix: pm_len=%d, max_duration=%.4e, time_resolution=%.4e, compare=%d\n",
         *pm_len, *max_duration, *time_resolution, 
         (2+EVOLVABLE_PARAMS) * (int)ceil((*max_duration) / (*time_resolution) + 1.0));
      *error = 2;
      return 0;
   }

   const gsl_odeiv_step_type *solver_type = gsl_odeiv_step_rk4;
   gsl_odeiv_step *s = gsl_odeiv_step_alloc (solver_type, EVOLVABLE_PARAMS);
 
   /* set up model parameters */
   model_t m;
   m.Y.a = m_in[0];
   m.Y.b = m_in[1];
   m.Y.n = m_in[2];
   m.Z.a = m_in[3];
   m.Z.b = m_in[4];
   m.Z.n = m_in[5];
   m.k = m_in[6];
   m.c[0] = m_in[7];
   m.c[1] = m_in[8];
   m.delta = m_in[9];
   m.max_step = *max_step;
   m.ep_bitmask = *evolvable_param_bitmask;

   gsl_odeiv_system ode_sys = {dFdt_ode, NULL, EVOLVABLE_PARAMS, &m};
 
   double t = 0.0, t_end = *max_duration;
   double y[EVOLVABLE_PARAMS], y_err[EVOLVABLE_PARAMS];
   double dfdt_in[EVOLVABLE_PARAMS], dfdt_out[EVOLVABLE_PARAMS];

   /* copy start position into y */
   y[0] = m.Y.a;
   y[1] = m.Y.b;
   y[2] = m.Z.a;
   y[3] = m.Z.b;
   y[4] = m.k;
 
   /* initialize dfdt_in from system parameters */
   GSL_ODEIV_FN_EVAL(&ode_sys, t, y, dfdt_in);
   *pm_rows_out = 1;

   /* store the initial record */
   store_row(pm, t, F(m), y);

   while (t < t_end) {
//      printf("df/dt=(%.4f, %.4f, %.4f, %.4f, %.4f)\n", dfdt_in[0], dfdt_in[1], dfdt_in[2], dfdt_in[3], dfdt_in[4]);
      *error = gsl_odeiv_step_apply(s, t, *time_resolution, y, y_err, dfdt_in, dfdt_out, &ode_sys);

      if (*error != GSL_SUCCESS) break;

      for (int k = 0; k < EVOLVABLE_PARAMS; k++) 
         dfdt_in[k] = dfdt_out[k];

      t += *time_resolution;

      /* make sure we haven't gone out of the meaninful range */
      if (y[0] < 0 || y[1] < 0 || y[2] < 0 || y[3] < 0 || y[4] < 0 < y[4] > 1) {
         fprintf(stderr, "C: out of bounds: y=(%.4f,%.4f,%.4f,%.4f,%.4f)\n",
            y[0], y[1], y[2], y[3], y[4]);
         break;
      }

      /* store the output in the matrix */
      store_row(pm + (*pm_rows_out) * (EVOLVABLE_PARAMS+2), t, F(m), y);
      (*pm_rows_out)++;

      /* if we've converged, then stop */
      if (*pm_rows_out > 20) {
         double sumsq = 0;
         for (int k = 0; k < EVOLVABLE_PARAMS; k++) 
            sumsq += pow(pm[(*pm_rows_out - 20) * (EVOLVABLE_PARAMS+2) + 2 + k] - y[k], 2);
         if (sumsq < (*tolerance) * (*time_resolution)) break;
      }
   }

   gsl_odeiv_step_free(s);
   return 0;
}

/* END */
