/* ************************************************************

   implicit evolution of the spins

   ************************************************************ */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <time.h>


// the inputs are: noise, initial spin configuration, interaction matrix
// external field on each spin, demagnetizing tensor for each spin
// should be able to generate the detailed time evolution data

// ------------------------------------------------------------------------------------

struct rparams // all the parameters at a given time slice, different parameters allowed for each spin
{
    int n; // number of spins
    gsl_vector* m; // the spin configuration at the previous time
    gsl_vector* hext; // external field
    gsl_vector* Nv; // demagnetizing tensor
    gsl_vector* dw; // noise 
    //gsl_vector* K1; // uniaxial anisotropy, alreday included in Nv
    gsl_vector* K2; // biaxial anisotropy (effective rescaled value)
    gsl_matrix* Cp; // the interaction matrix
    double dt; //time step
    double alpha; 
    double nu; // nouse intensity
};    

// ------------------------------------------------------------------------------------
void cross3d1 (double v1x, double v1y, double v1z,
	      double v2x, double v2y, double v2z,
	      double* rx, double* ry, double* rz)
{
    *rx = v1y*v2z - v1z*v2y;
    *ry = v1z*v2x - v1x*v2z;
    *rz = v1x*v2y - v1y*v2x;
}

// --------------------------------------------------------------------------------------
// the Landau-Lifshitz-Gilbert equation
int llg_f (const gsl_vector * x, void *params, gsl_vector* f)
{
    int n = ((struct rparams *) params)->n;
    gsl_vector* m = ((struct rparams *) params)->m;
    gsl_vector* hext = ((struct rparams *) params)->hext;
    gsl_vector* Nv = ((struct rparams *) params)->Nv;
    gsl_vector* dw = ((struct rparams *) params)->dw;
    //gsl_vector* K1 = ((struct rparams *) params)->K1;
    gsl_vector* K2 = ((struct rparams *) params)->K2;
    gsl_matrix* Cp = ((struct rparams *) params)->Cp;
    double dt = ((struct rparams *) params)->dt;
    double dt2 = sqrt (dt);
    double alpha = ((struct rparams *) params)->alpha;
    double nu = ((struct rparams *) params)->nu;

    gsl_vector* mid = gsl_vector_alloc (n*3); // mid-point
    gsl_vector* df = gsl_vector_alloc (n*3); // difference
    gsl_vector* h = gsl_vector_alloc (n*3); // effective field
    gsl_vector* Cpm = gsl_vector_alloc (n*3); // product of Cp and mid-point

    int i, j, j1, j2, j3;
    double a0, a1, midi;
    double v1x, v1y, v1z, v2x, v2y, v2z;
    double mh[3], mf[3];
    double tmpx, tmpy, tmpz;

    // mid-point and difference
    for (i = 0; i < n*3; i++) {
	a0 = gsl_vector_get (m, i);
	a1 = gsl_vector_get (x, i);
	gsl_vector_set (mid, i, (a0+a1)/2.);
	gsl_vector_set (df, i, a1-a0);
    }

    // effective field
    gsl_blas_dgemv (CblasNoTrans, 1., Cp, mid, 0., Cpm); // Cpm = Cp %*% mid
    
    for (i = 0; i < n*3; i++) {
	j = (int) (i / 3.);
	midi = gsl_vector_get (mid, i);
	gsl_vector_set (h, i, gsl_vector_get(hext,i) - gsl_vector_get(Nv,i)*midi + gsl_vector_get(Cpm,i) - 2*gsl_vector_get(K2,j)*midi*(1-midi*midi));
    }

    for (i = 0; i < n; i++) 
    {
	j1 = i * 3; j2 = j1 + 1; j3 = j2 + 1;

	v1x = gsl_vector_get (mid, j1);	v1y = gsl_vector_get (mid, j2);	v1z = gsl_vector_get (mid, j3);

	v2x = gsl_vector_get (h, j1) * dt + gsl_vector_get (dw, j1) * nu * dt2; 
	v2y = gsl_vector_get (h, j2) * dt + gsl_vector_get (dw, j2) * nu * dt2; 
	v2z = gsl_vector_get (h, j3) * dt + gsl_vector_get (dw, j3) * nu * dt2; 
	cross3d1 (v1x, v1y, v1z, v2x, v2y, v2z, mh, mh+1, mh+2);

	v2x = gsl_vector_get (df, j1); v2y = gsl_vector_get (df, j2); v2z = gsl_vector_get (df, j3); 
	cross3d1 (v1x, v1y, v1z, v2x, v2y, v2z, mf, mf+1, mf+2);

	tmpx = -mh[0] + alpha*mf[0] - gsl_vector_get(df,j1);
	tmpy = -mh[1] + alpha*mf[1] - gsl_vector_get(df,j2);
	tmpz = -mh[2] + alpha*mf[2] - gsl_vector_get(df,j3);

	gsl_vector_set (f, j1, tmpx); gsl_vector_set (f, j2, tmpy); gsl_vector_set (f, j3, tmpz);
    }
	
    gsl_vector_free (mid); gsl_vector_free (df); gsl_vector_free (h); gsl_vector_free (Cpm);

    return GSL_SUCCESS;
}

// -------------------------------------------------------------------------------------
// one step
// m1 stores the new values
int solve1 (int n, gsl_vector* m, gsl_vector* hext, gsl_vector* Nv, gsl_vector* dw, gsl_vector* K2, gsl_matrix* Cp, double dt, double alpha, double nu, gsl_vector* m1, double error_threshold, int iter_threshold)
{
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    
    int status;
    int i, iter = 0;
    
    struct rparams p = {n, m, hext, Nv, dw, K2, Cp, dt, alpha, nu};
    gsl_multiroot_function f = {&llg_f, n*3, &p};
    
    gsl_vector *x = gsl_vector_alloc (n*3);

    gsl_vector_memcpy (x, m);
    
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n*3);
    gsl_multiroot_fsolver_set (s, &f, x);
     
    //print_state (iter, s);

    do 
    {
	iter++;
	status = gsl_multiroot_fsolver_iterate (s);
	if (status) break;   // check if solver is stuck 
	status = gsl_multiroot_test_residual (s->f, error_threshold);
    } 
    while (status == GSL_CONTINUE && iter < iter_threshold);

    if (status == GSL_SUCCESS) {
	gsl_vector_memcpy (m1, s->x);
    } else {
	printf ("Cannot approach the required accuracy!\n");
    }
     
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);

    return 0;
}

// ------------------------------------------------------------------------------------

// implicit method for time evolution
void implicit_evol (
    double* m0, // initial spin configuration
    int* n, // number of spins, length(m0)=3*n
    double* Cp0, // interaction matrix store in vector form
    double* hext0, // external field on each spin
    double* Nv0, // demagnetizing tensor
    double* dW0, // noise
    double* K20, // biaxial aisotropy 
    double* alpha, 
    double* nu,
    double* dt, // time step
    int* T, // step number
    int* noiseMethod, // use pre-generated noise (1) or generate it in C (0)
    int* noiseType, // type of noise, default 0 is Gaussian, only effective when *noiseMethod==0
    int* record, // whether record detailed process, default 1 -record, 0 - do not record
    int* recordStep, // record step length
    int* recordList, // or record at given time slices
    int* recordLength, // length of recordList
    double* error_threshold, // used to judge whether to stop the solver iteration
    int* iter_threshold, // again judge whether to stop the solver iteration
    double* process // record the detailed info
    )
{
    int i, j, ct, ct1, ct2, l = (*n) * 3;
    gsl_matrix* Cp = gsl_matrix_alloc (l,l);
    gsl_vector* m = gsl_vector_alloc (l);
    gsl_vector* m1 = gsl_vector_alloc (l);
    gsl_vector* hext = gsl_vector_alloc (l);
    gsl_vector* Nv = gsl_vector_alloc (l);
    gsl_vector* dw = gsl_vector_alloc (l);
    gsl_vector* K2 = gsl_vector_alloc (l);
	
    // random generator ****************************--------------------------
    const gsl_rng_type* T1;
    gsl_rng* r;
    time_t t1;
    (void) time (&t1);

    gsl_rng_env_setup ();
    T1 = gsl_rng_default;
    r = gsl_rng_alloc (T1);
    gsl_rng_set (r, t1);
    // ***************************************---------------------------------

    ct = 0;
    for (i = 0; i < l; i++) {
	for (j = 0; j < l; j++) {
	    gsl_matrix_set (Cp, j, i, Cp0[ct]);
	    ct++;
	}
	gsl_vector_set (m, i, m0[i]);
	gsl_vector_set (hext, i, hext0[i]);
	gsl_vector_set (Nv, i, Nv0[i]);
    }

    for (i = 0; i < *n; i++) gsl_vector_set (K2, i, K20[i]);

    ct = 0; ct1 = 0; ct2 = 0;
    for (i = 0; i < *T; i++) 
    {
	if (*noiseMethod == 1) {
	    for (j = 0; j < l; j++) {
		gsl_vector_set (dw, j, dW0[ct]);
		ct++;
	    }
	} else {
	    // ******************************----------------------------------
	    // generate noise in C according to noiseType
	    // for now implement only Gaussian
	    // more possibilities
	    switch (*noiseType) 
	    {
	    case 0:
		for (j = 0; j < l; j++) 
		    gsl_vector_set (dw, j, gsl_ran_gaussian (r, 1));
		break;
	    default:
		break;
	    }
	    // ******************************--------------------------------------
	}

	solve1 (*n, m, hext, Nv, dw, K2, Cp, *dt, *alpha, *nu, m1, *error_threshold, *iter_threshold);
	gsl_vector_memcpy (m, m1);

	if (((*recordStep != 0 && i % (*recordStep) == 0) || (*recordStep == 0 &&  ct2 < *recordLength && i == recordList[ct2])) && *record == 1) {
	    for (j = 0; j < l; j++) {
		process[ct1] = gsl_vector_get (m, j);
		ct1++;
	    }
	    ct2++;
	}
    }

    // random generator **********************************---------------
    gsl_rng_free (r);
    // ***************************************************---------------
    
    gsl_vector_free (hext);
    gsl_vector_free (Nv);
    gsl_vector_free (dw);
    gsl_vector_free (K2);
    gsl_vector_free (m1);
    gsl_vector_free (m);
    gsl_matrix_free (Cp);
}
