/* ************************************************************

   explicit evolution of the spins

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
#include <stdio.h>
#include <time.h>


// the inputs are: noise, initial spin configuration, interaction matrix
// external field on each spin, demagnetizing tensor for each spin
// should be able to generate the detailed time evolution data

// ------------------------------------------------------------------------------------
void cross3d (double v1x, double v1y, double v1z,
	      double v2x, double v2y, double v2z,
	      double* rx, double* ry, double* rz)
{
    *rx = v1y*v2z - v1z*v2y;
    *ry = v1z*v2x - v1x*v2z;
    *rz = v1x*v2y - v1y*v2x;
}

// -------------------------------------------------------------------------------------
// one step
// m1 stores the new values
int solve (int n, gsl_vector* m, 
	   int ne, double* me, // fixed external magnets
	   double* hext, gsl_vector* Nv, 
	   gsl_vector* dw, gsl_vector* K2, 
	   gsl_matrix* Cp, 
	   double dt, double alpha, 
	   double nu, gsl_vector* m1)
{
    int i, j, j1, j2, j3;
    double dt2 = sqrt (dt), alpha2 = 1 + alpha*alpha;
    double v1x, v1y, v1z, v2x, v2y, v2z, H1, H2, H3, mi;
    double mh[3], mH[3];
    int l1 = 3*n, l2 = 3*(n+ne);

    gsl_vector* x = gsl_vector_alloc (3);
    gsl_vector* h = gsl_vector_alloc (3*n); // effective field
    gsl_vector* Cpm = gsl_vector_alloc (l2); // product of Cp and mid-point
    gsl_matrix* a = gsl_matrix_alloc (3,3);
    gsl_vector* b = gsl_vector_alloc (3);
    gsl_vector* am = gsl_vector_alloc (l2);

    // need to solve n group eqautions, each has 3 coupled equations
    if (ne == 0) {
	gsl_blas_dgemv (CblasNoTrans, 1., Cp, m, 0., Cpm); // Cpm = Cp %*% m
    } else {
	for (i = 0; i < l1; i++) gsl_vector_set (am, i, gsl_vector_get(m,i));
	for (i = l1; i < l2; i++) gsl_vector_set (am, i, me[i-l1]);
	gsl_blas_dgemv (CblasNoTrans, 1., Cp, am, 0., Cpm); // Cpm = Cp %*% am
    }

    for (i = 0; i < l1; i++) {
	j = (int) (i / 3.);
	mi = gsl_vector_get (m, i);
	gsl_vector_set (h, i, (hext[i] - gsl_vector_get(Nv,i)*mi + gsl_vector_get(Cpm,i) - 2*gsl_vector_get(K2,j)*mi*(1-mi*mi))*dt + gsl_vector_get(dw,i)*nu*dt2);
    }
    // printf ("%f\n", gsl_vector_get(dw,0)*nu*dt2);/////////
    for (i = 0; i < n; i++)
    {
	j1 = i * 3; j2 = j1 + 1; j3 = j2 + 1;
	
	v1x = gsl_vector_get (m, j1); v1y = gsl_vector_get (m, j2); v1z = gsl_vector_get (m, j3);
	v2x = gsl_vector_get (h, j1); v2y = gsl_vector_get (h, j2); v2z = gsl_vector_get (h, j3);
	cross3d (v1x, v1y, v1z, v2x, v2y, v2z, mh, mh+1, mh+2);

	H1 = (gsl_vector_get(h,j1)+alpha*mh[0])/alpha2;
	H2 = (gsl_vector_get(h,j2)+alpha*mh[1])/alpha2;
	H3 = (gsl_vector_get(h,j3)+alpha*mh[2])/alpha2;
	cross3d (v1x, v1y, v1z, H1, H2, H3, mH, mH+1, mH+2);

	gsl_vector_set (b, 0, gsl_vector_get(m,j1)-0.5*mH[0]);
	gsl_vector_set (b, 1, gsl_vector_get(m,j2)-0.5*mH[1]);
	gsl_vector_set (b, 2, gsl_vector_get(m,j3)-0.5*mH[2]);

	gsl_matrix_set (a, 0, 0, 4+H1*H1); gsl_matrix_set (a, 0, 1, H1*H2-2*H3); gsl_matrix_set (a, 0, 2, 2*H2+H1*H3);
	gsl_matrix_set (a, 1, 0, H1*H2+2*H3); gsl_matrix_set (a, 1, 1, 4+H2*H2); gsl_matrix_set (a, 1, 2, -2*H1+H2*H3);
	gsl_matrix_set (a, 2, 0, -2*H2+H1*H3); gsl_matrix_set (a, 2, 1, 2*H1+H2*H3); gsl_matrix_set (a, 2, 2, 4+H3*H3);
	gsl_matrix_scale (a, 1./(4+H1*H1+H2*H2+H3*H3));

	gsl_blas_dgemv (CblasNoTrans, 1., a, b, 0., x);

	gsl_vector_set (m1, j1, gsl_vector_get (x,0));
	gsl_vector_set (m1, j2, gsl_vector_get (x,1));
	gsl_vector_set (m1, j3, gsl_vector_get (x,2));
    }
    
    gsl_vector_free (h);
    gsl_vector_free (Cpm);
    gsl_matrix_free (a);
    gsl_vector_free (b);
    gsl_vector_free (x);
    gsl_vector_free (am);

    return 0;
}

// ------------------------------------------------------------------------------------

// implicit method for time evolution
// external magnets included
// time changing field included
//
void explicit_evol (
    double* m0, // initial spin configuration
    int* n, // number of spins, length(m0)=3*n
    double* me, // external magnets configuration, fixed absolutely, if not (for example, trapped in energy trap), can be given together in m0
    int* ne, // number of external spins
    double* Cp0, // interaction matrix store in vector form
    double* hext0, // external field on each spin
    int* htype, // external field type, 0 for constant, 1 for time-dependent
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
    int i, j, ct, ct1, ct2, ct3;
    int l = (*n) * 3, cpl = ((*n) + (*ne)) * 3;
 
    gsl_matrix* Cp = gsl_matrix_alloc (cpl,cpl);
    gsl_vector* m = gsl_vector_alloc (l);
    gsl_vector* m1 = gsl_vector_alloc (l);
    //gsl_vector* hext = gsl_vector_alloc (l);
    double* hext;
    gsl_vector* Nv = gsl_vector_alloc (l);
    gsl_vector* dw = gsl_vector_alloc (l);
    gsl_vector* K2 = gsl_vector_alloc (*n);
    
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
    for (i = 0; i < cpl; i++) {
	for (j = 0; j < cpl; j++) {
	    gsl_matrix_set (Cp, j, i, Cp0[ct]);
	    ct++;
	}
    }

    for (i = 0; i < l; i++) {
	gsl_vector_set (m, i, m0[i]);	
	gsl_vector_set (Nv, i, Nv0[i]);
    }

    hext = hext0;

    for (i = 0; i < *n; i++) gsl_vector_set (K2, i, K20[i]);

    ct = 0; ct1 = 0; ct2 = 0; ct3 = 0;
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

	if (*htype == 1 && i > 0) hext = hext + l;
	//printf ("%f\n", hext[0]);/////////////
	//printf ("%f\n", gsl_vector_get (dw, 1));///////
	solve (*n, m, *ne, me, hext, Nv, dw, K2, Cp, *dt, *alpha, *nu, m1);
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
 
    //gsl_vector_free (hext);
    hext = NULL; 
    gsl_vector_free (Nv); 
    gsl_vector_free (dw); 
    gsl_vector_free (K2); 
    gsl_vector_free (m1); 
    gsl_vector_free (m); 
    gsl_matrix_free (Cp); 
}
