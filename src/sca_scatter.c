#include "sca_scatter.h"
#include "sca_workqueue.h" /* Queue_Prepare */

/*****************************************************************************/

void jur_sca_bascoord(double *dz,
	      double *dy,
	      double *ex,
	      double *ey,
	      double *ez) {
  
  double dotp, norm;
  
  int i;
  
  /* Set first orthonormal vector... */
  norm=NORM(dz);
  for(i=0; i<3; i++)
    ez[i]=dz[i]/norm;
  
  /* Normalize second direction... */
  norm=NORM(dy);
  for(i=0; i<3; i++)
    ey[i]=dy[i]/norm;
  
  /* Avoid linear dependences... */
  dotp=DOTP(ey, ez);
  if(dotp==1) {
    ey[0]=1;
    ey[1]=0;
    ey[2]=0;
    dotp=ez[0];
    if(dotp==1) {
      ey[0]=0;
      ey[1]=1;
      dotp=ez[1];
    }
  }
  
  /* Set second orthonormal vector (Gram-Schmidt)... */
  for(i=0; i<3; i++)
    ey[i]-=dotp*ez[i];
  
  norm=NORM(ey);
  for(i=0; i<3; i++)
    ey[i]/=norm;
  
  /* Get third orthonormal vector (cross-product)... */
  ex[0]=ey[1]*ez[2]-ey[2]*ez[1];
  ex[1]=ey[2]*ez[0]-ey[0]*ez[2];
  ex[2]=ey[0]*ez[1]-ey[1]*ez[0];
}

/*****************************************************************************/

void jur_sca_bhmie(double x,
	   double n_real,
	   double n_imag,
	   double *phase,
	   double *qext,
	   double *qsca) {
  
  gsl_complex cxan, cxbn, cxxi, cxy, cxxi1, cxtemp, cxref, 
    cxs1[NTHETAMAX+2], cxs2[NTHETAMAX+2], cxd[10000];
  
  double apsi, apsi1, chi, chi0, chi1, dang, fn, p, rn, t, theta2,
    xstop, ymod, dn, dx, psi, psi0, psi1, amu[NTHETAMAX+2], pi[NTHETAMAX+2],
    pi0[NTHETAMAX+2], pi1[NTHETAMAX+2], tau[NTHETAMAX+2];
  
  int i, j, jj, n, nmx, nn, nstop, ntheta;
  
  /* Set scattering angles, ntheta=(NTHETAMAX+1)/2... */
  if((NTHETAMAX+1)%2!=0)
    ERRMSG("NTHETAMAX needs to be odd!");
  ntheta=(NTHETAMAX+1)/2;
  
  /* Bohren-Huffman Mie code... */
  cxref=gsl_complex_rect(n_real, n_imag);
  dx=x;
  cxy=gsl_complex_mul(gsl_complex_rect(x, 0.0), cxref);
  
  /* Series expansion terminated after NSTOP terms */  
  xstop=x+4.E0*pow(x,0.3333)+2.0;
  nstop=(int)xstop;
  ymod = gsl_complex_abs(cxy);
  nmx=(int)((xstop>ymod ? xstop : ymod)+15);
  if(nmx>10000)
    ERRMSG("Too many Mie terms!");
  dang=.5E0*M_PI/(double)(ntheta-1);
  for(j=1; j<=ntheta; j++) {
    theta2=(double)(j-1)*dang;
    amu[j]=cos(theta2);
  }

  /* Logarithmic derivative D(J) calculated by downward recurrence */
  /*  beginning with initial value (0.,0.) at J=NMX */ 
  cxd[nmx] = gsl_complex_rect(0.E0, 0.E0);
  nn=nmx-1;
  for(n=1; n<= nn; n++) {
    rn=nmx-n+1;
    cxtemp=gsl_complex_add(cxd[nmx-n+1],
			   gsl_complex_div(gsl_complex_rect(rn, 0.0), cxy));
    cxtemp=gsl_complex_div(gsl_complex_rect(1.0, 0.0), cxtemp);
    cxd[nmx-n]
      =gsl_complex_sub(gsl_complex_div(gsl_complex_rect(rn, 0.0), cxy), cxtemp);
  }
  
  for(j=1; j<=ntheta; j++) {
    pi0[j]=0.E0;
    pi1[j]=1.E0;
  }
  
  nn=2*ntheta-1;
  for(j=1; j<=nn; j++) {
    cxs1[j] = gsl_complex_rect(0.E0, 0.E0);
    cxs2[j] = gsl_complex_rect(0.E0, 0.E0);
  }

  /* Riccati-Bessel functions with real argument X */
  /*  calculated by upward recurrence */
  psi0=cos(dx);
  psi1=sin(dx);
  chi0=-sin(x);
  chi1=cos(x);
  apsi1=psi1;
  cxxi1 = gsl_complex_rect(apsi1,-chi1);
  *qsca=0.E0;

  for(n=1; n<=nstop; n++) {  
    dn=n;
    rn=n;
    fn=(2.E0*rn+1.E0)/(rn*(rn+1.E0));
    psi=(2.E0*dn-1.E0)*psi1/dx-psi0;
    apsi=psi;
    chi=(2.E0*rn-1.E0)*chi1/x-chi0;
    cxxi=gsl_complex_rect(apsi,-chi);

    cxan=gsl_complex_div(cxd[n], cxref);
    cxan=gsl_complex_add(cxan, gsl_complex_rect(rn/x, 0.0));
    cxan=gsl_complex_mul(cxan, gsl_complex_rect(apsi, 0.0));
    cxan=gsl_complex_sub(cxan, gsl_complex_rect(apsi1, 0.0));

    cxtemp=gsl_complex_div(cxd[n], cxref);
    cxtemp=gsl_complex_add(cxtemp, gsl_complex_rect(rn/x, 0.0));
    cxtemp=gsl_complex_mul(cxtemp, cxxi);
    cxtemp=gsl_complex_sub(cxtemp, cxxi1);
    cxan=gsl_complex_div(cxan, cxtemp);

    cxbn=gsl_complex_mul(cxref, cxd[n]);
    cxbn=gsl_complex_add(cxbn, gsl_complex_rect(rn/x, 0.0));
    cxbn=gsl_complex_mul(cxbn, gsl_complex_rect(apsi, 0.0));
    cxbn=gsl_complex_sub(cxbn, gsl_complex_rect(apsi1, 0.0));

    cxtemp=gsl_complex_mul(cxref, cxd[n]);
    cxtemp=gsl_complex_add(cxtemp, gsl_complex_rect(rn/x, 0.0));
    cxtemp=gsl_complex_mul(cxtemp, cxxi);
    cxtemp=gsl_complex_sub(cxtemp, cxxi1);
    cxbn=gsl_complex_div(cxbn, cxtemp);
 
    *qsca=*qsca+(2.*rn+1.)
      *(gsl_complex_abs(cxan)*gsl_complex_abs(cxan)
	+gsl_complex_abs(cxbn)*gsl_complex_abs(cxbn));
    for(j=1; j<=ntheta; j++) {
      jj=2*ntheta-j;
      pi[j]=pi1[j];
      tau[j]=rn*amu[j]*pi[j]-(rn+1.E0)*pi0[j];
      p=pow(-1.0,n-1);
      
      cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(pi[j], 0.0));
      cxtemp
	=gsl_complex_add(cxtemp, gsl_complex_mul(cxbn, gsl_complex_rect(tau[j],
									0.0)));
      cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
      cxs1[j]=gsl_complex_add(cxs1[j], cxtemp);
      
      t=pow(-1.0,n);
      
      cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(tau[j], 0.0));
      cxtemp
	=gsl_complex_add(cxtemp, gsl_complex_mul(cxbn, gsl_complex_rect(pi[j],
									0.0)));
      cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
      cxs2[j]=gsl_complex_add(cxs2[j], cxtemp);
      
      if(j!=jj) {
	cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(pi[j]*p, 0.0));
	cxtemp
	  =gsl_complex_add(cxtemp,
			   gsl_complex_mul(cxbn,
					   gsl_complex_rect(tau[j]*t, 0.0)));
	cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
	cxs1[jj]=gsl_complex_add(cxs1[jj], cxtemp);
	cxtemp=gsl_complex_mul(cxan, gsl_complex_rect(tau[j]*t, 0.0));
	cxtemp
	  =gsl_complex_add(cxtemp,
			   gsl_complex_mul(cxbn,
					   gsl_complex_rect(pi[j]*p, 0.0)));
	cxtemp=gsl_complex_mul(gsl_complex_rect(fn, 0.0), cxtemp);
	cxs2[jj]=gsl_complex_add(cxs2[jj], cxtemp);
      }
    }
    
    psi0=psi1;
    psi1=psi;
    apsi1=psi1;
    chi0=chi1;
    chi1=chi;
    cxxi1=gsl_complex_rect(apsi1, -chi1);
    
    /* For each angle J, compute pi_n+1 from PI = pi_n , PI0 = pi_n-1 */
    for(j=1; j<=ntheta; j++) {
      pi1[j]=((2.*rn+1.)*amu[j]*pi[j]-(rn+1.)*pi0[j])/rn;
      pi0[j]=pi[j];
    }
  }
  
  /* Compute efficiencies... */
  *qsca=(2.E0/(x*x))*(*qsca);
  *qext=(4.E0/(x*x))*cxs1[1].dat[0];
  
  /* Compute phase function from scattering amplitudes... */
  /* calculate phase function following Liou: p192 eq 5.2.111a */
  /* P_11 = 4*PI*(i_1+i_2)/(2*k^2*sigma_s)) */
  /* i_1(theta), i_2(theta) = abs(S_1(theta))^2, abs(S_2(theta))^2 */
  /* intensity functions for perpendicular and parallel components */
  /* Q_s = sigma_s/(PI *rad^2) - scattering efficiency (here qsca) */
  /* sigma_s = scattering cross section */
  /* x = k*rad - size parameter (x=2*PI*rad/lambda) */

  for(i=0; i<2*ntheta-1; i++)
    phase[i]=2/(x*x*(*qsca))
      *(cxs1[i+1].dat[0]*cxs1[i+1].dat[0]+cxs1[i+1].dat[1]*cxs1[i+1].dat[1]
	+cxs2[i+1].dat[0]*cxs2[i+1].dat[0]+cxs2[i+1].dat[1]*cxs2[i+1].dat[1]);
}

/*****************************************************************************/
void jur_sca_copy_aero(ctl_t *ctl,
	       aero_t *aero_dest,
	       aero_t *aero_src,
	       int init) {

  int id, ia, im, il;
  
  /* Copy data... */
  memcpy(aero_dest, aero_src, sizeof(aero_t));
    
  /* Initialize... */
  if(init){
    for(im=0; im<aero_dest->nm; im++) {
      aero_dest->top_mod[im]=0;
      aero_dest->bottom_mod[im]=0;
      aero_dest->trans_mod[im]=0;
      aero_dest->nn[im]=0;
      aero_dest->rr[im]=0;
      aero_dest->ss[im]=1;
    }
    for(il=0; il<aero_dest->nl; il++) {
      aero_dest->nmod[il]=0;
      aero_dest->top[il]=0;
      aero_dest->bottom[il]=0;
      aero_dest->trans[il]=0;
      for(id=0; id<ctl->nd; id++){
	aero_dest->beta_e[il][id]=0;
	aero_dest->beta_s[il][id]=0;
	aero_dest->beta_a[il][id]=0;
	for(ia=0; ia<NTHETAMAX; ia++)
	  aero_dest->p[il][id][ia]=0;
      }
    }
  }
}

/*****************************************************************************/

void jur_sca_gauher(double *x,
	    double *w){

  /* Calculate abcissas (x) and weights (w) for Gauss-Hermite quadrature. */
  /* This routine basically follows the Numerical Recipes. */
  /* W. H. Press and S. A. Teukolsky and W. T. Vetterling and B. P. Flannery: */
  /* Numerical Recipes in C, 3rd edition, Cambridge University Press, 2007 */

  const double EPS=1.0e-14, PIM4=0.7511255444649425;

  const int MAXIT=10000;

  int i,its,j,m;

  double p1,p2,p3,pp,z=0,z1;

  int n=NRADMAX;

  m = (n+1)/2;
  for (i=0; i<m; i++) {
    if (i == 0) {
      z = sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if (i == 1) {
      z -= 1.14*pow((double)n,0.426)/z;
    } else if (i == 2) {
      z = 1.86*z-0.86*x[0];
    } else if (i == 3) {
      z = 1.91*z-0.91*x[1];
    } else {
      z = 2.0*z-x[i-2];
    }
    for (its=0; its<MAXIT; its++) {
      p1 = PIM4;
      p2 = 0.0;
      for (j=0; j<n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = z*sqrt(2.0/((double)j+1))*p2-sqrt((double)j/(j+1.))*p3;
      }
      pp = sqrt(2.*(double)n)*p2;
      z1 = z;
      z = z1-p1/pp;
      if (fabs(z-z1) <= EPS) break;
    }
    if (its >= MAXIT) {
      ERRMSG("Too many iterations in jur_sca_gauher.");
    }
    x[i] = z;
    x[n-1-i] = -z;
    w[i] = 2.0/(pp*pp);
    w[n-1-i] = w[i];
  }
}

/*****************************************************************************/
void jur_sca_get_opt_prop(ctl_t *ctl,
		  aero_t *aero){

  int nl=1, nm=1, count=0;
  int ii,jj, id, itheta;  
  double mbeta_e[NDMAX], mbeta_s[NDMAX], mp[NDMAX][NTHETAMAX];
 
  /* check input data */
  /* ToDo: improve check and sort data */
  if(aero->top_mod[0]<=aero->bottom_mod[0])
    ERRMSG("Aerosol top altitude is smaller than bottom altitude. Please check aero.tab.");

  aero->top[0] = aero->top_mod[0];
  aero->bottom[0] = aero->bottom_mod[0];
  aero->trans[0] = aero->trans_mod[0];
  
  for (ii=1; ii<aero->nm; ii++){
    if(aero->top_mod[ii]>aero->top_mod[ii-1] ||
       aero->top_mod[ii]<=aero->bottom_mod[ii] ||
       aero->bottom_mod[ii]>aero->bottom_mod[ii-1])
      ERRMSG("Aerosol/Cloud altitudes and/or transition layers are wrong. Please check aero.tab.");
    /* Identify number of aerosol/cloud layers and set top, */
    /* bottom and transition layer */
    /* Check if input is already sorted list from top to bottom; */
    if(aero->top_mod[ii]!=aero->top_mod[ii-1]){
      aero->top[nl] = aero->top_mod[ii];
      aero->bottom[nl] = aero->bottom_mod[ii];
      aero->trans[nl] = aero->trans_mod[ii];
      aero->nmod[nl-1] = nm;
      nl++;
      nm=0;
    }
    nm++;
  }
  aero->nmod[nl-1] = nm;
  aero->nl = nl;
  
  /* Get optical properties for each layer. */
  count=0;
  for (ii=0; ii<aero->nl; ii++){

    /* Initialise each layer */
    for(id=0; id<ctl->nd; id++){
      aero->beta_e[ii][id] = 0.;
      aero->beta_a[ii][id] = 0.;
      aero->beta_s[ii][id] = 0.;
      for(itheta=0; itheta<NTHETAMAX; itheta++){
	aero->p[ii][id][itheta] = 0.;
	mp[id][itheta] = 0.;
      }
      mbeta_e[id] = 0.;
      mbeta_s[id] = 0.;
    }

    /* Get optical properties for each mode. */
    for (jj=0; jj<aero->nmod[ii]; jj++){

      if(strcasecmp(aero->type[count], "MIE")==0){
    	/* Get optical properties for log-normal mode using Mie theory. */ 
	/* Gauss-Hermite integration */
	jur_sca_opt_prop_mie_log(ctl, aero, count, mbeta_e, mbeta_s, mp);
      } 
      else if(strcasecmp(aero->type[count], "Ext")==0){ 
	/* Get optical properties from external data base. Selects properties from closest wavenumber in data base file. */
	jur_sca_opt_prop_external(ctl, aero, count, mbeta_e, mbeta_s, mp);
      } else if(strcasecmp(aero->type[count], "Const")==0){ 
    	printf("Using constant extinction [1/km]: %g\n", aero->nn[count]);
    	ERRMSG("Implement me!");
      }
      else {
 	ERRMSG("Please give valid scattering model (MIE, Ext, Const)!");
      }
      
      /* Sum up optical properties for each layer */
      for(id=0; id<ctl->nd; id++){
    	aero->beta_e[ii][id] += mbeta_e[id];
    	aero->beta_s[ii][id] += mbeta_s[id];
     	aero->beta_a[ii][id] += (mbeta_e[id] - mbeta_s[id]);
    	for(itheta=0; itheta<NTHETAMAX; itheta++){
    	  aero->p[ii][id][itheta] += mp[id][itheta];
	  mp[id][itheta] = 0.;
	}
	mbeta_e[id] = 0.;
	mbeta_s[id] = 0.;
      }
      count++;
    }
  }
}

/*****************************************************************************/

void jur_sca_opt_prop_mie_log(ctl_t *ctl,
		     aero_t *aero,
		     int count,
		     double *beta_ext,
		     double *beta_sca,
		     double phase[NDMAX][NTHETAMAX]){

  FILE *in;
  
  char line[LENMAX];

  static int init=0;

  static double nu[REFMAX], nr[REFMAX], ni[REFMAX], n_imag[NDMAX], n_real[NDMAX], 
    rad_min=0.001, rad_max=1000., weights[NRADMAX], zs[NRADMAX];

  int npts=0, id, idx, nn, jj;

  double K1, rad, lambda, x, qext, qsca, qphase[NTHETAMAX];

  /* Read and interpolate refractive indices... */
  /* Check if previous mode has the same refractive index */
  if(count==0 || strcmp(aero->filepath[count], aero->filepath[count-1])!=0) { 
    
    /* Read data... */
    printf("Read refractive indices: %s\n", aero->filepath[count]);
    if(!(in=fopen(aero->filepath[count], "r")))
      ERRMSG("Cannot open file!");
    while(fgets(line, LENMAX, in))
      if(sscanf(line, "%lg %lg %lg", &nu[npts], &nr[npts], &ni[npts])==3)
	if((++npts)>REFMAX)
	  ERRMSG("Too many data points!");
    fclose(in);
    
    /* Interpolate... */
    for(id=0; id<ctl->nd; id++) {
      idx=jur_locate(nu, npts, ctl->nu[id]);
      n_real[id]=LIN(nu[idx], nr[idx], nu[idx+1], nr[idx+1], ctl->nu[id]);
      n_imag[id]=LIN(nu[idx], ni[idx], nu[idx+1], ni[idx+1], ctl->nu[id]);
    }
  } 

  /* Check log-normal parameters... */
  if(aero->nn[count]<=0 || aero->rr[count]<=0 || aero->ss[count]<=1)
    ERRMSG("The log-normal parameters are nonsense. ((?_?)) ");
    
  /* Integrate Mie parameters over log-normal mode */
  if(!init) {
    init=1;  
    /* get abcissas and weights for Gauss-Hermite quadrature */
    jur_sca_gauher(zs, weights); 
  }

  /* set coefficient */
  K1 = aero->nn[count] * 1e-3 * sqrt(M_PI);

  /* sum up Gaussian nodes */
  for (nn=0; nn<NRADMAX; ++nn) {
    rad = exp(sqrt(2) * log(aero->ss[count]) * zs[nn] + log(aero->rr[count]));
    if (rad >= rad_min && rad <= rad_max && K1 > 0.) {
 
     for(id=0; id<ctl->nd; id++){

	/* size parameter */
	lambda = 1./(ctl->nu[id])* pow(10,4.);
	x = 2*M_PI*rad/lambda;

	/* evaluate Mie Code at the nodes */
	jur_sca_bhmie(x, n_real[id], n_imag[id], qphase, &qext, &qsca);
	/* jur_sca_bhmie(rad, wavn, nang); */

	beta_ext[id] += K1 * pow(rad,2) *  qext * weights[nn];
	beta_sca[id] += K1 * pow(rad,2) *  qsca * weights[nn];
	    
	for (jj=0; jj<NTHETAMAX; ++jj)
	  phase[id][jj] += K1 * qsca * pow(rad,2) * qphase[jj] * weights[nn];

      } 
    } 
  }
  
  /* Weight phase function with beta_s */
  for(id=0; id<ctl->nd; id++){
    for (jj=0; jj<NTHETAMAX; ++jj)
      phase[id][jj] /= beta_sca[id];
  }
}

/*****************************************************************************/

void jur_sca_opt_prop_external(ctl_t *ctl,
		       aero_t *aero,
		       int count,
		       double *beta_ext,
		       double *beta_sca,
		       double phase[NDMAX][NTHETAMAX]){

  FILE *in;
  
  char line[LENMAX], *tok; 

  static double nu[REFMAX], n_ext[REFMAX], n_sca[REFMAX], n_phase[REFMAX][NTHETAMAX];

  int npts=0, ia, id, im;

  /* evtl. sp??ter Interpolation zw. Wellenl??ngen als Option einbauen */

  /* Read optical properties and find closest match to each wavenumber */
  if(aero->nn[count]==0) {
  
    /* Check if previous mode has the same optical properties */
    if(count==0 || strcmp(aero->filepath[count], aero->filepath[count-1])!=0) {
      
      /* Check for file... */
      printf("Read non-spherical optical properties: %s\n", aero->filepath[count]);
      if(!(in=fopen(aero->filepath[count], "r")))
	ERRMSG("Cannot open file!");
      
      /* Read data... */
      while(fgets(line, LENMAX, in)) {
    	
	TOK(line, tok, "%lg", nu[npts]);  
	TOK(NULL, tok, "%lg", n_ext[npts]);
	TOK(NULL, tok, "%lg", n_sca[npts]);
	for(ia=0; ia<NTHETAMAX; ia++)
	  TOK(NULL, tok, "%lg", n_phase[npts][ia]);
	
	if((++npts)>REFMAX)
	  ERRMSG("Too many data points!");
      }
      
      /* Close file... */
      fclose(in);
      
      /* Check number of points... */
      if(npts<1)
	ERRMSG("Could not read any data!");
    }

    /* Find closest match in wavenumber for each spectral point */
    for(id=0; id<ctl->nd; id++){
      im=jur_locate(nu, npts, ctl->nu[id]);
      if(im != npts && fabs(nu[im] - ctl->nu[id]) > fabs(nu[im+1] - ctl->nu[id]))
      /* if(im != npts && (nu[im] - ctl->nu[id]) > (nu[im+1] - ctl->nu[id])) */
	im=im+1;
      
      beta_ext[id] = n_ext[im];
      beta_sca[id] = n_sca[im];
      for (ia=0; ia<NTHETAMAX; ++ia)
      	phase[id][ia] = n_phase[im][ia];
    }
  }

  /* Read optical properties and interpolate to wavenumber */
  if(aero->nn[count]==1) {
    ERRMSG("Implement me!");
  }
}

/*****************************************************************************/

void jur_sca_read_aero(const char *dirname,
	       const char *filename,
	       ctl_t *ctl,
	       aero_t *aero){

  FILE *in;
  
  char file[LENMAX], line[LENMAX], *tok;
  
  /* Init... */
  aero->nm=0;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Read aerosol data: %s\n", file);
  
  /* Open file... */
  if(!(in=fopen(file, "r")))
    ERRMSG("Cannot open file!");
  
  /* Read data... */
  while(fgets(line, LENMAX, in)) {
    
    /* Read data... */
    TOK(line, tok, "%lg", aero->top_mod[aero->nm]);
    TOK(NULL, tok, "%lg", aero->bottom_mod[aero->nm]);
    TOK(NULL, tok, "%lg", aero->trans_mod[aero->nm]);
    TOK(NULL, tok, "%s",  aero->type[aero->nm][0]);
    TOK(NULL, tok, "%s",  aero->filepath[aero->nm][0]); 
    TOK(NULL, tok, "%lg", aero->nn[aero->nm]); 
    TOK(NULL, tok, "%lg", aero->rr[aero->nm]);
    TOK(NULL, tok, "%lg", aero->ss[aero->nm]);
    
    /* Increment counter... */
    if((++aero->nm)>SCAMODMAX)
      ERRMSG("Too many aerosol models!");
  }
  
  /* Close file... */
  fclose(in);
  
  /* Check number of points... */
  if(aero->nm<1)
    ERRMSG("Could not read any data!");

  /* Check consistency */
  if(ctl->sca_n != aero->nm)
    ERRMSG("Number of scattering models in control file and aerosol file does not match."); 
  /* printf("\nWARNING (%s, %s, l%d): %s\n\n",				\ */
  /* 	 __FILE__, __FUNCTION__, __LINE__, "Number of scattering models in control file and aerosol file does not match."); */

}

/*****************************************************************************/

void jur_sca_srcfunc_sca(ctl_t *ctl,
		 atm_t *atm,
		 aero_t *aero,
		 double sec,
		 double *x,
		 double *dx,
		 int il,
		 double *src_sca,
		 int scattering,
     queue_t *q) {
  
  /* Compute scattering of thermal radiation... */
  if(ctl->ip==1)
    jur_sca_srcfunc_sca_1d(ctl, atm, aero, sec, x, dx, il, src_sca, scattering, q); 
  else
    jur_sca_srcfunc_sca_3d(ctl, atm, aero, sec, x, dx, il, src_sca, scattering, q);
  
  /* Compute scattering of solar radiation... */
  if(TSUN>0)
    jur_sca_srcfunc_sca_sun(ctl, atm, aero, sec, x, dx, il, src_sca, q);
}

/*****************************************************************************/

void jur_sca_srcfunc_sca_1d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
        double sec,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering,
        queue_t *q) {
  
  obs_t *obs2;
  
  double alpha[NTHETAMAX], alpha2, dnorth[3], ek[3], lx[3], ly[3], lz[3], phi,
    phase2, rad, sx[3], sy[3], sz[3], theta[NTHETAMAX], theta2, w=0, wsum[NDMAX], xv[3];
  
  int i, id, idp, idx, iphi, ir, itheta, nalpha=28, nphi=180, ntheta2=180;

  int n1=1, n2=2;
  double midang=83, up=92, down=81, step=1+2;
  
  /* Allocate... */
  ALLOC(obs2, obs_t, 1);
  
  /* Set scattering phase function angles... */
  for(itheta=0; itheta<NTHETAMAX; itheta++)
    theta[itheta]=(double)itheta*M_PI/180.;
  
  /* Get local coordinate system... */
  dnorth[0]=-x[0];
  dnorth[1]=-x[1];
  dnorth[2]=2*RE-x[2];
  jur_sca_bascoord(x, dnorth, lx, ly, lz);
  
  /* Set angles - tested version with nalpha=28 and Fibonacci Numbers */
  alpha[0] = 0;
  for (ir=nalpha/2-4; ir<nalpha/2+4; ++ir){
    alpha[ir] = midang*M_PI/180.;
    midang++;
  }
  for (ir=0; ir<nalpha/2-5; ++ir){
    alpha[nalpha/2+ir+4] = up*M_PI/180.;
    alpha[nalpha/2-ir-5] = down*M_PI/180.;
    up+=step;
    down-=step;
    if (step < 13){
      step = (double)n1 + (double)n2;
      n1 = n2;
      n2 = (int)step;
    } 
  }
  alpha[nalpha-1] = M_PI;

  /* Get incident radiation... */
  /* nalpha=181; */
  for(ir=0; ir<nalpha; ir++) {
    
    /* Set angle... */
    /* traditional 0-180 deg in 1 deg steps with nalpha=181 */
    /* alpha[ir] = ir*M_PI/180.; */

    /* initial version with nalpha=21 */
    /* alpha[ir]=acos(2*(double)ir/(nalpha-1.0)-1.0); */
  
    /* Set view point... */
    for(i=0; i<3; i++)
      /* xv[i]=x[i]+10.*cos(alpha[ir])*lz[i]+10.*sin(alpha[ir])*ly[i];  */
      xv[i]=x[i]+10.*cos(alpha[ir])*(-1)*lz[i]+10.*sin(alpha[ir])*ly[i];
    
    /* Set observation geometry... */
    obs2->nr=nalpha;
    jur_cart2geo(x, &obs2->obsz[ir], &obs2->obslon[ir], &obs2->obslat[ir]);
    jur_cart2geo(xv, &obs2->vpz[ir], &obs2->vplon[ir], &obs2->vplat[ir]);
    obs2->time[ir] = sec; 
    /* Get pencil beam radiance... */
    jur_sca_formod_pencil(ctl, atm, obs2, aero, scattering-1, ir, q); 
  }
  if (Queue_Prepare == ctl->queue_state) return; /* prepare work queue items only */
  
  /* Get orthonormal basis (with respect to LOS)... */
  jur_sca_bascoord(dx, x, sx, sy, sz);  
  
  /* Initialize... */
  for(id=0; id<ctl->nd; id++){
    src_sca[id]=0;
    wsum[id]=0;
  }
    
  /* Loop over phase function angles... */
  for(itheta=0; itheta<ntheta2; itheta++) {
    
    /* Set phase function angle in 1?? steps (avoid 0 and 180??)... */
    theta2=(0.5+itheta)*M_PI/180.;
    
    /* Loop over azimuth angles... */
    for(iphi=0; iphi<nphi; iphi++) {
      
      /* Set azimuth angle in 2?? steps... */
      phi=2.*(0.5+iphi)*M_PI/180.;
      
      /* Get unit vector on sphere... */
      for(i=0; i<3; i++)
	ek[i]
	  =sin(theta2)*sin(phi)*sx[i]
	  +sin(theta2)*cos(phi)*sy[i]
	  +cos(theta2)*sz[i];

      /* Get phase function index */
      idp=jur_locate(theta, NTHETAMAX, theta2);

      /* Get zenith angle... */
      alpha2=ANGLE(-1.*lz, ek);

      /* Get source ray angle index */
      idx=jur_locate(alpha, nalpha, alpha2);
	
      /* Loop over channels... */
      for(id=0; id<ctl->nd; id++) {
	
	/* Interpolate phase function... */
	phase2=LIN(theta[idp], aero->p[il][id][idp],
		   theta[idp+1], aero->p[il][id][idp+1], theta2);

	/* Get weighting factor (area of surface element * phase function)... */
	w=sin(theta2)*phase2;
	
	/* Interpolate radiance to particular angle... */
	rad=LIN(alpha[idx], obs2->rad[idx][id], //CHANGED
	 	alpha[idx+1], obs2->rad[idx+1][id], alpha2); //CHANGED
	
	/* Integrate... */
	src_sca[id]+=w*rad;
	wsum[id]+=w;
      }
    }
  }
  
  /* Normalize... */
  for(id=0; id<ctl->nd; id++)
    src_sca[id]/=wsum[id]; 
  
  /* Free... */
  free(obs2);
}

/*****************************************************************************/

void jur_sca_srcfunc_sca_3d(ctl_t *ctl,
		    atm_t *atm,
		    aero_t *aero,
        double sec,
		    double *x,
		    double *dx,
		    int il,
		    double *src_sca,
		    int scattering,
        queue_t *q) {
  
  obs_t *obs2;
  
  double phi, phase2, sx[3], sy[3], sz[3], theta[NTHETAMAX], theta2,
    w, wsum=0, xv[3];
  
  int i, id, idx, iphi, itheta, nphi=180, ntheta2=180;
  
  /* Allocate... */
  ALLOC(obs2, obs_t, 1);
  
  /* Set scattering phase function angles... */
  for(itheta=0; itheta<NTHETAMAX; itheta++)
    theta[itheta]=M_PI*(double)itheta/(NTHETAMAX-1);
  
  /* Get orthonormal basis (with respect to LOS)... */
  jur_sca_bascoord(dx, x, sx, sy, sz);  
  
  /* Initialize... */
  for(id=0; id<ctl->nd; id++)
    src_sca[id]=0;
  
  /* Loop over phase function angles... */
  for(itheta=0; itheta<ntheta2; itheta++) {
    
    /* Set phase function angle... */
    theta2=(0.5+itheta)/ntheta2*M_PI;
    
    /* Loop over azimuth angles... */
    for(iphi=0; iphi<nphi; iphi++) {
      
      /* Set azimuth angle... */
      phi=(0.5+iphi)/nphi*2*M_PI;
      
      /* Set view point... */
      for(i=0; i<3; i++)
	xv[i]=x[i]
	  +10*sin(theta2)*sin(phi)*sx[i]
	  +10*sin(theta2)*cos(phi)*sy[i]
	  +10*cos(theta2)*sz[i];
      
      /* Set observation geometry... */
      obs2->nr=1;
      jur_cart2geo(x, &obs2->obsz[0], &obs2->obslon[0], &obs2->obslat[0]);
      jur_cart2geo(xv, &obs2->vpz[0], &obs2->vplon[0], &obs2->vplat[0]);
      obs2->time[0] = sec; 
      
      /* Get incident radiation... */
      jur_sca_formod_pencil(ctl, atm, obs2, aero, scattering-1, 0, q);
      
      /* Get phase function index */
      idx=jur_locate(theta, NTHETAMAX, theta2);

      /* Loop over channels... */
      for(id=0; id<ctl->nd; id++) {

        /* Interpolate phase function... */
        phase2=LIN(theta[idx], aero->p[il][id][idx],
            theta[idx+1], aero->p[il][id][idx+1], theta2);

        /* Get weighting factor (surface element area * phase function)... */
        w=M_PI/ntheta2*2*M_PI*sin(theta2)/nphi*phase2;

        /* Integrate... */
        src_sca[id]+=w*obs2->rad[0][id]; //CHANGED
        wsum+=w;
      }
    }
  }
  if (Queue_Prepare == ctl->queue_state) return; /* prepare work queue items only */
  
  /* Normalize... */
  for(id=0; id<ctl->nd; id++)
    src_sca[id]/=wsum;
  
  /* Free... */
  free(obs2);
}

/*****************************************************************************/

void jur_sca_srcfunc_sca_sun(ctl_t *ctl,
		     atm_t *atm,
		     aero_t *aero,
		     double sec,
		     double *x,
		     double *dx,
		     int il,
		     double *src_sca,
         queue_t *q) {

  pos_t *los;
  int np;
  double tsurf;
  
  obs_t *obs;
  
  double azi, dnorth[3], dout[3], ek[3], lx[3], ly[3], lz[3], phase2, sza,
    sza_beam, sza_cor, theta[NTHETAMAX], theta2, x0[3], x1[3];
  
  int i, i2, id, idx, itheta;
 
  /* Allocate... */
  los = (pos_t*) malloc((NLOSMAX) * sizeof(pos_t));
  ALLOC(obs, obs_t, 1);

  /* Set scattering phase function angles... */
  for(itheta=0; itheta<NTHETAMAX; itheta++)
    theta[itheta]=M_PI*(double)itheta/(NTHETAMAX-1);
  
  /* Get local coordinate system... */
  dnorth[0]=-x[0];
  dnorth[1]=-x[1];
  dnorth[2]=2*RE-x[2];
  jur_sca_bascoord(x, dnorth, lx, ly, lz);
  
  /* Get geometric coordinates of the Sun... */
  obs->nr=1;
  jur_cart2geo(x, &obs->obsz[0], &obs->obslon[0], &obs->obslat[0]);
  jur_sca_suncoord(sec, obs->obslon[0], obs->obslat[0], &azi, &sza);
  obs->time[0] = sec; // this was missing! 
  
  /* Find true elevation angle of Sun... */
  sza_cor=sza;
  for(i2=0; i2<10; i2++) {

    /* Set observation geometry... */
    for(i=0; i<3; i++) {
      ek[i]
        =sin(sza_cor*M_PI/180)*sin(azi*M_PI/180)*lx[i]
        +sin(sza_cor*M_PI/180)*cos(azi*M_PI/180)*ly[i]
        +cos(sza_cor*M_PI/180)*lz[i];
      x1[i]=x[i]+10*ek[i];
    }
    jur_cart2geo(x1, &obs->vpz[0], &obs->vplon[0], &obs->vplat[0]);

    /* Get zenith angle at end of beam... */
    np = jur_traceray(ctl, atm, obs, 0, los, &tsurf, aero, 1); // with scattering

    if(np<2)
      break;
    jur_geo2cart(los[np-2].z, los[np-2].lon, los[np-2].lat, x0);
    jur_geo2cart(los[np-1].z, los[np-1].lon, los[np-1].lat, x1);
    for(i=0; i<3; i++)
      dout[i]=x1[i]-x0[i];
    sza_beam=ANGLE(x, dout)*180/M_PI;

    /* Test for convergence... */
    if(fabs(sza_beam-sza)<0.01)
      break;

    /* Adapt geometric solar zenith angle (0.61803 golden ratio)... */
    sza_cor-=0.61803*(sza_beam-sza);
    sza_cor=GSL_MIN(GSL_MAX(sza_cor, 0), 180);
  }
 
  /* Check that LOS doesn't hit the ground... */
  if(tsurf<0) {
    /* Compute path transmittance... */
    jur_sca_formod_pencil(ctl, atm, obs, aero, 0, 0, q);
    
    // this was also missing!
    if (Queue_Prepare == ctl->queue_state) { /* prepare work queue items only */
      free(los);
      return;
    }
    
    /* Get phase function position... */
    theta2=ANGLE(ek, dx);
    idx=jur_locate(theta, NTHETAMAX, theta2);

    /* Loop over channels... */
    for(id=0; id<ctl->nd; id++) {
      
      /* Get phase function... */
      phase2=LIN(theta[idx], aero->p[il][id][idx],
		 theta[idx+1], aero->p[il][id][idx+1], theta2);
      
      /* Add solar radiance (6.764e-5 solid angle of the sun)... */
      src_sca[id] += 6.764e-5 * phase2 * jur_planck(TSUN, ctl->nu[id]) * obs->tau[0][id]; //CHANGED
    }
  }

  /* Free... */
  free(los);
  free(obs);
}

/*****************************************************************************/

void jur_sca_suncoord(double sec,
	      double lon,
	      double lat,
	      double *azi,
	      double *sza) {
  
  double D, dec, e, g, GMST, h, L, LST, q, ra;
  
  /* Number of days and fraction with respect to 2000-01-01T12:00Z... */
  D=sec/86400-0.5;
  
  /* Geocentric apparent ecliptic longitude [rad]... */
  g=(357.529+0.98560028*D)*M_PI/180;
  q=280.459+0.98564736*D;
  L=(q+1.915*sin(g)+0.020*sin(2*g))*M_PI/180;
  
  /* Mean obliquity of the ecliptic [rad]... */
  e=(23.439-0.00000036*D)*M_PI/180;
  
  /* Declination [rad]... */
  dec=asin(sin(e)*sin(L));
  
  /* Right ascension [rad]... */
  ra=atan2(cos(e)*sin(L), cos(L));
  
  /* Greenwich Mean Sidereal Time [h]... */
  GMST=18.697374558+24.06570982441908*D;
  
  /* Local Sidereal Time [h]... */
  LST=GMST+lon/15;
  
  /* Hour angle [rad]... */
  h=LST/12*M_PI-ra;
  
  /* Convert latitude... */
  lat*=M_PI/180;
  
  /* Azimuth [deg]... */
  *azi=atan2(-sin(h), cos(lat)*tan(dec)-sin(lat)*cos(h))*180/M_PI;
  
  /* Solar zenith angle (90 deg - elevation) [deg]... */
  *sza=acos(sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(h))*180/M_PI;
}

/*****************************************************************************/

void jur_sca_write_aero(const char *dirname,
		const char *filename,
		aero_t *aero) {
  
  FILE *out;
  
  char file[LENMAX];
  
  int im;
  
  /* Set filename... */
  if(dirname!=NULL)
    sprintf(file, "%s/%s", dirname, filename);
  else
    sprintf(file, "%s", filename);
  
  /* Write info... */
  printf("Write particle data: %s\n", file);
  
  /* Create file... */
  if(!(out=fopen(file, "w")))
    ERRMSG("Cannot create file!");
  
  /* Write header... */
  fprintf(out,
	  "# $1 = aerosol layer top altitude [km]\n"
	  "# $2 = aerosol layer bottom altitude [km]\n"
	  "# $3 = transition layer thickness [km]\n"
	  "# $4 = source for optical properties\n"
	  "# $5 = refractive index file\n"
	  "# $6 = particle concentration of log-normal mode [cm-3]\n"
	  "# $7 = median radius of log-normal mode [mum]\n"
	  "# $8 = width of log-normal mode\n");
  
  /* Write data... */
  for(im=0; im<aero->nm; im++) {
    if(im==0)
      fprintf(out, "\n");
    fprintf(out, "%g %g %g %s %s %g %g %g", aero->top_mod[im], 
	    aero->bottom_mod[im], aero->trans_mod[im], 
	    aero->type[im], aero->filepath[im], 
	    aero->nn[im], aero->rr[im], aero->ss[im]);
    fprintf(out, "\n");
  }
  
  /* Close file... */
  fclose(out);
}

/*****************************************************************************/
