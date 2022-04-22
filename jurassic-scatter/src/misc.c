#include "misc.h"

/*****************************************************************************/

// removing tbl_t
// void init_tbl(ctl_t *ctl,
//         tbl_t *tbl) {
//   
//   FILE *in;
//   
//   char filename[2*LEN], line[LEN];
//   
//   double eps, eps_old, press, press_old, temp, temp_old, u, u_old,
//     f[NSHAPE], fsum, nu[NSHAPE], tmin=100, tmax=400;
//   
//   int i, id, ig, ip, it, n;
//   
//   /* Loop over trace gases and channels... */
//   for(ig=0; ig<ctl->ng; ig++)
//     for(id=0; id<ctl->nd; id++) {
//       
//       /* Set filename... */
//       sprintf(filename, "%s_%.4f_%s.bin",
//         ctl->tblbase, ctl->nu[id], ctl->emitter[ig]);
//       
//       /* Try to open binary file... */
//       if((in=fopen(filename, "r"))) {
//   
//   /* Write info... */
//   LOGMSG(2, printf("Read emissivity table: %s\n", filename));
//   
//   /* Read data... */
//   FREAD(&tbl->np[ig][id], int, 1, in);
//   if(tbl->np[ig][id]>TBLNPMAX)
//     ERRMSG("Too many pressure levels!");
//   FREAD(&tbl->p[ig][id], double, tbl->np[ig][id], in);
//   FREAD(tbl->nt[ig][id], int, tbl->np[ig][id], in);
//   for(ip=0; ip<tbl->np[ig][id]; ip++) {
//     if(tbl->nt[ig][id][ip]>TBLNTMAX)
//       ERRMSG("Too many temperatures!");
//     FREAD(&tbl->t[ig][id][ip], double, tbl->nt[ig][id][ip], in);
//     FREAD(tbl->nu[ig][id][ip], int, tbl->nt[ig][id][ip], in);
//     for(it=0; it<tbl->nt[ig][id][ip]; it++) {
//       FREAD(&tbl->u[ig][id][ip][it], float,
//       GSL_MIN(tbl->nu[ig][id][ip][it], TBLNUMAX), in);
//       FREAD(&tbl->eps[ig][id][ip][it], float,
//       GSL_MIN(tbl->nu[ig][id][ip][it], TBLNUMAX), in);
//     }
//   }
//   
//   /* Close file... */
//   fclose(in);
//       }
//       
//       /* Try to read ASCII file... */
//       else {
//   
//   /* Initialize... */
//   tbl->np[ig][id]=-1;
//   eps_old=-999;
//   press_old=-999;
//   temp_old=-999;
//   u_old=-999;
//   
//   /* Try to open file... */
//   sprintf(filename, "%s_%.4f_%s.tab",
//     ctl->tblbase, ctl->nu[id], ctl->emitter[ig]);
//
//    if(!(in=fopen(filename, "r"))) {
//      LOGMSG(2, printf("Missing emissivity table: %s\n", filename));
//      continue;
//    }
//    LOGMSG(2, printf("Read emissivity table: %s\n", filename));
//    
//    /* Read data... */
//    while(fgets(line, LEN, in)) {
//      
//      /* Parse line... */
//      if(sscanf(line,"%lg %lg %lg %lg", &press, &temp, &u, &eps)!=4)
//        continue;
//      
//      /* Determine pressure index... */
//      if(press!=press_old) {
//        press_old=press;
//        if((++tbl->np[ig][id])>=TBLNPMAX)
//          ERRMSG("Too many pressure levels!");
//        tbl->nt[ig][id][tbl->np[ig][id]]=-1;
//      }
//      
//      /* Determine temperature index... */
//      if(temp!=temp_old) {
//        temp_old=temp;
//        if((++tbl->nt[ig][id][tbl->np[ig][id]])>=TBLNTMAX)
//          ERRMSG("Too many temperatures!");
//        tbl->nu[ig][id][tbl->np[ig][id]]
//          [tbl->nt[ig][id][tbl->np[ig][id]]]=-1;
//      }
//      
//      /* Determine column density index... */
//      if((eps>eps_old && u>u_old) || tbl->nu[ig][id][tbl->np[ig][id]]
//        [tbl->nt[ig][id][tbl->np[ig][id]]]<0) {
//        eps_old=eps;
//        u_old=u;
//        if((++tbl->nu[ig][id][tbl->np[ig][id]]
//      [tbl->nt[ig][id][tbl->np[ig][id]]])>=TBLNUMAX) {
//          tbl->nu[ig][id][tbl->np[ig][id]]
//      [tbl->nt[ig][id][tbl->np[ig][id]]]--;
//          continue;
//        }
//      }
//      
//      /* Store data... */
//      tbl->p[ig][id][tbl->np[ig][id]]=press;
//      tbl->t[ig][id][tbl->np[ig][id]][tbl->nt[ig][id][tbl->np[ig][id]]]
//        =temp;
//      tbl->u[ig][id][tbl->np[ig][id]][tbl->nt[ig][id][tbl->np[ig][id]]]
//        [tbl->nu[ig][id][tbl->np[ig][id]]
//        [tbl->nt[ig][id][tbl->np[ig][id]]]]=(float)u;
//      tbl->eps[ig][id][tbl->np[ig][id]][tbl->nt[ig][id][tbl->np[ig][id]]]
//        [tbl->nu[ig][id][tbl->np[ig][id]]
//        [tbl->nt[ig][id][tbl->np[ig][id]]]]=(float)eps;
//    }
//    
//    /* Increment counters... */
//    tbl->np[ig][id]++;
//    for(ip=0; ip<tbl->np[ig][id]; ip++) {
//      tbl->nt[ig][id][ip]++;
//      for(it=0; it<tbl->nt[ig][id][ip]; it++)
//        tbl->nu[ig][id][ip][it]++;
//   }
//   
//   /* Close file... */
//   fclose(in);
//       }
//     }
//   
//   /* Write info... */
//   LOGMSG(2, printf("Initialize source function table...\n"));
//   
//   /* Loop over channels... */
//   for(id=0; id<ctl->nd; id++) {
//     
//     /* Read filter function... */
//     sprintf(filename, "%s_%.4f.filt", ctl->tblbase, ctl->nu[id]);
//     read_shape(filename, nu, f, &n);
//     
//     /* Compute source function table... */
//     for(it=0; it<TBLNSMAX; it++) {
//       
//       /* Set temperature... */
//       tbl->st[it]=LIN(0.0, tmin, TBLNSMAX-1.0, tmax, (double)it);
//       
//       /* Integrate Planck function... */
//       fsum=0;
//       tbl->sr[id][it]=0;
//       for(i=0; i<n; i++) {
//   fsum+=f[i];
//   tbl->sr[id][it]+=f[i]*planck(tbl->st[it], nu[i]);
//       }
//       tbl->sr[id][it]/=fsum;
//     }
//   }
// }

/*****************************************************************************/
