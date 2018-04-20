/*
* =============================================================
* DTWave.c
*
* Mex C function to calculate the dynamic time warp distance between
* two arrays. The arrays may be of different length. Also computes a 
* weighted average vector sequence.
*
*
* =============================================================
*/

#include <mex.h>
/* #include <matrix.h> */
#include <math.h>
/* #include <mcheck.h> for debugging */
/* #include <stdio.h> for debugging */

/* PROTOTYPES */
double warpav(double *traceb, double *d, int numrows, int numcols, double *indelc, int fe, double *refend);
double warpav_ce(double *traceb, double *d, int numrows, int numcols, double *indelc, int fe, double *refend);
double *averseq(int *duree, double *traceb, double *mat1, double *mat2, double weight, int numrows, int numcols, int numv, double *refstart);
double *averseq_fe(int *duree, double *traceb, double *mat1, double *mat2, double weight, int numrows, int numcols, int numv1, int numv2, double *refstart, double indel_limit);
int mincost(double cost[], int a, int b);
double louiround(double x);
void euclidist(double *d, double *mat1, double *mat2, double *indelc, int numrows, int numcols, int numv1, int numv2, double *cow, int fe);
void vecdist(double *d, double *mat1, double *mat2, double *indelc, int numrows, int numcols, int numv1, int numv2, double *cow, int fe);
char printbase(double *mat, int position);
double *windel(double *mataverage, int numv1, int *duree, int win_size, double indel_limit);

double warpav_ce(double *traceb, double *d, int numrows, int numcols, double *indelc, int fe, double *refend)
{
    int i,j,k = 0 ; /* i rows ; j cols ; k traceback */
    double cost[5] ;
    double *dsf,dist ;

    dsf = (double *)mxMalloc((numrows+1) * (numcols+1) * sizeof(double));

    if (fe){
        printf("warpav_ce fe\n");
        /* first cell (1,1) */
        *dsf = 0;
        *traceb = 0;
        /* first column */
       for (i=1;i<=numcols;i++) {
           *(dsf+((numrows+1)*i)) = 0;
           *(traceb+((numrows+1)*i)) = 1;
       }
        /* first row */
       for (j=1;j<=numrows;j++) {
           *(dsf+j) = 0;
           *(traceb+j) = 3;
       }
    }else{
        /* first cell (1,1) */
        *dsf = *d * 2;
        *traceb = 0;
       /* first column */
       for (i=1;i<=numcols;i++) {
           *(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
           *(traceb+((numrows+1)*i)) = 1;
       }
        /* first row */
       for (j=1;j<=numrows;j++) {
           *(dsf+j) = *(dsf+(j-1)) + *indelc;
           *(traceb+j) = 3;
       }
    }

    /* cell (2,2) */
   cost[1] = *(dsf+1) + *indelc;
   cost[2] = *dsf + *d;
   cost[3] = *(dsf+numrows+1) + *indelc;
   k = mincost(cost,1,3);
   *(traceb+numrows+2) = k;
   *(dsf+numrows+2) = cost[k];

    /* second column */
   j=1;
   for (i=2;i<=numcols;i++) {
       cost[0] = *(dsf+((numrows+1)*(i-2))+(j-1)) + *(d+((numrows)*(i-2)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-1)))*0.5;
       cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
       cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
       cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
       k = mincost(cost,0,3);
       *(traceb+((numrows+1)*i)+1) = k;
       *(dsf+((numrows+1)*i)+1) = cost[k];
   }

    /* second row */
   i=1;
   for (j=2;j<=numrows;j++) {
       cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
       cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
       cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
       cost[4] = *(dsf+((numrows+1)*(i-1))+(j-2)) + *(d+((numrows)*(i-1)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-2)))*0.5;
       k = mincost(cost,1,4);
       *(traceb+((numrows+1))+j) = k;
       *(dsf+(numrows+1)+j) = cost[k];
   }

    /* rest of the matrix */
   for (i=2;i<=numcols;i++) {
       for (j=2;j<=numrows;j++) {
           cost[0] = *(dsf+((numrows+1)*(i-2))+(j-1)) + *(d+((numrows)*(i-2)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-1)))*0.5;
           cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
           cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
           cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
           cost[4] = *(dsf+((numrows+1)*(i-1))+(j-2)) + *(d+((numrows)*(i-1)+(j-1)))*0.5 + *(d+((numrows)*(i-1)+(j-2)))*0.5;
           k = mincost(cost,0,4);
           *(traceb+((numrows+1)*i)+j) = k;
           *(dsf+((numrows+1)*i)+j) = cost[k];
       }
   }
/* print matrices
for (i=0;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       printf(" %f", *(traceb+((numrows+1)*(i))+j));
   }
   printf("\n");
}
printf("\n");
for (i=0;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       printf(" %f", *(dsf+((numrows+1)*(i))+j));
   }
   printf("\n");
} */
   /* compute distance */
   if (fe) {
       /* end free space ***TO BE TESTED*** */
       double mini=*(dsf+((numrows+1)*(numcols+1))-1), minj=*(dsf+((numrows+1)*(numcols+1))-1) ;
       for (i=0;i<=numcols;i++) { if ((*(dsf+((numrows+1)*(i))+numrows))<mini) {mini=*(dsf+((numrows+1)*(i))+numrows);} }
       for (j=0;j<=numrows;j++) { if ((*(dsf+((numrows+1)*(numcols))+j))<minj) {minj=*(dsf+((numrows+1)*(numcols))+j);} }
       if (minj<mini) {dist = minj;} else {dist = mini;}
   } else {
       /* last cell */
       dist = (*(dsf+((numrows+1)*(numcols+1))-1)) / (double)(numrows+numcols);
   }

   mxFree(dsf);
   return dist ;
}

double warpav(double *traceb, double *d, int numrows, int numcols, double *indelc, int fe, double *refend)
{
    int i,j,k = 0 ; /* i rows ; j cols ; k traceback */
    double cost[4] ; /* need to keep 4 cells for compatibility with warpav_ce */
    double *dsf,dist ;

    dsf = (double *)mxMalloc((numrows+1) * (numcols+1) * sizeof(double));

    
    if (fe){/* column is pattern and row is reference */
        /*printf("warpav fe\n");*/
        /* first cell (1,1) */
        *dsf = 0;
        *traceb = 0;
        /* first column */
       for (i=1;i<=numcols;i++) {
           *(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
           *(traceb+((numrows+1)*i)) = 1;
       }
        /* first row */
       for (j=1;j<=numrows;j++) {
           *(dsf+j) = 0;
           *(traceb+j) = 3;
       }
    }else{
        /* first cell (1,1) */
        *dsf = *d * 2;
        *traceb = 0;
       /* first column */
       for (i=1;i<=numcols;i++) {
           *(dsf+((numrows+1)*i)) = *(dsf+((numrows+1)*(i-1))) + *indelc;
           *(traceb+((numrows+1)*i)) = 1;
       }
        /* first row */
       for (j=1;j<=numrows;j++) {
           *(dsf+j) = *(dsf+(j-1)) + *indelc;
           *(traceb+j) = 3;
       }
    }

    /* rest of the matrix */
   for (i=1;i<=numcols;i++) {
       for (j=1;j<=numrows;j++) {
           cost[1] = *(dsf+((numrows+1)*(i-1))+j) + *indelc;
           cost[2] = *(dsf+((numrows+1)*(i-1))+(j-1)) + *(d+((numrows)*(i-1)+(j-1)));
           cost[3] = *(dsf+((numrows+1)*(i))+(j-1)) + *indelc;
           k = mincost(cost,1,3);
           *(traceb+((numrows+1)*i)+j) = k;
           *(dsf+((numrows+1)*i)+j) = cost[k];
       }
   }

   /* compute distance */
    if (fe){/* end free space, look for minimum in the last row */
       /* double mini=*(dsf+((numrows+1)*(numcols+1))-1); */
       double minj=*(dsf+((numrows+1)*(numcols+1))-1);
       int lowest_cost_position=numrows;
       /* last column */ 
       /*for (i=0;i<=numcols;i++) { 
           if ((*(dsf+((numrows+1)*(i))+numrows))<mini) {
               mini=*(dsf+((numrows+1)*(i))+numrows);
           } 
       }*/   
       /* last row */
       for (j=0;j<=numrows;j++) { 
           if ((*(dsf+((numrows+1)*(numcols))+j))<minj) {
               minj=*(dsf+((numrows+1)*(numcols))+j);
               lowest_cost_position=j;
           } 
       }
       *refend = lowest_cost_position;
       /* if (minj<mini) {dist = minj;} else {dist = mini;} */
       dist = minj / (double)(numrows+numcols);
       /* dist = minj ; TO BE TESTED + NORMALIZE BY LENGTH OF ALIGNMENT */
       /* update the path at the end of the last row for deletions only, until pattern is found */
       /*printf("pos:%i,cost=%f\n",lowest_cost_position,minj);  */
       for (j=lowest_cost_position+1;j<=numrows;j++) {
           *(traceb+((numrows+1)*(numcols))+j) = 3;
       }
    }else{
       /* last cell */
       dist = (*(dsf+((numrows+1)*(numcols+1))-1)) / (double)(numrows+numcols);
    }

/* print matrices
for (i=-1;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       if (i==-1) {printf("\t|%d|", j);} 
       else {printf("\t%.0f", *(traceb+((numrows+1)*(i))+j));}
   }
   printf("\n");
}
printf("\n");
for (i=-1;i<=numcols;i++) {
   for (j=0;j<=numrows;j++) {
       if (i==-1) {printf("\t|%d|", j);} 
       else {printf("\t%.2f", *(dsf+((numrows+1)*(i))+j));}
   }
   printf("\n");
}
printf("\n"); */

   mxFree(dsf);
   return dist ;
}

/* compute the average sequence given the alignment */
double *averseq(int *duree, double *traceb, double *mat1, double *mat2, double weight, int numrows, int numcols, int numv, double *refstart)
{
    int i,j,n=0;
    int l=0;
    int maxt=0; /* length of the matrix */
    int ta,tb,temp;
    double *mat4,*t,*mata;

    mat4 = (double *)mxMalloc((numcols+numrows)*numv*sizeof(double)) ;
    t = (double *)mxMalloc((numcols+numrows)*sizeof(double)) ;

    i = numcols;
    j = numrows;

    /* trace back the alignment and implement the average vector sequence */
    while ( i>0 || j>0 ) {
       *(t+l) = weight*j + (1-weight)*i;           /*printf(" %f->%f\n",*(t+l),louiround(*(t+l)) );*/
       temp = (int)(louiround(*(t+l)));
       if ( temp>maxt ) { maxt=temp; } 
       switch ( (int) *(traceb+((numrows+1)*(i))+(j)) ) {
           case 0:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = weight*(*(mat1+n+(numv*(j-1)))) + 0.5*(1-weight)*(*(mat2+n+(numv*(i-1)))) + 0.5*(1-weight)*(*(mat2+n+(numv*(i-2)))) ; }
               i-=2;
               j--;
               break;
           case 1:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = (*(mat2+n+(numv*(i-1)))); }
               i--;
               break;
           case 2:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = weight*(*(mat1+n+(numv*(j-1)))) + (1-weight)*(*(mat2+n+(numv*(i-1)))); }
               i--;
               j--;
               break;
           case 3:
               for (n=0;n<numv;n++) { 
                    *(mat4+n+(numv*l)) = (*(mat1+n+(numv*(j-1)))); }
               j--;
               break;
           case 4:
               for (n=0;n<numv;n++) {
                    *(mat4+n+(numv*l)) = 0.5*weight*(*(mat1+n+(numv*(j-1)))) + 0.5*weight*(*(mat1+n+(numv*(j-2)))) + (1-weight)*(*(mat2+n+(numv*(i-1)))); }
               i--;
               j-=2;
               break;
           default:
               printf("error: traceback no match\n");
               break;
       }
       l++;
    }
    /*printf("\nnumcols=%i,numrows=%i,numv=%i,maxt=%i,l=%i,n_del=%i,n_ins=%i\n",numcols,numrows,numv,maxt,l,n_del,n_ins); */
    
    /* compute average with correct time direction and number of time points */
    mata = (double *)mxMalloc(maxt*numv*sizeof(double)) ;
    for ( ta=1 ; ta<=maxt ; ta++ ) {
       for ( tb=l-1 ; tb>=0 ; tb-- ) {
           if ( *(t+tb)==(float)ta ) {
               for ( n=0 ; n<numv ; n++ ) {
                    *(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
               break;
           } else if ( *(t+tb)>(float)ta ) {
               if ( tb==l-1 ) {
                    for ( n=0 ; n<numv ; n++ ) {
                            *(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
               } else {
                    for ( n=0 ; n<numv ; n++ ) {
                            *(mata+n+(numv*(ta-1))) = 0.5*(*(mat4+n+(numv*(tb+1)))) + 0.5*(*(mat4+n+(numv*(tb)))) ; }
               }
               break;
           } else if ( tb==0 ) {
               for ( n=0 ; n<numv ; n++ ) { 
                    *(mata+n+(numv*(ta-1))) = (*(mat4+n+(numv*(tb)))) ; }
           }
       }
    }
    *duree = maxt ;
    
    mxFree(mat4);
    mxFree(t);
    return mata ;
}

/* compute the average sequence given the alignment */
double *averseq_fe(int *duree, double *traceb, double *mat1, double *mat2, double weight, int numrows, int numcols, int numv1, int numv2, double *refstart, double indel_limit)
{
    int i,j,l,n=0;
    int maxt=0; /* length of the matrix */
    int ta,tb,temp;
    double *mat4,*mata,vsum;
    int *t; /* map path to reference position */
    int n_del=0, n_ins=0; /* count number of insertions and deletions into mat1 */
    int it; /* DEBUG */
    int refstart_saved=0;
    /* double indel_limit; determine how far the persistence vector needs to go before an INDEL is performed in reference (is now an argument) */
    
    /* indel_limit = 0.99 ;learning pull persistence toward 0 (deletion) or 2 (insertion) */
    /* indel_limit = weight/2 ; "/2" because it can go either up or down from 1 */
    /* indel_limit = weight*1.5 ;*/
    /* indel_limit = 0.03 ; */
    /* indel_limit = 0.5 ; DEBUG: set 9999999999999 to prevent any indels to be performed */
    /* printf( "%.2f\n", indel_limit ); */
    /* mat4 = (double *)mxMalloc((numcols+numrows)*(numv1)*sizeof(double)) ; FASTER */
    mat4 = (double *)mxCalloc((numcols+numrows)*(numv1),sizeof(double)) ; /* allocate and initialise to 0, needed to initialise coverage to 0 */
    t = (int *)mxMalloc((numcols+numrows)*sizeof(int)) ; /* store position index */

    i = numcols;/* mat2=read */
    j = numrows;/* mat1=reference */
    
    /* debug variables
    int db_n=0; */
    /*PRINT*/
    char read[numcols+numrows];
    char ref[numcols+numrows];
    char new[numcols+numrows];

    /* trace back the alignment and implement the average vector sequence */
    l=0;
    while ( i>0 || j>0 ) {
       *(t+l) = j;
       /*printf("i=%d,j=%d,l=%d,*(t+l)=%d,traceb=%d\n",i,j,l,*(t+l),(int)*(traceb+((numrows+1)*(i))+(j)));*/
       temp = (int)(louiround(*(t+l))); /* round needed? */
       if ( temp>maxt ) { maxt=temp; } 
       *(mat4+(numv1-1)+(numv1*l)) = *(mat1+(numv1-1)+(numv1*(j-1))) ; /* copy coverage from the reference */
       switch ( (int) *(traceb+((numrows+1)*(i))+(j)) ) {
           case 1:
               for (n=0;n<numv2;n++) { 
                   /**(mat4+n+(numv1*l)) = (*(mat1+n+(numv1*(j-1))));
                   *(mat4+n+(numv1*l)) = (*(mat2+n+(numv2*(i-1))));*/
                   *(mat4+n+(numv1*l)) = weight*( *(mat1+n+(numv1*(j-1))) ) + (1-weight)*( *(mat2+n+(numv2*(i-1))) ); /*INSERTION:weight the inserted base value with read value, it will be averaged with other bases mapping to the same position in the next loop (see counter)*/
               }
               if (j<numrows) { /* we are inside the alignment to reference (not in the free ends) */
                   *(mat4+n+(numv1*l)) = *(mat1+n+(numv1*(j-1))) + (1-weight);
                   /**(mat4+n+(numv1*l)) = weight*( *(mat1+n+(numv1*(j-1))) ) + (1-weight)*2;persistence*/ /*printf("ins %f,%f,%f\n",*(mat4+n+(numv1*l)),*(mat1+n+(numv1*(j-1))),1-weight);*/
                   if ( *(mat4+n+(numv1*l))>(1+indel_limit) ) {/**/
                       n_ins++;
                   }
               }else{
                   *(mat4+n+(numv1*l)) = *(mat1+n+(numv1*(j-1)));
               }
               /*printf("\ncase 1,n=%i,numv1=%i,numv2=%i,l=%i,pos=%i,val=%f\n",n,numv1,numv2,l,n+1+(numv1*l),*(mat4+n+1+(numv1*l)));*/
               read[l]=printbase(mat2,numv2*(i-1)); ref[l]='-'; new[l]=printbase(mat4,numv1*l);/* */
               i--;
               break;
           case 2:
               /*vsum=0; store the sum to normalise proba so that they sum to 1*/
               for (n=0;n<numv2;n++) { 
                    *(mat4+n+(numv1*l)) = weight*( *(mat1+n+(numv1*(j-1))) ) + (1-weight)*( *(mat2+n+(numv2*(i-1))) ); 
                    /* *(mat4+n+(numv1*l)) = roundf(*(mat4+n+(numv1*l)) * 10000) / 10000 ; */
                    /*vsum+=*(mat4+n+(numv1*l)); */
                    /*if ( *(mat1+n+(numv1*(j-1))) != *(mat2+n+(numv2*(i-1))) ) {db_n=1;}*/
               }
               /*if ( (round(vsum*10000)/10000)!=1 ) printf("%.16f\n",vsum);*/
               /*for (n=0;n<numv2;n++) {  add a normlisation step to ensure that the sum of the frequency vector sum to 1
                    *(mat4+n+(numv1*l)) = *(mat4+n+(numv1*l))/vsum; } */
               /* print alignment */
               /*if (db_n>0){
                   printf(" %d",j);
                   printf("\n%d  1\t  2\t  A\n",j);
                   for (n=0;n<numv2;n++) {printf("%.4f\t%.4f\t%.4f\n",*(mat1+n+(numv1*(j-1))),*(mat2+n+(numv2*(i-1))),*(mat4+n+(numv1*l)));}*/
               /*}*/
               /* persistence vector */
               *(mat4+n+(numv1*l)) = weight*( *(mat1+n+(numv1*(j-1))) ) + (1-weight); /* pull back persistence toward 1 */
               read[l]=printbase(mat2,numv2*(i-1)); ref[l]=printbase(mat1,numv1*(j-1)); new[l]=printbase(mat4,numv1*l); /* */
               n++; *(mat4+n+(numv1*l))+=1; /* coverage */
               i--;
               j--;
               break;
           case 3:
               for (n=0;n<numv2;n++) { 
                   *(mat4+n+(numv1*l)) = (*(mat1+n+(numv1*(j-1)))); /*DELETION:keep the reference value*/
               }
               if (i<numcols && i>0) { /* we are inside the alignment to reference (not in the free ends) */
                   *(mat4+n+(numv1*l)) = *(mat1+n+(numv1*(j-1))) - (1-weight) ;
                   /**(mat4+n+(numv1*l)) =  weight*( *(mat1+n+(numv1*(j-1))) ) ;persistence*/ /*printf("\ndel");*/
                   if ( *(mat4+n+(numv1*l))<(1-indel_limit) ) {/**/
                       n_del++;
                   }/**/
               }else {
                   *(mat4+n+(numv1*l)) = *(mat1+n+(numv1*(j-1)));
               }
               read[l]='-'; ref[l]=printbase(mat1,numv1*(j-1)); new[l]=printbase(mat4,numv1*l); /* */
               j--;
               break;
           default:
               printf("error: traceback no match\n");
               break;
       }
       if (refstart_saved==0 && i==0) { *refstart=(float)j; refstart_saved=1; } 
       /*printf("l=%d t=%d |",l,*(t+l));for (n=0;n<numv1;n++) { printf("%f ",l,*(mat4+n+(numv1*l))); } printf("|\n");*/
       l++;
    }
    /*if (n_ins>0 || n_del>0) {printf("n_ins=%i ; n_del=%i\n",n_ins,n_del);}*/
    /*printf("\nnumcols=%i(%i),numrows=%i(%i),numv1=%i,numv2=%i,l=%i,maxt+n_ins-n_del=%i(%i,%i,%i),refstart=%.0f\n",numcols,i,numrows,j,numv1,numv2,l,maxt+n_ins-n_del,maxt,n_ins,n_del,*refstart);*/
    /*PROBLEM HERE?:printf("t:\n"); for ( it=0 ; it<=(l-1) ; it++ ) { printf("%d\t",*(t+it)); } printf("\n");*/
    
    /*PRINT*/
    for (n=0;n<=(maxt+n_ins-n_del);n++){ printf("%d\t",*(t+n)); } printf("\n");/**/
    for (n=l-1;n>=0;n--){ printf("%d\t",*(t+n)); } printf("\n");/**/
    /* print a bar every 10 nucleotides */
    /* */ for (n=l-1;n>=0;n--){ if (n%10==0) {printf("|");} else {printf("\t");} } printf("\n"); 
    for (n=l-1;n>=0;n--){ printf("%c\t",ref[n]); } printf("\n");
    for (n=l-1;n>=0;n--){ printf("%c\t",read[n]); } printf("\n");
    for (n=l-1;n>=0;n--){ printf("%c\t",new[n]); } printf("\n");
    /*PRINT mat4*/
    printf("\n");
    for ( tb=0 ; tb<=l-1 ; tb++ ) { printf( "%d\t",tb); } printf("\n");
    for ( n=0 ; n<numv1 ; n++ ) {
        for ( tb=l-1 ; tb>=0 ; tb-- ) {
            printf( "%.2f\t", *(mat4+n+(numv1*(tb))) ); }
        printf("\n");
    }
    
    /* compute average with correct time direction and number of time points */
    int c,m,last_t=-1,counter=1 ;
    *duree = maxt+n_ins-n_del ;
    /*printf("allocate %d\n",*duree); */
    /* mata = (double *)mxMalloc((*duree)*numv1*sizeof(double)) ; FASTER */
    mata = (double *)mxCalloc((*duree)*numv1,sizeof(double)) ; /* allocate and initialise to 0, needed to initialise coverage to 0 */
    tb=l-1 ;
    /*printf("-> ta=%d, tb=%d, t+tb=%d\n",ta,tb,*(t+tb));*/
    for ( ta=1 ; ta<=(*duree) ; ta++ ) {
        /*printf("0 ta=%d, tb=%d, t+tb=%d mat4(4,tb)=%.2f \n",ta,tb,*(t+tb),(*(mat4+4+(numv1*(tb)))) ); */
        /*printf("mata(ta=%d), mat4(tb=%d), t+tb=%d, maxt=%d\n",ta,tb,*(t+tb),maxt); */
        /* PRINT ALIGNMENT
        for ( n=0 ; n<numv2 ; n++ ) {printf("") ;} */
        /* perform INDELS in reference */
        if ( *(mat4+numv2+(numv1*(tb)))>(1+indel_limit) ){ /* printf("INSERTION\n"); *//* */
            /*printf("INSERTION mata(ta=%d), mat4(tb=%d), t+tb=%d, last_t=%d\n",ta,tb,*(t+tb),last_t);*/
            for ( n=0 ; n<numv2 ; n++ ) { /* insert new base */
                *(mata+n+(numv1*(ta-1))) = *(mat4+n+(numv1*(tb))) ;
            }
            *(mata+n+(numv1*(ta-1))) = 1; /* new persistence value is 1 */
            n++; *(mata+n+(numv1*(ta-1))) = 1; /* new coverage value is 1 */
            last_t = *(t+tb) ;
            tb-- ; /* */
            for ( m=tb ; m>=0 ; m-- ) { *(t+m) = *(t+m)+1 ; } /* change the mapping position in the alignment for all successive positions */
        }else if ( *(mat4+numv2+(numv1*(tb)))<(1-indel_limit) ){ /* printf("DELETION\n"); *//* */
            /*printf("DELETION mata(ta=%d), mat4(tb=%d), t+tb=%d, last_t=%d\n",ta,tb,*(t+tb),last_t);*/
            ta-- ;
            tb-- ;
        }else{ /*printf("SUBSTITUTION\n");*/
            /* if (*t+tb<0) {printf("t must be positive");break;} */
            /* 1st option: keep the reference's vector value */
            /* 2nd option: compute mean of all position mapping to the same position: */
            if ( *(t+tb)==last_t ){ /* still aligning to the same position, compute the running average */
                /*printf("A ta(mata)=%d, tb(mat4,t)=%d, *(t+tb)=%d, last_t=%d, persist=%.4f\n",ta,tb,*(t+tb),last_t,*(mat4+numv1+(numv1*(tb)))); */
                counter++;
                ta--;
                for ( n=0 ; n<numv1-1 ; n++ ) { /* "counter" loop, do running average for persistence value too */
                    /*printf("mata=%.4f, mat4=%.4f, counter=%d\n",*(mata+n+(numv1*(ta-1))),*(mat4+n+(numv1*(tb))),counter); */
                    *(mata+n+(numv1*(ta-1))) = *(mata+n+(numv1*(ta-1))) + ( *(mat4+n+(numv1*(tb))) - *(mata+n+(numv1*(ta-1)))  ) / counter ;
                }
            }else{
                /*printf("B ta(mata)=%d, tb(mat4,t)=%d, *(t+tb)=%d, last_t=%d\n",ta,tb,*(t+tb),last_t);*/
                counter=1;
                for ( n=0 ; n<numv1-1 ; n++ ) { /* all bases + persistence value */
                    *(mata+n+(numv1*(ta-1))) = *(mat4+n+(numv1*(tb))) ;
                }
                /* coverage (not averaged with other positions aligning here) */
                *(mata+n+(numv1*(ta-1))) = *(mat4+n+(numv1*(tb))) ;
            }
            last_t = *(t+tb) ;
            tb-- ;
        }
    }
    
    /*PRINT mata
    for ( ta=0 ; ta<*duree ; ta++ ) { printf( "%d\t",ta); } printf("\n");
    for ( n=0 ; n<numv1 ; n++ ) {
        for ( ta=0 ; ta<*duree ; ta++ ) { printf( "%.2f\t", *(mata+n+(numv1*(ta))) ); }
        printf("\n");
    }*/
    
    mxFree(mat4);
    mxFree(t);
    return mata ;
}

/* find the minimum value in an array between indexes a and b
int mincost(double cost[], int a, int b)
{
    int j=a;
    while (a<b+1) {
        if (cost[a]<cost[j]) {j=a;}
        a++;
    }
    return j;
}*/
/* find the minimum value in an array between indexes a and b
and returns its index
if different values are equals, randomly chose one using time, quite slower */
int mincost(double cost[], int a, int b)
{
   int j=a;
   while (a<=b) {
       if (cost[a]<cost[j]) {j=a;}
       else if (cost[a]==cost[j]) {
            if (time(NULL)%2==0) {j=a;} }
       a++;
   }
   return j;
}

/* no round() function in math.h so here it is, SOMETIMES THIS FUNCTION IS BUILT IN */
double louiround(double x)
{
    double intpart;
    if(modf(x,&intpart)>=0.5){
      return ceil(x);}
    else{
      return floor(x);}
}

/* compute Euclidean distance between vectors of 2 matrices 
   the insertion deletion cost is set as (half) the average distance between vectors */
void euclidist(double *d, double *mat1, double *mat2, double *indelc, int numrows, int numcols, int numv1, int numv2, double *cow, int fe)
{
    int i,j,n ;

    for (i=0;i<numcols;i++) {
       for (j=0;j<numrows;j++) {
           *(d+((numrows)*i)+j)=0;
           if (fe){
               for ( n=0 ; n<numv2 ; n++ ) {
                   /* if (i==0) printf("%f\n",(*(mat1+n+((numv+1)*(j))))); */
                   *(d+((numrows)*i)+j) = *(d+((numrows)*i)+j) + ((*(mat2+n+(numv2*(i))))-(*(mat1+n+(numv1*(j))))) * ((*(mat2+n+(numv2*(i))))-(*(mat1+n+(numv1*(j))))) ;
               }
               *(d+((numrows)*i)+j) = sqrt( *(d+((numrows)*i)+j) ) ;
               /* printf("%f\n",*(d+((numrows)*i)+j)); */
           }else{
               for ( n=0 ; n<numv2 ; n++ ) {/* if (i==0) printf("%f\n",(*(mat1+n+((numv+1)*(j))))); */
                   *(d+((numrows)*i)+j) = *(d+((numrows)*i)+j) + *(cow+n) * ((*(mat2+n+(numv2*(i))))-(*(mat1+n+(numv1*(j))))) * ((*(mat2+n+(numv2*(i))))-(*(mat1+n+(numv1*(j))))) ; }
               *(d+((numrows)*i)+j) = sqrt( *(d+((numrows)*i)+j) ) ;
               *indelc += *(d+(numrows*i)+j) ;
               /* printf("%f\n",*(d+((numrows)*i)+j)); */
           }
       }
    }
    
    if (fe){ /* hard code the cost as twice the mismatch of fully defined base position (e.g. 2*d([1 0 0 0],[0 1 0 0]) )*/
        *indelc = 2*sqrt(2) ; /* */
        /* *indelc = sqrt(2) ; DEBUG force insertion to test average sequence computation */
        /* *indelc = numcols*sqrt(2) ; DEBUG to completely forbid indels */
    }else{
        *indelc = *indelc/(numrows*numcols) ;
    }
    
    return ;
}

/* compute distance between vectors of 2 matrices, mat2 is {0,1}* and looks only at position with {1} in mat2 in both matrices.
   UPDATE: look at position (base) in the read that has maximum probability (in case quality is taken into account) 
   the insertion deletion cost is set to 2 */
void vecdist(double *d, double *mat1, double *mat2, double *indelc, int numrows, int numcols, int numv1, int numv2, double *cow, int fe)
{
    int i,j,n,m ;

    for (i=0;i<numcols;i++) {
        for (j=0;j<numrows;j++) {
            *(d+((numrows)*i)+j)=0;
            m=0;
            for ( n=0 ; n<numv2 ; n++ ) {
                /* if (i==0) printf("%f\n",(*(mat1+n+((numv+1)*(j))))); */
                /* if ( (*(mat2+n+(numv2*(i))))==1 ){ only consider the base of the read that have 1 */
                if ( (*(mat2+n+(numv2*(i))))>m ){ /* UPDATE: only consider the base of the read with maximum probability */
                    m=(*(mat2+n+(numv2*(i)))) ;
                    *(d+((numrows)*i)+j) = fabs( (*(mat2+n+(numv2*(i))))-(*(mat1+n+(numv1*(j)))) ) ;
                }
            }
            /* printf("%f\n",*(d+((numrows)*i)+j)); */
        }
    }
    
    /* hard code the cost as twice the mismatch of fully defined base position (e.g. 2*d([1 0 0 0],[0 1 0 0]) ) */
    /* *indelc = 2 ; */
    /* *indelc = 99 ; */
    return ;
}

char printbase(double *mat, int position)
{
    double maxvalue=0.0;
    int n,maxbase=0;
    char c;
    
    for (n=0;n<4;n++) {
        if (*(mat+n+position)>maxvalue) { 
            maxbase=n;
            maxvalue=*(mat+n+position);
        }
    }
    switch ( maxbase ) {
           case 0:
               c='A'; break;
           case 1:
               c='C'; break;
           case 2:
               c='T'; break;
           case 3:
               c='G'; break;
           default:
               c='N';
    }
    return c;
}


double *windel(double *mataverage, int numv1, int *duree, int win_size, double indel_limit)
{
    int n,t;
    double *mata;
    
    mata = (double *)mxMalloc((*duree)*numv1*sizeof(double)) ;
    for ( n=0 ; n<numv1 ; n++ ) {
        for ( t=0 ; t<*duree ; t++ ) { printf( "%.2f\t", *(mataverage+n+(numv1*(t))) ); }
        printf("\n");
    }
    
    return mata;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *indelc, *d, *mat1, *mat2, *traceb, *mataverage, *distance, *cow ;
    int numrows, numcols, numv, numv1, numv2, *duree, ce, fe, cw ;
    /* mwSize numrows, numcols, numv, *duree ; */
    double weight ;
    double *refend=0 ;
    double *refstart=0 ;
    double x ;
    double indel_limit, indel_weight ;
    
    /* memory allocation */
    indelc = (double *)mxMalloc(sizeof(double));
    duree = (int *)mxMalloc(sizeof(int));
    refend = (double *)mxMalloc(sizeof(double));
    refstart = (double *)mxMalloc(sizeof(double));

    /* read in arguments */
    mat1 = mxGetPr(prhs[0]); /* reference */
    mat2 = mxGetPr(prhs[1]); /* read */
    cow = mxGetPr(prhs[2]); /* weight per base */
    weight = mxGetScalar(prhs[3]); /* weight on reference versus read */
    ce = mxGetScalar(prhs[4]); /* compression/expansion */
    fe = mxGetScalar(prhs[5]); /* free-end alignment */
    cw = mxGetScalar(prhs[6]); /* flag for computing weight (1) or using argumented value (0) */
    *indelc = mxGetScalar(prhs[7]); /* cost for indels */
    indel_weight = mxGetScalar(prhs[8]); /* factor on weight to compute the indel limit */
    
    numrows = mxGetN(prhs[0]); /* Number of columns in mxArray */
    numcols = mxGetN(prhs[1]);
    numv1 = mxGetM(prhs[0]); /* Number of rows in mxArray */
    numv2 = mxGetM(prhs[1]); /* Number of rows in mxArray */
    /*printf("weight=%f, ce=%d, fe=%d, numrows=%d, numcols=%d, numv1=%d, numv2=%d\n\n",weight,ce,fe,numrows,numcols,numv1,numv2);*/
    
    /* if free ends alignment, second matrix must be pattern to find in first matrix */
    if (fe){
        /*if (numcols>numrows) {
            mexErrMsgIdAndTxt("DTWave:mexFunction","Matrix 2 must be shorter than matrix 1");
        }*/
        if (numv1<5 || numv1>6) { /* 6 when coverage is saved */
            mexErrMsgIdAndTxt("DTWave:mexFunction","Reference matrix dimension error (numv1)");
        }
        if (numv2<4 || numv2>4) {
            mexErrMsgIdAndTxt("DTWave:mexFunction","Matrix dimension error (numv2)");
        }
        if (ce) {
            mexErrMsgIdAndTxt("DTWave:mexFunction","Free ends alignment is not implemented with compression/expansion yet");
        }
    }

    /* more memory allocation */
    d = (double *)mxMalloc((numrows) * (numcols) * sizeof(double));
    traceb = (double *)mxMalloc((numrows+1) * (numcols+1) * sizeof(double));
    
    /* Create variable for return values and get pointers */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL); /* DISTANCE */
    distance = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL); /* AVERAGE VECTOR SEQUENCE */
    plhs[2] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL); /* POSITION IN MATRIX 1 END */
    plhs[3] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL); /* POSITION IN MATRIX 1 START */

    /* compute pairwise vector distance */
    /* *indelc = 2 ; set indel cost (now an argument) */
    if (fe){
        vecdist(d,mat1,mat2,indelc,numrows,numcols,numv1,numv2,cow,fe);
        /* euclidist(d,mat1,mat2,indelc,numrows,numcols,numv1,numv2,cow,fe); */
    }else{
        euclidist(d,mat1,mat2,indelc,numrows,numcols,numv1,numv2,cow,fe);
    }
    
    /* Call the warpav subroutine */
    *refend = numrows ; /* by default the alignment length is the length of mat1 */
    if (ce){
        *distance = warpav_ce(traceb,d,numrows,numcols,indelc,fe,refend) ;
    }else{
        *distance = warpav(traceb,d,numrows,numcols,indelc,fe,refend) ;
    }
    
    numv = numv2 ; /* TO BE FIXED? */
    /* Call the traceback subroutine to get the average sequence */
    *refstart = 0 ; /* by default the alignment starts at 0 */
    if ( weight!=0 || cw==1) {
        if (fe){
            /* normalise distance by read length, the distance becomes a base alignment rate: average distance per base of the read */
            *distance = *distance/numv2 ;
            if (cw==1){
            /* modify weight according to alignment distance */ /* printf("weight=%f ",weight) ; printf("dist=%.16f, weight=%f\n",*distance,weight) ; */
            /* weight = weight * fmin( (*distance/0.0015)*((1/weight)-1)+1 , 1/weight ); */ /* printf("after weight, dist=%.16f, new weight=%f\n",*distance,weight) ; */
            /* x = 0.3 ;
            weight = 0.5*(*distance*x*1000-5)/(1+abs(*distance*x*10000-5)) + 0.5 ; approximation */
            /* weight = *distance*660 ; if (weight>0.99) {weight=0.99;} linear approximation: y = 0.99*x/1.5e-3 */
              if (*distance<.001) {weight=0.05;}
              else if (*distance>0.002) {weight=0.95;}
              else {weight=*distance*900-0.85;} /* approximation of logistic function by "cut" function */
            }
            indel_limit = weight*indel_weight ;
            printf("distance=%f - weight=%.3f - indel_weight=%.3f - indel_limit=%.3f\n", *distance, weight, indel_weight, indel_limit);/*  */
            mataverage = averseq_fe(duree,traceb,mat1,mat2,weight,numrows,numcols,numv1,numv2,refstart,indel_limit);
            /* normalise distance by length of alignment */
            /* *distance = *distance/(*refend-*refstart) ; NEED TESTING... */
            /* deal with indels using sliding windows
            mataverage = windel(mataverage,numv1,duree,5,0.03); */
            /* set the pointer for return matrix */
            mxSetPr( plhs[1], mataverage );
            mxSetM( plhs[1], numv1 );
            mxSetN( plhs[1], *duree );
            mexMakeMemoryPersistent(mataverage); /* memory deallocation bug with octave */
        }else{
            mataverage = averseq(duree,traceb,mat1,mat2,weight,numrows,numcols,numv,refstart);
            /* set the pointer for return matrix */
            mxSetPr( plhs[1], mataverage );
            mxSetM( plhs[1], numv );
            mxSetN( plhs[1], *duree );
            mexMakeMemoryPersistent(mataverage); /* memory deallocation bug with octave */
        }
    }
    mxSetPr( plhs[2], refend );
    mxSetPr( plhs[3], refstart );
    /* printf("dist=%f, mata[1,1]=%d, refend=%d, refstart=%d\n",plhs[0],plhs[1],plhs[2],plhs[3]) ;*/
    
    mxFree(duree);
    mxFree(d);
    mxFree(traceb);
    mxFree(indelc);
}
