#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

/* ********** */

void 
fwtnpperm(input,f,nkeep,intercept,initboundhandl,neighbours,closest,LocalPred,n,coeff,lengthsremove,lengths,lca,pointsin,nc,traj,doW,W,varonly,v)

double *input,*f,*coeff,*lca,*lengthsremove,*lengths;
int 
*nkeep,*intercept,*initboundhandl,*neighbours,*closest,*LocalPred,*n,*pointsin,*nc, *traj;  
double *W,*v;
int *doW, *varonly;
{

int i,j,k=0,N=*n,nn,scheme,r,remove,nr=0,max,dim,dim1,dim2,dim2sq,dummy=*n-*nkeep, nnmax=2**neighbours,one=1, ex=0;
int *po;     
int *nbrs2;    
int *index2;   
int *nbrs;       
int *index;      

double *X=malloc(*n*sizeof(double));    
double *I=malloc((*n+1)*sizeof(double));
double *sX=malloc(*n*sizeof(double));

double *weights2;   
double *len2;
double *newline;                             
double *tmplca;                              
double *alpha, *weights;                

double *Wnew=0,*Wtmp=0;

void adaptneigh();
void adaptpred();
void cubicpred();
void delrow();
void getnbrs();
void getridd();
void getridi();
void intervals();
void linearpred();
void makelcaline();
void mycpyd();
void mycpyi();
void mysortd();
void mymind();
void pointsupdate();
void pts();
void quadpred();
void updatelca();

mycpyd(input,n,X);
intervals(X,initboundhandl,n,I);
for(i=0;i<*n;i++){
    *(lengths+i)=*(I+i+1)-*(I+i);
}
free(I);

mysortd(X,n,sX,pointsin,&one);
mycpyd(f,n,coeff);
free(sX);

/* (doW,varonly) should be (0,0),(1,0),or(0,1) */

ex=*doW+*varonly;

if(ex==1){
        Wnew=calloc(*n**n,sizeof(double));
        for(i=0;i<*n;i++){
                *(Wnew+(i**n)+i)=1;
        }
}

/*	if(*varonly==1){
		*(v+i)=1;
	}
*/	
	

if (*nkeep!=*n) {
    for (j=1;j<=dummy;j++) {
/*	if((j%100)==0){
	printf("j:%d\n",j); 
	}
*/
    	remove=*(traj+j-1);
/*	printf("remove: %d\n",remove);        */
        nbrs=calloc(nnmax,sizeof(int));    /* set up as maximal */
        index=calloc(nnmax,sizeof(int));   /* ... */
               
       
        if(*LocalPred==5){          
        
            nbrs2=calloc(nnmax,sizeof(int));
            index2=calloc(nnmax,sizeof(int));

/*            mycpyi(nbrs,&nnmax,nbrs2);                   
            mycpyi(index,&nnmax,index2); 	wierd?         

		mycpyi(nbrs,&nn,nbrs2);
		mycpyi(index,&nn,index2);*/

            weights2=calloc(nnmax,sizeof(double));
        
        }
        else{                       /* known nn given by getnbrs */
        getnbrs(X, &remove, pointsin, &N,neighbours,closest,nbrs,index,&nn);
        
            nbrs2=calloc(nn,sizeof(int));    
            index2=calloc(nn,sizeof(int));   
        
            mycpyi(nbrs,&nn,nbrs2);                   /*   fill to proper size */
            mycpyi(index,&nn,index2);                 /*   ... */
        
            weights2=calloc(nn,sizeof(double));

        }
            free(nbrs);
            free(index);

       
        switch(*LocalPred){
            case 1: 
            scheme=1;
            linearpred(pointsin,X,coeff,nbrs2,&remove,intercept,&nn,weights2,&one);
            break;
            
            case 2: 
            scheme=2;
            quadpred(pointsin,X,coeff,nbrs2,&remove,intercept,&nn,weights2,&one);
            break;
            
            case 3: 
            scheme=3;
            cubicpred(pointsin,X,coeff,nbrs2,&remove,intercept,&nn,weights2,&one);
            break;
            
            case 4: 
            scheme=1;   
            adaptpred(pointsin,X,coeff,nbrs2,&remove,intercept,&nn,weights2,&scheme,&one);
            break;
            
            case 5: 
            scheme=1;   
            adaptneigh(pointsin,X,coeff,nbrs2,&remove,intercept,&nn,weights2,&scheme,closest,index2,neighbours,&N,&one);
            break;
        }
        
        nbrs=calloc(nn,sizeof(int));                /* nn should be known in all cases now */
        index=calloc(nn,sizeof(int));
        weights=calloc(nn,sizeof(double));
        alpha=calloc(nn,sizeof(double));
        
        mycpyi(nbrs2,&nn,nbrs);
        mycpyi(index2,&nn,index);
        mycpyd(weights2,&nn,weights);
        
        free(nbrs2);
        free(index2);
        free(weights2);
        
        pointsupdate(X,coeff,&nn,index,&remove,pointsin,weights,lengths,&N,alpha,&r);
                
        *(lengthsremove+j-1)=*(lengths+r-1);
        
        newline=calloc((3*nn+5),sizeof(double));
        makelcaline(&remove,&nn,nbrs,alpha,weights,&scheme,intercept,closest,newline);
        max=(*nc>=(3*nn+5))? *nc: (3*nn+5);         
        tmplca=calloc(max*j,sizeof(double));     
        updatelca(lca,&nr,nc,newline,tmplca);                     

	if(ex==1){
        	if(*varonly==1){
                  	for(i=0;i<*n;i++){
                                for(k=0;k<nn;k++){
                                        *(Wnew+(r-1)**n+i)-=*(weights+k)**(Wnew+(*(index+k)-1)**n+i);
                                }
                        }
                        for(i=0;i<*n;i++){
                                for(k=0;k<nn;k++){
                                        *(Wnew+(*(index+k)-1)**n+i)+=*(alpha+k)**(Wnew+(r-1)**n+i);
                                }
                                *(v+remove-1)+=pow(*(Wnew+(r-1)**n+i),2);       
                        }       
                	
                        dim=*n-j+1;
                	dim1=dim-1;
                	dim2=dim**n;
                	Wtmp=calloc(dim2,sizeof(double));
                	mycpyd(Wnew,&dim2,Wtmp);
			free(Wnew);
                	Wnew=calloc(dim1**n,sizeof(double));
                	delrow(Wtmp,&dim,n,&r,Wnew);
            	        free(Wtmp);
                }
                else{	
                      for(i=0;i<*n;i++){
                                for(k=0;k<nn;k++){
                                        *(Wnew+(remove-1)**n+i)-=*(weights+k)**(Wnew+(*(nbrs+k)-1)**n+i);  
                                }
                        }
                        for(i=0;i<*n;i++){
                                for(k=0;k<nn;k++){
                                        *(Wnew+(*(nbrs+k)-1)**n+i)+=*(alpha+k)**(Wnew+(remove-1)**n+i);         
                                }
                        }
                }	
        }	
                        
        free(nbrs);
        free(alpha);
        free(weights);
        free(newline);
        free(index);              
                                                      
        len2=calloc(N,sizeof(double));
        po=calloc(N,sizeof(int));
    
        mycpyd(lengths,&N,len2);      
        mycpyi(pointsin,&N,po);
     
        getridd(len2,&N,&r,lengths);
        getridi(po,&N,&r,pointsin);
        free(len2);
        free(po);
    
        dim=nr**nc;
        mycpyd(tmplca,&dim,lca);
        free(tmplca);
        N-=1;
   }  	/* j */
}	/* if */
   
if(ex==1){
        if(*varonly==1){
                for(i=0;i<N;i++){
                        for(k=0;k<*n;k++){
                                *(v+*(pointsin+i)-1)+=pow(*(Wnew+(i**n)+k),2);
                        }       
                }
                dim1=*nkeep**n;
                mycpyd(Wnew,&dim1,W);   
        }
        else{
                dim2sq=pow(*n,2);
                mycpyd(Wnew,&dim2sq,W);
        }       

        free(Wnew);
}       
/* }        */

free(X);   
}

