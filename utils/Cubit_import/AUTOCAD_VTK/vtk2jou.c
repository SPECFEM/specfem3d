/* This program converts ASCII VTK file consisting of triangular/quadrilateral meshes into CUBIT Journal file

AUTHOR:
  Hom Nath Gharti
  NORSAR
  homnath_AT_norsar_DOT_no
DEPENDENCY:
  stringmanip.c: string manipulation routines
COMPILE:
  gcc vtk2jou.c -o vtk2jou
USAGE: 
  vtk2jou <inputfile>
  Example: vtk2jou ore_surface_mesh.vtk
HISTORY: 
  HNG,Dec 07,2010
-------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define OFF 0
#define ON 1

/* auxiliary routines */
void removeExtension(char *, char *);

/* main routine */
int main(int argc,char **argv){
int i,itmp,j,k;
int ndim;	/* geometry dimension */ 
int inode,nenode,nnode,nelmt; /* number of nodes, number of elements */
int point_stat,cell_stat; /* status */
int dumi,n1,n2,n3;
float x,y,z;

char token[62],dumc[10];
char fonly[62],infname[62],outfname[62];

FILE *inf,*outf;

if(argc<2){
  fprintf(stderr,"ERROR: input file not entered!\n");
  exit(-1);
}

printf("input file: %s\n",argv[1]);
printf("--------------------------------\n");
/* default input file name is argv[1]*/
strcpy(infname,argv[1]);
removeExtension(argv[1],fonly);

/* open input file */
inf=fopen(argv[1],"r");
if(inf==NULL){
  fprintf(stderr,"ERROR: file \"%s\" not found!",argv[1]);
  exit(-1);
}
/*printf("--------------------------------\n");*/

    
/* initialize some variables to 0 */
nnode=0; nelmt=0;

/* set default status to OFF */
point_stat=OFF; cell_stat=OFF;

sprintf(outfname,"%s.jou",fonly);			  
outf=fopen(outfname,"w");

point_stat=OFF;
cell_stat=OFF;

while(!feof(inf)){
  fscanf(inf,"%s",token);
  /* read dimensions */
  if(point_stat!=ON && strcmp(token,"POINTS")==0){
    printf("reading/writing POINTS...");
	fscanf(inf,"%d",&nnode);
	fscanf(inf,"%s\n",dumc);
	for(i=0;i<nnode;i++){
		fscanf(inf,"%f %f %f",&x,&y,&z);
		fprintf(outf,"create vertex %.4f %.4f %.4f\n",x,y,z);
	}
	printf("complete\n");
	continue;
  }

  if(cell_stat!=ON && strcmp(token,"CELLS")==0){
    printf("reading/writing CELLS...");
	fscanf(inf,"%d",&nelmt);
	fscanf(inf,"%d\n",&dumi);
	for(i=0;i<nelmt;i++){
		fprintf(outf,"create surface vertex ");
		fscanf(inf,"%d",&nenode);
		/* all nodes except last */
		for(j=0;j<nenode-1;j++){
			fscanf(inf,"%d",&inode);
			fprintf(outf,"%d ",inode+1); /* ParaView indices start from 0 whereas CUBIT indices start from 1!) */
		}
		/* Last node */
		fscanf(inf,"%d",&inode);
		fprintf(outf,"%d\n",inode+1);
	}
	printf("complete\n");
	fclose(inf);
	fclose(outf);
	break;
  }
}
printf("--------------------------------\n");
return(0);
}
/*======================================*/
