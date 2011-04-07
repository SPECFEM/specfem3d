/* 	Collection of string manipulation functions
	Hom Nath Gharti, NORSAR
	History:
	Apr 23,2010 (NORSAR)
	Mar 18, 2009 (Princeton University)
	Mar 13, 2008; Mar 19, 2008; HNG (NORSAR) */
#include <stdio.h>	
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void extractFileonly(char *filename, char *fileonly)
{
	int i, nchar, nslash;
	/* Extract file name without directory path */
	nchar=strlen(filename);
	nslash=0; /* Default value */
	for(i=0; i<nchar; i++){
		if(filename[nchar-i]=='/' || filename[nchar-i]=='\\'){ /* Integer epression! */
			nslash=nchar-i+1;
			break;
		}
	}

	strcpy(fileonly,&filename[nslash]);
}
/*======================================*/

/* 	This function returns the file name removing the extension if any
	In case of more than one '.' it will remove only last extension */
void removeExtension(char *filename, char *noextfile)
{
	int i, nchar, ndot;
	/* Extract file name without directory path */
	nchar=strlen(filename);/*printf("%s\n",filename);*//*printf("%s\n",noextfile);*/
	ndot=nchar; /* Default value */
	for(i=0; i<nchar; i++){
		if(filename[nchar-i]=='.'){ /* Integer epression! */
			ndot=nchar-i;
			break;
		}
	}
	/* printf("%d %d\n",nchar,ndot); */
	strncpy(noextfile,filename,ndot);
	noextfile[ndot]='\0'; /* This is usually requred for srncpy and srncat */
	/* printf("%s\n",noextfile); */
}
/*======================================*/

/* This function had been imported from e3d by shawn larsen */
void getFileName(filename, head, tail, number, max)
char *filename, *head, *tail; 
int number, max;
{
	int digits;
	char ext[10], format[20];
	
	/* printf("Hi\n"); */
	/* determine number of digits */
	digits = 1;
	if (max < 1) max = 1;
	while(max /= 10) digits++;

	/* determine format of number */
	/* sprintf(format, "%%0%dd\0", digits);
	following two lines are equivalent to previous line */
	sprintf(format, "%%0%dd", digits);
	strcat(format,"\0");
	sprintf(ext, format, number);

	/* get file name */
	strcpy(filename, head);
	strcat(filename,  ".");
	strcat(filename,  ext);
	strcat(filename,  ".");
	strcat(filename, tail);

}
/*======================================*/

/* This function had been imported from e3d by shawn larsen */
/* Function to open file name inputted through command line */
FILE *
getFile(argc, argv) 
int    argc;
char **argv;
{
	FILE *file;
	if (argc == 1) return(stdin);
	file = fopen(argv[1], "r");
	if (file == NULL) {
		fprintf(stderr, "Error: Can't open input file \"%s\"\n", argv[1]);
		exit(-1);
		}
	return(file);
}
/*======================================*/

/*	This function returns the position of last character of string s2 if s2 is found in string s1 */
int stringpos(char *s1, char *s2)
{
	int nchar1,nchar2,i1,i2,ipos,spos,stat;
	nchar1=strlen(s1);
	nchar2=strlen(s2);
	ipos=0;spos=0;
	if(nchar2>nchar1)return (spos);
	for (i1=0; i1<nchar1; i1++){
		if(i1+nchar2>nchar1)break;
		if(s1[i1]==s2[0]){
			ipos=i1;			
			stat=1;
			for (i2=0; i2<nchar2; i2++){
				if(s1[i2+ipos] != s2[i2]){
					stat=0;
					break;
				}
			}
			if(stat==1){
				spos=ipos+nchar2;
				break;
			}
		}
	}
	return(spos);
}
/*======================================*/

/* 	This function returns 1 if string s contains '#' as the first non-white space character
	otherwise 0 */ 
int commentline(char *s)
{
	int i,nchar,stat; /* stat = 1:yes, 0: No */
	nchar=strlen(s);
	stat=0; /* Default initialization to No */
	for(i=0; i<nchar; i++){
		if(s[i] != ' ' && s[i] != '\0' && s[i] != '\n' && s[i] != '\t'){
			if(s[i]=='#'){
				stat=1;
			}else{
				stat=0;
			}
			break;
		}
	}
	return(stat);
}
/*======================================*/

/* This function returns 1 if line s is blank, 0 if not blank */
int blankline(char *s)
{
	int i,nchar,stat; /* stat = 1:yes, 0: No */
	
	nchar=strlen(s);

	stat=1; /* Default is yes */
	for(i=0; i<nchar; i++){
		if(s[i] != ' ' && s[i] != '\t' && s[i] != '\n' && s[i] != '\0'){
			stat=0;
			break;
		}
	}	
	return(stat);
}

/* get integer value from the bulk string */
int get_int(int *var, char *arg, char *src)
{
int pos; /* position of matched string */

pos=stringpos(src,arg);

if(pos == 0){
  printf("ERROR: variable \"%s\" not found!\n",arg);
  exit(-1);
}

*var = atoi(&src[pos]); /* convert to integer value */
return(0);
}
/*======================================*/

/* look for integer value. if found return intger value and 
function value as 0 otherwise return -1 as a function value */
int look_int(int *var, char *arg, char *src)
{
int pos; /* position of matched string */

pos=stringpos(src,arg);

if(pos == 0){		
  return(-1);
}

*var = atoi(&src[pos]); /* convert to integer value */
return(0);
}
/*======================================*/

/* look for float value. if found return float value and 
function value as 0 otherwise return -1 as a function value */
int look_float(float *var, char *arg, char *src)
{
int pos; /* position of matched string */

pos=stringpos(src,arg);

if(pos == 0){		
  return(-1);
}

*var = atof(&src[pos]); /* convert to float value */
return(0);
}
/*======================================*/

/* look for double value. if found return double value and 
function value as 0 otherwise return -1 as a function value */
int look_double(double *var, char *arg, char *src)
{
int pos; /* position of matched string */

pos=stringpos(src,arg);

if(pos == 0){		
  return(-1);
}

*var = atof(&src[pos]); /* convert to float value */
return(0);
}
/*======================================*/

/* 	This function assigns the corresponding value immediately after the '=' or ':' sign following 
	the string arg to var, and exits the execution if no such arg is found in string s 

	May 16,2008,HNG: Now the argument name can be a part of other word in the line
	for eg vfile and file can not be problem!*/
#define nsymb 4
int getvalue(char *s, char *arg, char *type, int *var)
{
	int pos,inum; /* i,nchar,apos */
	double *dbl;
	float  *flt;
	char argt[10],*str,*symb[nsymb];
	
	symb[0]="="; /* x= */
	symb[1]=":"; /* x: */
	symb[2]=" ="; /* x = */
	symb[3]=" :"; /* x : */
	
	pos=0;inum=0;
	while(pos==0 && inum<nsymb){
		strcpy(argt,arg);
		strcat(argt,symb[inum]);
		pos=stringpos(s,argt);
		inum++;
	}
	
	if(pos == 0){
		printf("Variable \"%s\" not found!\n",arg);
		exit(-1);
	}
	
	/*nchar=strlen(s);
	for(i=pos; i<nchar; i++){
		if(s[i] == '=' || s[i] == ':'){
			apos=i+1;
			break;
		}
	}
	apos=pos;*/
	switch(type[0]){
		case 'd':
			*var = atoi(&s[pos]);
			return (0);
		case 'f':
			flt = (float *) var;
			*flt = atof(&s[pos]);
			return(0);
		case 'F':
			dbl = (double *) var;
			*dbl = atof(&s[pos]);
			return(0);
		case 's':
			str = (char *) var;
			strcpy(str, &s[pos]);
			return(0);
		default :
			printf("WARNING: bad argument type \"%s\" found!\n",type);
			exit(-1);
	}
}
/*======================================*/

/* 	This function assigns the corresponding value immediately after the '=' or ':' sign following 
	the string arg to var, and returns (-1) if no such arg is found in string s, otherwise this 
	function is exactly same as getvalue */
#define nsymb 4
int getvaluestat(char *s, char *arg, char *type, int *var)
{
	int pos,inum; /* i,nchar,apos */
	double *dbl;
	float  *flt;
	char argt[10],*str,*symb[nsymb];
	
	symb[0]="="; /* x= */
	symb[1]=":"; /* x: */
	symb[2]=" ="; /* x = */
	symb[3]=" :"; /* x : */
	
	pos=0;inum=0;
	while(pos==0 && inum<nsymb){
		strcpy(argt,arg);
		strcat(argt,symb[inum]);
		pos=stringpos(s,argt);
		inum++;
	}
	
	if(pos == 0){
		/* Variable not found */
		/* printf("Variable \"%s\" not found!\n",arg); */
		return(-1);
	}
	
	/*nchar=strlen(s);
	for(i=pos; i<nchar; i++){
		if(s[i] == '=' || s[i] == ':'){
			apos=i+1;
			break;
		}
	}
	apos=pos;*/
	switch(type[0]){
		case 'd':
			*var = atoi(&s[pos]);
			return (0);
		case 'f':
			flt = (float *) var;
			*flt = atof(&s[pos]);
			return(0);
		case 'F':
			dbl = (double *) var;
			*dbl = atof(&s[pos]);
			return(0);
		case 's':
			str = (char *) var;
			strcpy(str, &s[pos]);
			return(0);
		default :
			printf("WARNING: bad argument type \"%s\" found!\n",type);
			exit(-1);
	}
}
/*======================================*/

/*	This function returns the position of last character of string s2 in s1 if s2 is found in string s1 */
int matchfirstword(char *s1, char *s2)
{
	int nchar1,nchar2,i1,i2,ipos,stat,match;
	nchar1=strlen(s1);
	nchar2=strlen(s2);
	/*printf("%d %d\n",nchar1,nchar2);*/
	match=0;
	if(nchar2>nchar1)return(match); /* No match */
	
	for (i1=0; i1<nchar1; i1++){
		if(s1[i1] != ' ' && s1[i1] != '\0' && s1[i1] != '\n' && s1[i1] != '\t'){
			if(i1+nchar2>nchar1)break;
			if(s1[i1]==s2[0]){
				ipos=i1;			
				stat=1;
				for (i2=0; i2<nchar2; i2++){
					if(s1[i2+ipos] != s2[i2]){
						stat=0;
						break;
					}
				}
				if(stat==1)match=1;
			}
			break;
		}
	}
	return(match);
}
/*======================================*/

/* 	This function extracts and returns the string in the first quotes */
int getfirstquote(char *s1, char *s2)
{
	int i, nchar, ipos1,ipos2,nquote,num;
	char temp[62];

	nchar=strlen(s1);
	nquote=0; /* Default value */
	num=0;
	ipos1=ipos2=0;
	for(i=0; i<nchar; i++){
		if(s1[i]=='"'){ /* Integer exression! */
			nquote=i;
			num=num+1;
			if(num==1){
				ipos1=i;
			}
			else if(num==2){
				ipos2=i;
				break;
			}	
		}	
	}
	num=0;
	for(i=ipos1+1; i<ipos2;i++){
		temp[num]=s1[i];
		num=num+1;
	}
	
	temp[ipos2-1]='\0';
	strcpy(s2,temp);
	return(0);
}
/*======================================*/

/* 	This function extracts and returns the string after the given string */
int getintegervect(char *s1, char *arg, int n, int ivect[n])
{
	int i,nchar,slen;
	char temp[62];

	int pos; /* position of matched string */

    pos=stringpos(s1,arg);
   
	stringafterstring(s1,"=",temp);	
	for(i=0;i<n;i++){		
		slen=strlen(temp);
		if(sscanf(temp,"%d%n",&ivect[i],&nchar)!= 1){
			printf("ERROR: wrong values for -xmat option!\n");
			exit(-1);
		}		
		strncpy(temp,temp+nchar+1,slen);		
	}
	return(0);
}
/*======================================*/

/* 	This function extracts and returns the string after the given string */
int stringafterstring(char *s1, char *arg, char *s2)
{
	int i,nchar,num;
	char temp[62];

	int pos; /* position of matched string */

    pos=stringpos(s1,arg);

	nchar=strlen(s1);
	
	num=0;
	for(i=pos; i<nchar;i++){
		temp[num]=s1[i];
		num=num+1;
	}
	
	temp[nchar-pos]='\0';
	strcpy(s2,temp);
	return(0);
}
/*======================================*/

/* 	This function converts all the letters of string to lower case
	HNG,May 23,2008 */
int lowercase(char *s){
	int i,nchar;
	nchar=strlen(s);
	for(i=0;i<nchar;i++){
		s[i]=tolower(s[i]);
	}
	return(0);
}

/* 	This function converts all the letters of string to upper case
	HNG,May 23,2008 */
int uppercase(char *s){
	int i,nchar;
	nchar=strlen(s);
	for(i=0;i<nchar;i++){
		s[i]=toupper(s[i]);
	}
	return(0);
}
/*======================================*/

/* This function determines the byte order of the processor architecture
source: http://www.ibm.com/developerworks/aix/library/au-endianc/index.html?ca=drs-
Use a character pointer to the bytes of an int and then check its first byte to see if it is 0 or 1.
Mar 18,2009 (Princeton University) */
#define LE 0 /* Little Endian */
#define BE 1 /* Big Endian */

int getEndian() {
    int i = 1;
    char *p = (char *)&i;

    if (p[0] == 1)
        return LE;
    else
        return BE;
}
/*======================================*/


	
