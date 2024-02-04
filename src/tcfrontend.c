/* Program to process the input for  tcenum  Coset Enumeration,
   using the Praeger-Soicher book format for Coxeter relators.
   Copyright L.H.Soicher, 1992-2024.

   This version expects ASCII input.

   The input presentation is read from tcfein, 
   error messages go to stderr, and the output for the 
   FORTRAN enumerator goes to tcfeout. 

   This program was heavily modified by LHS from the output of p2c
   applied to tcfronend4.p. */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define maxexpr  1000
#define maxword 10000

typedef enum {
  badmaxexpr, badmaxword, illchar, illgen, syntax 
} errortype;

static FILE *tcfein, *tcfeout;

static int num[256]; 
/* num[ch]=i if input generator character ch is generator number i, 
   with counting of generators starting from 1. Generator characters
   must be alphabetic (with upper and lower case distinguished). */

static int inv[105]; /* 105 == 2*(2*26)+1 */
/* inv[i]=j if the inverse of generator i is generator j */

static int instring(char *s, char ch)
{
/* Returns  1  (i.e. true)  if ch is a character in the string s,
   and  0  (i.e. false)  if not. */
return strchr(s,ch)!=NULL;
}

static void error(errortype err, char ch)
{
putchar('\n');
switch (err) {
   case badmaxexpr:
      fprintf(stderr,"*** constant maxexpr too small\n");
      break;
   case badmaxword:
      fprintf(stderr,"*** constant maxword too small\n");
      break;
   case illchar:
      fprintf(stderr,"*** illegal character '%c'\n", ch);
      break;
   case illgen:
      fprintf(stderr,"*** undeclared generator '%c'\n", ch);
      break;
   case syntax:
      fprintf(stderr,"*** syntax error - did not expect '%c'\n", ch);
      break;
   }
exit(EXIT_FAILURE);   /* error exit */
}

static void readexpr(FILE *f, char *e, char *valid, char *stop, char *ignore)
/* Reads an expression from file  *f  into the string  e,  up to and 
   including when a character in the string  stop  is read into  e
   and all the brackets are matched. 

   The valid characters are those in the the string  valid,  and the
   characters to be ignored are those in the string  ignore. */
{
int j;
int brackets;
char ch;
j = 0;
brackets = 0;   /* brackets == 0 iff all brackets are matched */
do {
   ch = getc(f);
   if (instring(ignore,ch)) 
      continue;
   if (!instring(valid,ch)) {
      /* ch is an invalid character */
      if (isalpha(ch))
         error(illgen, ch);
      else
         error(illchar, ch);
      }
   if (j > maxexpr)
       error(badmaxexpr, ' ');
   e[j] = ch;
   if (ch == '[' || ch == '(')
      brackets--;
   else if (ch == ']' || ch == ')')
      brackets++;
   j++; 
   } while (!instring(stop,ch) || brackets != 0);
}

static int value(char *e, int front, int *back)
{
/* Returns the value of the unsigned integer in e[front]...e[*back],
   where front is given, and e[*back+1] is the first non-digit in 
   e[front]... */
int mpr, val, j;
j = front;
while (isdigit(e[j]))
   j++;
*back = j - 1;
mpr = 1;
val = 0;
for (j = *back; j >= front; j--) {
   val += mpr*(e[j] - '0');
   if (j > front)
      mpr *= 10;
   }
return val;
}

static void invert(int *w, int front, int back)
{
/* inverts the word w[front]...w[back] */
int temp;
while (front <= back) {
    temp = w[front];
    w[front] = inv[w[back]];
    w[back] = inv[temp];
    front++;
    back--;
    }
}

static void power(int *w, int front, int *back, int n)
{
/* Where n is a non-negative integer, this function 
   puts (w[front]...w[*back])^n into the word w[front]... 
   and updates *back. */
int i, j, k;
if (n == 0) {
   *back = front - 1;
   return;
   }
k = *back;
for (i = 2; i <= n; i++) {
   for (j = front; j <= *back; j++) {
      k++;
      if (k > maxword)
         error(badmaxword, ' ');
      w[k] = w[j];
      }
   }
*back = k;
}

static void writeword(FILE *f, int *w, int front, int back)
/* Writes to file  *f  the word in w[front]...w[back] in the enum FORTRAN 
   coset enumerator format, starting with the length of the word, and
   then the word itself, with generators represented by integers 1,2,.... */
{
int j;
fprintf(f, "%d\n", back - front + 1);
for (j = front; j <= back; j++)
   fprintf(f, "%d\n", w[j]);
}

static void commutate(char *e, int *last, int *w, int front, int *back);
/* Prototype for function commutate. */ 

static void process(char *e, int *last, int *w, int front, int *back)
{
/* Calculates the word (as a sequence of generator numbers)
   defined by the expression  e[*last+1],...,e[stophere],  where
   e[stophere]  is the first occurrence in  e[*last+1],... of 
   "=", ".", ";", or ",", or is a right parenthesis ")" or "]" 
   matching a left parenthesis "(" or "[" in e[*last].

   The result is put in  w[front]...w[*back],  where front is
   given and *back is determined. 

   The function updates  *last  to  stophere. */ 

*back = front-1;
(*last)++;
if (instring("=.;,])",e[*last]))
   /* empty word */
   return;
if (e[*last] == '1') {
   process(e, last, w, front, back);
   return;
}
if (isalpha(e[*last])) {
   (*back)++;
   if (*back > maxword)
      error(badmaxword,' ');
   w[*back] = num[e[*last]];
   } 
else 
   if (e[*last] == '[' || e[*last] == '(') {
      process(e, last, w, front, back);
      if (e[*last] == ';' || e[*last] == ',')
         commutate(e, last, w, front, back);
      } 
   else
      error(syntax, e[*last]);
if (e[*last+1] == '-') {
   invert(w, front, *back);
   (*last)++;
}
if (isdigit(e[*last+1]))
   power(w, front, back, value(e, *last+1, last));
process(e, last, w, *back + 1, back);
}

static void commutate(char *e, int *last, int *w, int front, int *back)
{
/* Calculates the commutator [ w[front]...w[*back] , x ], 
   where  x  is the word or word sequence defined by  e[*last+1]...
   The function puts the result in  w,  and updates *last and *back. 
   
   Commutators are left normed, so that [a,b,c,...] means [[a,b],c...]. */

int i, backsave;
backsave = *back;
process(e, last, w, backsave + 1, back);
for (i = *back+1; i <= *back + (*back-front+1);  i++) {
   if (i > maxword)
      error(badmaxword, ' ');
   w[i]=w[i-(*back-front+1)];
   }
invert(w, front, backsave);
invert(w, backsave + 1, *back);
*back = *back + (*back-front+1);
if (e[*last] == ';' || e[*last] == ',')
   commutate(e, last, w, front, back);
}

static void mainproc(FILE *tcfein, FILE *tcfeout)
{
char e[maxexpr+1]; /* holds an input expression */
int w[maxword+1],ww[maxword+1]; /* words */
char ignore[256],valid[256],stop[256]; /* strings storing sets of characters */
char gens[2*26+1]; /* string listing the (alphabetic) input generators */
int a[53][53];   /* stores Coxeter rels */
int i, j, k, ngen, x, y, v, min, imin, jmin, back, last, firstback;
int flag;
/* read generators */
strcpy(valid,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
strcpy(ignore,",; \r\n\t");
strcat(valid,ignore); 
strcpy(stop,"."); 
strcat(valid,stop); 
readexpr(tcfein, e, valid, stop, ignore); 
j = 0;
while (e[j] != '.') {
   gens[j] = e[j]; 
   num[gens[j]] = j+1; /* since generators are indexed starting from 1 */
   j++;
   }
gens[j]='\0';  /* adding null terminator to gens */  
ngen = strlen(gens);
for (j = 1; j <= ngen; j++)
   inv[j] = j; /* this may change */
strcpy(valid,gens); 
strcat(valid,ignore); 
strcat(valid,stop); 
/* read non-involutions */
readexpr(tcfein, e, valid, stop, ignore);
j = 0;
while (e[j] != '.') {
   inv[num[e[j]]] = ngen+j+1;
   inv[ngen+j+1] = num[e[j]];
   j++;
   }
fprintf(tcfeout, "%d\n", ngen + j);
for (k = 1; k <= ngen + j; k++)
   fprintf(tcfeout, "%d\n", inv[k]);
/* read subgroup generators */
strcpy(valid,gens);
strcat(valid,"0123456789()[]-"); 
strcpy(ignore," \r\n\t+");
strcat(valid,ignore); 
strcpy(stop,",;."); 
strcat(valid,stop); 
do {
   readexpr(tcfein, e, valid, stop, ignore);
   last = -1;
   process(e, &last, w, 0, &back);
   writeword(tcfeout, w, 0, back);
   } while (e[last] != '.');
fprintf(tcfeout, "%d\n", -1);  /* to indicate start of relators */
/* check for Coxeter relators */
strcpy(valid,gens);
strcat(valid,"0123456789"); 
strcpy(ignore,",; \r\n\t");
strcat(valid,ignore); 
strcpy(stop,"."); 
strcat(valid,stop); 
readexpr(tcfein, e, valid, stop, ignore); 
if (e[0] != '.') {
   /* process Coxeter relators */
   for (i = 1; i <= ngen; i++) {
      for (j = 1; j <= ngen; j++)
         a[i][j] = 2;
      }
   /* Use the Praeger-Soicher book format for Coxeter relators. */
   j = 0;
   do {
      /* process a path in the Coxeter graph */
      if (!instring(gens,e[j]))
         error(syntax, e[j]);
      x = num[e[j]];
      while (isdigit(e[j+1])) {
         v = value(e, j+1, &j);
	 j++;
         if (!instring(gens,e[j]))
            error(syntax, e[j]);
	 y = num[e[j]];
         if (x==y)
            error(syntax, e[j]);
	 a[x][y] = v;
	 a[y][x] = v;
	 x = y;
         }
      j++;
      } while (e[j] != '.');
   min = 0;
   do {
      flag = 1;
      for (i = 1; i <= ngen; i++) {
	 for (j = i + 1; j <= ngen; j++) {
            if (a[i][j] > 0 && (a[i][j] < min || flag)) {
	       flag = 0;
	       imin = i;
	       jmin = j;
	       min = a[i][j];
               }  
            }
         }
      if (!flag) {
         a[imin][jmin] = 0;
	 fprintf(tcfeout, "%d\n", min*2);
	 for (j = 1; j<= min; j++) {
	    fprintf(tcfeout, "%d\n", imin);
	    fprintf(tcfeout, "%d\n", jmin);
            }
         }
      } while (!flag);
   }
strcpy(valid,gens);
strcat(valid,"0123456789()[]-"); 
strcpy(ignore," \r\n\t+");
strcat(valid,ignore); 
strcpy(stop,"=,;."); 
strcat(valid,stop); 
/* read other relators */
do {
   readexpr(tcfein, e, valid, stop, ignore);
   last = -1;
   process(e, &last, w, 0, &firstback);
    if (e[last] == '=')
       invert(w, 0, firstback);
    else
       writeword(tcfeout, w, 0, firstback);
    while (e[last] == '=') {
      readexpr(tcfein, e, valid, stop, ignore);
      last = -1;
      process(e, &last, ww, 0, &back);
      fprintf(tcfeout, "%d\n", firstback + back + 2);
      for (j = 0; j <= firstback; j++)
         fprintf(tcfeout, "%d\n", w[j]);
      for (j = 0; j <= back; j++)
         fprintf(tcfeout, "%d\n", ww[j]);
      }
   } while (e[last] != '.');
/* Now include relators of the form  x*x^{-1}  and  x^{-1}*x  for 
   each generator  x.  */
for (j = 1; j <= ngen; j++) {
   fprintf(tcfeout, "%d\n", 2);
   fprintf(tcfeout, "%d\n", j);
   fprintf(tcfeout, "%d\n", inv[j]);
   if (j != inv[j]) {
      fprintf(tcfeout, "%d\n", 2);
      fprintf(tcfeout, "%d\n", inv[j]);
      fprintf(tcfeout, "%d\n", j);
      }
   }
fprintf(tcfeout, "%d\n", -2); /* end of relators */
}

int main(int argc, char *argv[])
{  
if(argc != 3) {
   fprintf(stderr,"\n*** error: usage is 'tcfrontend <infile> <outfile>'\n");
   exit(EXIT_FAILURE);
   }
tcfein = fopen(argv[1],"r");
if (tcfein == NULL) {
   fprintf(stderr,"\n*** error opening file %s \n",argv[1]);
   exit(EXIT_FAILURE);
   }
tcfeout = fopen(argv[2],"w");
if (tcfeout == NULL) {
   fprintf(stderr,"\n*** error opening file %s \n",argv[2]);
   exit(EXIT_FAILURE);
   }
rewind(tcfein);
rewind(tcfeout);

mainproc(tcfein,tcfeout);

fclose(tcfein);
fclose(tcfeout);
exit(EXIT_SUCCESS); 
}
