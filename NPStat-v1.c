/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------- NPStat ---------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/* Code to extract pool statistics (theta, neutrality tests) 
from pileup data and a fasta file representing the outgroup sequence */

/* Compile with
gcc -O3 -o npstat NPStat-vXX.c -lgsl -lgslcblas -lm
substituting XX with version number
*/

/* Arguments: 
run "npstat" to see the arguments
*/ 

/* Output stats:
window number, bases (after filtering), bases with known outgroup allele (after filtering), 
average read depth, number of segregating sites, Watterson theta, Tajima's Pi, Tajima's D, Fay and Wu's H, divergence per base (from outgroup), HKA etc...
*/
  


/* Include libraries */ 

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>


/* Define substitutions */

#define PRINTSNPS 0
#define DEB(x) //x

#define print_comb 0 /* put 1 if you want to print combinatorial data, 0 otherwise. NOT IMPLEMENTED */
#define maxcheck 10000000 /* NOT IMPLEMENTED */
#define minimum_coverage 2 /* this fixes the minimum acceptable coverage in terms of reads aligned */
#define maximum_coverage 1000 /* this fixes the maximum acceptable coverage in terms of reads aligned */
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

/* Declare structures */

struct tests 
{
  unsigned long cov;
  unsigned long l;
  unsigned long l_out;
  unsigned long s;
  double num_t;
  double num_p;
  double num_hl;
  double num_hq;  
  double den_t;
  double den_p;
  double den_hl;
  double den_hq;  
};

struct fst_calc 
{
  unsigned long l;
  double gen_diff;
  double c_s;
};

struct combinatorial
{
  double *d_t;
  double *d_p;
  double *d_hl;
  double *d_hq;
};

struct combinatorial_fst
{
  double **c_s;
};

/* Declare functions */

int generate_pool_covariance_matrix(double *covmat, double *covmatpool, int na, int nb, int n);
int generate_covariance_matrix(double *covmat, int n);




int base_to_num(char base);

int read_line_pileup(FILE * bam_file, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int *  alt_base) ;

int extract_outgroup_base(FILE * fasta_out, unsigned long pos, unsigned long oldpos, int fasta_length);

void extract_stats(struct tests * test, struct combinatorial * comb, int n0, unsigned long n_ref, unsigned long n_alt_allele, unsigned long rd, unsigned long * n_alt, int ref_base, int alt_base, int out_base, int mb);

void extract_fst(struct fst_calc * fst, struct combinatorial_fst * combfst, int n01, unsigned long n_ref1, unsigned long n_alt_allele1, unsigned long rd1, unsigned long * n_alt_1, int ref_base1, int alt_base1, int n02, unsigned long n_ref2, unsigned long n_alt_allele2, unsigned long rd2, unsigned long * n_alt_2, int ref_base2, int alt_base2, int out_base, int mb);

//SNPS
void extract_snps(unsigned long pos, FILE * output_snps, unsigned long * n_alt_1, unsigned long * n_alt_2, int mb);

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*-------------------------- Functions -------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


int generate_covariance_matrix(double *covmat, int n) 
{
	double s,*a,*b,bb,bz,zz,a2,*xi_null2,*xi_b2;/*a[n+1],b[n]*/
	int i,j,n2,ca,cb;
	double theta; /*temporary*/
	gsl_matrix *c,*invc,*sigma,*ct,*invct;
	/*struct test_r0 *opt;*/
	
	if(n<2) return -10000;
	
	n2=ceil((double)n/2);
	//debug
	//assert(n2*2==n);
	n=(int)n2*2;
	
	a = (double *)calloc((unsigned long int)n+1,sizeof(double));
	b = (double *)calloc((unsigned long int)n,sizeof(double));
	
	a[0]=0;
	a2=0;
	for(i=2;i<=n+1;i++){
		a[i-1]=a[i-2]+1/(double)(i-1);
	}
	for(i=1;i<n;i++){
		a2+=1/(double)(i*i);
		b[i-1]=2*(double)n*(a[n]-a[i-1])/(double)((n-i+1)*(n-i))-2/(double)(n-i);
	}
	b[n-1]=0;
	sigma=gsl_matrix_alloc(n-1,n-1);
	gsl_matrix_set_zero(sigma);
	for (i=1;i<n;i++){
		if (2*i<n) { 
			gsl_matrix_set(sigma,i-1,i-1,b[i]); 
		}  else {
			if (2*i>n) { 
				gsl_matrix_set(sigma,i-1,i-1,b[i-1]-1/gsl_pow_2((double)i)); 
			}  else { 
				gsl_matrix_set(sigma,i-1,i-1,(a[n-1]-a[i-1])*2/(double)(n-i)-1/gsl_pow_2((double)i));
			}
		}
	}
	for (i=1;i<n;i++){
		for (j=1;j<i;j++){  
			if (i+j<n) { 
				gsl_matrix_set(sigma,i-1,j-1,(b[i]-b[i-1])/2); gsl_matrix_set(sigma,j-1,i-1,(b[i]-b[i-1])/2);
			}  else {
				if (i+j>n) { 
					gsl_matrix_set(sigma,i-1,j-1,(b[j-1]-b[j])/2-1/(double)(i*j));  
					gsl_matrix_set(sigma,j-1,i-1,(b[j-1]-b[j])/2-1/(double)(i*j)); 
				}  else { 
					gsl_matrix_set(sigma,i-1,j-1,(a[n-1]-a[i-1])/(double)(n-i)+(a[n-1]-a[j-1])/(double)(n-j)-(b[i-1]+b[j])/2-1/(double)(i*j));  
					gsl_matrix_set(sigma,j-1,i-1,(a[n-1]-a[i-1])/(double)(n-i)+(a[n-1]-a[j-1])/(double)(n-j)-(b[i-1]+b[j])/2-1/(double)(i*j)); 
				}
			}
		}
	}

        for(ca=1;ca<n;ca++){
		for(cb=1;cb<n;cb++){
			covmat[(cb-1)*(n-1)+ca-1]=gsl_matrix_get(sigma,ca-1,cb-1);
		}
	}


//gsl_ran_poisson_pdf(k,mu)

	free(a);
	free(b);
	
	gsl_matrix_free(sigma);
	return 1;
}



int base_to_num(char base)
{
  switch(base)
    {
    case 'A': return 1; break;
    case 'a': return 1; break;
    case 'C': return 2; break;
    case 'c': return 2; break;
    case 'G': return 3; break;
    case 'g': return 3; break;
    case 'T': return 4; break;
    case 't': return 4; break;
    default: return 0; 
    };
};


/*--------------------------------------------------------------*/
//int read_line_pileup(FILE * bam_file, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base)

//read_line_pileup(bam_file1, min_qual, min_mqual, &pos_base1, &n_ref1, &n_alt_allele1, &rd1, n_alt_1, &ref_base1, &alt_base1);
//int read_line_pileup(FILE * bam_file, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int *  alt_base) ;

int read_line_pileup(FILE * bam_file, unsigned long min_qual, unsigned long min_mqual, unsigned long * pos_base, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base) {

  int count_i; 
  char *cline, *cchrom, *cpileup, *cqual, *cmqual, crefbase;   
  unsigned long count_j, n_ins;
  unsigned long nseq;
  size_t nline;

  DEB(printf("entering read routine\n")); //debug

  cchrom=(char *) malloc(100);
  cpileup=(char *) malloc(40);
  cqual=(char *) malloc(20);
  cmqual=(char *) malloc(20); 
  cline=(char *) malloc(1);  

  nline=1;

  getdelim(&cline,&nline,9,bam_file);
  sscanf(cline,"%s\t",cchrom);	  
  getdelim(&cline,&nline,9,bam_file);
  sscanf(cline,"%lu\t",pos_base);	  
  getdelim(&cline,&nline,9,bam_file);
  sscanf(cline,"%c\t",&crefbase);	  
  getdelim(&cline,&nline,9,bam_file);
  sscanf(cline,"%lu\t",&nseq);	  
  getdelim(&cline,&nline,9,bam_file);
  if (strlen(cline)>=40) cpileup=(char *)realloc(cpileup,strlen(cline)+1);
  sscanf(cline,"%s\t",cpileup);	  
  //getdelim(&cline,&nline,9,bam_file);
  getdelim(&cline,&nline,10,bam_file);
  if (strlen(cline)>=20) cqual=(char *)realloc(cqual,strlen(cline)+1);
  ////if (strlen(cline)>=20) cmqual=(char *)realloc(cmqual,strlen(cline)+1);
  sscanf(cline,"%s\t",cqual);
  ////sscanf(cline,"%s\t",cmqual);	  
  //getdelim(&cline,&nline,10,bam_file);
  //if (strlen(cline)>=20) cmqual=(char *)realloc(cmqual,strlen(cline)+1);
  //sscanf(cline,"%s\n",cmqual);

  DEB(printf("read data %s\t%lu\t%c\t%lu\t%s\t%s\t%s\n",cchrom,*pos_base,crefbase,nseq,cpileup,cqual,cmqual)); //debug 

  //for (;(ct!="\n")&&(ct!=EOF);ct=fgetc(bam_file)) {};
  //ct=fgetc(bam_file); printf("%c",ct);
  //ct=fgetc(bam_file); printf("%c",ct);
  //ct=fgetc(bam_file); printf("%c\n",ct);

  for(count_i=0;count_i<5;count_i++){
    n_alt[count_i]=0;
  };
  count_j=0;
  *n_ref=0;
  *n_alt_allele=0;
  *n_tot_allele=0;

  if((nseq>=minimum_coverage)&&(nseq<=maximum_coverage)) {

  for(count_i=0;count_i<strlen(cpileup);count_i++){

    switch (cpileup[count_i]) {
      case '^': count_i++; break;
      case '$': break;
      case '*': count_j++; break;
      case '+': { 
	count_i++;
	for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
	  n_ins=n_ins*10+(cpileup[count_i]-48);  //printf(">"); //debug
	}; //printf("%lu ",n_ins);  //debug
	for(n_ins=n_ins-1;n_ins>0;n_ins--) {
	  count_i++; //printf("<"); printf("%lu ",n_ins); if(count_i>20) break; //debug
	}; //printf("%c ",cpileup[count_i]); //debug
      }; break;
      case '-': { 
	count_i++;
	for(n_ins=0;isdigit(cpileup[count_i])!=0;count_i++){
	  n_ins=n_ins*10+(cpileup[count_i]-48);   
	};
	for(n_ins=n_ins-1;n_ins>0;n_ins=n_ins-1) {
	  count_i++;
	}; 
      }; break;
      case 'N': count_j++; break;
      case 'n': count_j++; break;
      case '.': {if((cqual[count_j]>=min_qual+33)) (n_alt[0])++; count_j++;}; break;
      case ',': {if((cqual[count_j]>=min_qual+33)) (n_alt[0])++; count_j++;}; break;
      case 'A': {if((cqual[count_j]>=min_qual+33)) (n_alt[1])++; count_j++;}; break;
      case 'C': {if((cqual[count_j]>=min_qual+33)) (n_alt[2])++; count_j++;}; break;
      case 'G': {if((cqual[count_j]>=min_qual+33)) (n_alt[3])++; count_j++;}; break;
      case 'T': {if((cqual[count_j]>=min_qual+33)) (n_alt[4])++; count_j++;}; break;
      case 'a': {if((cqual[count_j]>=min_qual+33)) (n_alt[1])++; count_j++;}; break;
      case 'c': {if((cqual[count_j]>=min_qual+33)) (n_alt[2])++; count_j++;}; break;
      case 'g': {if((cqual[count_j]>=min_qual+33)) (n_alt[3])++; count_j++;}; break;
      case 't': {if((cqual[count_j]>=min_qual+33)) (n_alt[4])++; count_j++;}; break;
      }; 
      

  };
  
  //  *ref_base=base_to_num(crefbase);
  *n_alt_allele = 0;
  *n_ref=0; //n_alt[0];
  for(count_i=0;count_i<5;count_i++){
    if (*n_ref < n_alt[count_i]) {
      *n_ref = n_alt[count_i]; *ref_base=count_i;
    };
  };
  for(count_i=0;count_i<5;count_i++){
    if ((*n_alt_allele < n_alt[count_i])&&(count_i!=*ref_base)) {
      *n_alt_allele = n_alt[count_i]; *alt_base=count_i;
    };
  };  
  if (*ref_base==0) { *ref_base=base_to_num(crefbase); };
  if (*alt_base==0) { *alt_base=base_to_num(crefbase); };  
  *n_tot_allele = *n_ref + *n_alt_allele;


  };
  
  DEB(printf("using data %lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%u\t%u\n", *n_ref, *n_alt_allele, *n_tot_allele, n_alt[0], n_alt[1],n_alt[2],n_alt[3],n_alt[4], *ref_base, *alt_base)); //debug

  DEB(printf("exit read routine\n")); //debug
  free(cchrom);
  free(cpileup);
  free(cqual);
  free(cmqual); 
  free(cline);  

};

















/* Gipo's input and ms2pileup functions
void read_line_ms(FILE *bam_file, unsigned long w_s, int n03, long int * n_pos, unsigned long int int_positions[], int array_counts[])
{
  int i = 0;
  char *line, * word, * pEnd;;
  line=(char *)malloc(10000);
  size_t nline;
  nline = 1;
  // Four lines that won't be used
  for (i=1; i<=4; i++)
    {
      getline(&line, &nline, bam_file);
    }
  // Line n° 5, number of positions
  getline(&line, &nline, bam_file);
  word = strtok(line, " ");
  word = strtok(NULL, " ");
  *n_pos = strtol(word, &pEnd, 10);
  // The sixth line, that with positions
  getline(&line, &nline, bam_file);
  word = strtok(line, " ");
  for (i=0; i<*n_pos; i++)
    {
      word = strtok(NULL, " ");
      int_positions[i] = ceil(w_s * strtof(word, &pEnd));
      if (i>0)
	{
	  while (int_positions[i] <= int_positions[i - 1])
	    {
	      (int_positions[i])++;
	    }
	}
    }
  // Variation matrix
  int var_line, column;
  for (var_line=1; var_line<=n03; var_line++)
    {
      getline(&line, &nline, bam_file);
      char character;
      for (column=1; column<=*n_pos; column++)
	{
	  character = line[column - 1];
	  if (character == '1')
	    {
	      (array_counts[column - 1])++;
	    }
	}
    }
}

void base_repetition(unsigned long w_s, int n_ind, unsigned long * pos_base, long int * n_pos, unsigned long int int_positions[], int *n_variation, int array_counts[], double input_coverage, double error_rate, const gsl_rng * s, unsigned long * n_ref, unsigned long * n_alt_allele, unsigned long * n_tot_allele, unsigned long * n_alt, int * ref_base, int * alt_base)
{
  unsigned long int cur_phys_pos= (unsigned long int)(*pos_base)+1;
  unsigned long int pos_with_variation = int_positions[*n_variation - 1];
  int * cur_array_counts = &array_counts[*n_variation - 1];
//  for (cur_phys_pos=1; cur_phys_pos<=w_s; cur_phys_pos++)
//    {
      if (cur_phys_pos == pos_with_variation)
	{
	  *n_tot_allele = gsl_ran_poisson(s, input_coverage);
	  double ones_mean = (double) (*cur_array_counts) / (double) (n_ind);
	  *n_alt_allele = gsl_ran_binomial (s, ones_mean, *n_tot_allele);
	  *n_ref = *n_tot_allele - *n_alt_allele;
	  *ref_base = 1;
	  *alt_base = 2;
	  n_alt[0] = 0;
	  n_alt[1] = *n_ref;
	  n_alt[2] = *n_alt_allele;
	  n_alt[3] = 0;
	  n_alt[4] = 0;
	  cur_array_counts++;
	  (*n_variation)++; 
	}
      else
	{
	  *n_tot_allele = gsl_ran_poisson(s, input_coverage);
	  double argument = (double) (*n_tot_allele) * error_rate;
	  *n_alt_allele = min(gsl_ran_poisson(s, argument), *n_tot_allele);
	  *n_ref = *n_tot_allele - *n_alt_allele;
	  *ref_base = 1;
	  *alt_base = 2;
	  n_alt[0] = 0;
	  n_alt[1] = *n_ref;
	  n_alt[2] = *n_alt_allele;
	  n_alt[3] = 0;
	  n_alt[4] = 0;
	}
//    }
(*pos_base)++;
}
*/

/*--------------------------------------------------------------*/


int extract_outgroup_base(FILE * fasta_out, unsigned long pos, unsigned long oldpos, int fasta_length)
{  
  unsigned long diff;
  ldiv_t pos_newlines, oldpos_newlines;
  pos_newlines=ldiv(pos-1,(unsigned long)fasta_length);
  oldpos_newlines=ldiv(oldpos-1,(unsigned long)fasta_length);
  diff=pos-oldpos+pos_newlines.quot-oldpos_newlines.quot;
  fseek(fasta_out,diff-1,SEEK_CUR);
  return base_to_num(fgetc(fasta_out));
};




/*
int dmau_extract_outgroup_base(FILE * fasta_out, unsigned long pos, unsigned long oldpos, unsigned long * fasta_length, char ref_base)
{ 
  char *line;
  size_t nline;
  char ref;
  char out;
  long fpos;
  nline=1000;
  line=(char *)malloc(nline*sizeof(char));
  fpos=ftell(fasta_out);
  getline(&line,&nline,fasta_out);
  sscanf(line,"%lu %s %s\n",fasta_length,&ref,&out);
  
  if(pos<(*fasta_length)){
    fseek(fasta_out, fpos, SEEK_SET);
    return ref_base;
  } else {
    return base_to_num(out);
  }
};

*/





/*--------------------------------------------------------------*/


void extract_stats(struct tests * test, struct combinatorial * comb, int n0, unsigned long n_ref, unsigned long n_alt_allele, unsigned long rd, unsigned long * n_alt, int ref_base, int alt_base, int out_base, int mb)
{  
  char is_out=0;
  
  if (out_base>0)
    {
      if (ref_base==out_base) { is_out=1; } 
      else if (n_alt[out_base]==n_alt_allele) { is_out=2; };
        // else is_out=3;
    };
  
  test->cov+=rd;
  test->l+=1;
  test->den_t+=comb->d_t[rd-1];
  test->den_p+=comb->d_p[rd-1];
  if (is_out>0)
    {
      test->l_out+=1;
      test->den_hl+=comb->d_hl[rd-1];
      test->den_hq+=comb->d_hq[rd-1];     
    };
  // else printf("outgroup bases %u %u %u\n",is_out,out_base,ref_base);//debug
  
  if ((n_ref>mb)&&(n_ref<rd-mb))
    {
    test->s+=1;
    test->num_t+=1;
    test->num_p+=(double)(2*n_alt_allele*n_ref)/(double)(rd*(rd-1));
    if (is_out==1){   
      test->num_hl+=(double)(n_alt_allele)/(double)(rd-1);
      test->num_hq+=(double)(n_alt_allele*n_alt_allele)/(double)(rd*(rd-1));
    };
    if (is_out==2)
      {   
	test->num_hl+=(double)(n_ref)/(double)(rd-1);
	test->num_hq+=(double)(n_ref*n_ref)/(double)(rd*(rd-1));
      };
    };
  
  
  DEB(printf("using data %u\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%u\t%u\t%u\n", n0, n_ref, n_alt_allele, rd, n_alt[0], n_alt[1],n_alt[2],n_alt[3],n_alt[4], ref_base, alt_base, out_base)); //debug
  
};



/*--------------------------------------------------------------*/

unsigned long extract_pos_snpinput( FILE * snpinput){
  char *line;
  int c,i,j;
  size_t l_line;
  unsigned long pos;
  line = (char *) malloc (100*sizeof(char));
  l_line=100;
  char test; test=fgetc(snpinput); if (test!=EOF) {
  ungetc(test, snpinput); 
  //getline(&line,&l_line,snpinput);
  //for(i=1;i<=column;i++){
  getline(&line,&l_line,snpinput);
  //};
  sscanf(line,"%lu",&pos);//sscanf(line,"%lu\t",&pos); // or sscanf(line,"%lu ",pos);
  //i=0;
  //for(c=1;c<column;c++){
  //  for(;(i<l_line)&&((line[i]==9)||(line[i]==32))&&((line[i+1]!=9)&&(line[i+1]!=32));i++){};    
  //};
  //j=i;
  //for(;(i<l_line)&&((line[i]!=9)&&(line[i]!=32))&&((line[i+1]==9)||(line[i+1]==32));i++){};    
  //getline(&line,&l_line,snpinput);
  } else { pos=2000000000; };
  free(line); 
  return pos;
};  


/*
void extract_fst(struct fst_calc * fst, struct combinatorial_fst * combfst, int n01, unsigned long n_ref1, unsigned long n_alt_allele1, unsigned long rd1, unsigned long * n_alt_1, int ref_base1, int alt_base1, int n02, unsigned long n_ref2, unsigned long n_alt_allele2, unsigned long rd2, unsigned long * n_alt_2, int ref_base2, int alt_base2, int out_base, int mb)
{  
  int i,j;
  
  fst->l+=1;
  for (i=1;i<5;i++)
    {
      for (j=0;j<i;j++)
	{
	  fst->gen_diff+=(double)(n_alt_1[i]*n_alt_2[j])/(double)(rd1*rd2);
	};
    };
  fst->c_s+=combfst->c_s[rd1-1][rd2-1];  
};
*/

//SNPS
void extract_snps(unsigned long pos, FILE * output_snps, unsigned long * n_alt_1, unsigned long * n_alt_2, int mb)
{
  int i, c1, c2;
  
  c2=0;
  c1=0;
  for(i=1;i<5;i++)
    {
      if ((n_alt_1[i]+n_alt_2[i])>(n_alt_1[c1]+n_alt_2[c1]))
	{
	  c1=i;
	};
    };
  if (c1==0) c2=1;
  for(i=c2+1;i<5;i++)
    {
      if (((n_alt_1[i]+n_alt_2[i])>(n_alt_1[c2]+n_alt_2[c2]))&&(i!=c1))
	{
	  c2=i;
	};
    };
  if((n_alt_1[c2]+n_alt_2[c2])>1)
    {
      fprintf(output_snps, "%lu\t%lu\t%lu\t%lu\t%lu\n", pos, n_alt_1[c2], n_alt_2[c2], n_alt_1[c1], n_alt_2[c1]);
    };
};

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*------------------------ Main program ------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/



int main(int argc, char *argv[])
{
  // arguments: bam file 1, bam file 2, fasta outgroup file, window size, haploid sample size, minimum coverage, maximum coverage, minimum base quality, minimum mapping quality, low frequency alleles removed 

  /* Input flags */
  double input_coverage1;
  /* Variables */
  
  int count_i,count_j;
  FILE *bam_file, *bam_file1, *bam_file2, *fasta_ref, *fasta_out;
  //char outgroup_available;
    char ct, ct1, ct2;
    int lt;
  unsigned long line, line1, line2, pos, oldpos, pos_base, pos_base1, pos_base2, nseq;
  unsigned long n_ref1, n_alt_allele1, n_tot_allele1; 
  unsigned long n_alt_1[5];
  unsigned long n_ref2, n_alt_allele2, n_tot_allele2; 
  unsigned long n_alt_2[5];
  unsigned long max_cov, min_cov, min_qual, min_mqual;
  int m_bar, m1t, m2t;
  unsigned long window_size;
  int n0max, n01, n02;
  unsigned long rd, rd1, rd2;
  FILE *output_stat1, *output_stat2, *output_fst;
  char *temp_file_name1, *temp_file_name2;
  //unsigned long 
  int fasta_length;
  int ref_base, ref_base1, ref_base2, alt_base1, alt_base2, out_base;
  /* Gipo's variables */
  long int n_pos1;
  int array_counts1[2000] = {0};
  unsigned long int int_positions1[1000] = {0};
  int n_variation1;
  n_variation1=1;
  double error_rate1 = 0.001;
  
  // Random number generator
  const gsl_rng_type * type;
  gsl_rng * r;
  gsl_rng_env_setup();
  type = gsl_rng_default;
  r = gsl_rng_alloc (type);
  gsl_rng_env_setup();
  
  unsigned long n_window;
  unsigned long div;
  
  unsigned long *vec_rd; 
  double *vec_s, *vec_p, *vec_h, *vec0_s, *vec0_d, *vec0_h;
  double *covmat; 
  
  struct combinatorial comb1, comb2;
  struct combinatorial_fst combfst;
  struct tests test1, test2;
  struct fst_calc fst;
  
  FILE * list_snps;
  unsigned long pos_snp;
  
  FILE * gff;
  char if_gff;
  unsigned long cds_start, cds_end;
  int phase_cds, frame; 
  char strand, feature[256], *gff_file_name;
  unsigned long psyn, pnon, dsyn, dnon;
  
  
  
  /* Main part of the program */
  
  /* Read arguments */
  window_size=0;
  n01=0;
  n02=0;
  min_cov=4;
  max_cov=100;
  min_qual=10;
  min_mqual=10;
  m_bar=1;
  
  char from_stdin, outgroup_available, compute_fst, ext_snps, map_qual_pileup;
  char *pileup_file_name, *pileup2_file_name, *outgroup_file_name, *snp_file_name;
  int arg_i, col_chrom, col_pos;
  from_stdin=2;
  outgroup_available=0;
  compute_fst=0;
  ext_snps=0;
  map_qual_pileup=0;
  if_gff=0;
  
  for(arg_i=1;arg_i<argc-1;arg_i++)
    {
      if (strcmp(argv[arg_i], "-n") == 0) {arg_i++; sscanf(argv[arg_i], "%u", &n01); }
      else if (strcmp(argv[arg_i], "-l") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &window_size); }
      else if (strcmp(argv[arg_i], "-cov") == 0) {arg_i++; sscanf(argv[arg_i], "%lf", &input_coverage1); }
//      else if (strcmp(argv[arg_i], "-n2") == 0) {arg_i++; sscanf(argv[arg_i], "%u", &n02); }
      else if (strcmp(argv[arg_i], "-mincov") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &min_cov); }
      else if (strcmp(argv[arg_i], "-maxcov") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &max_cov); }
      else if (strcmp(argv[arg_i], "-minqual") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &min_qual); }
//      else if (strcmp(argv[arg_i], "-mapqual") == 0) {map_qual_pileup=1; }
//      else if (strcmp(argv[arg_i], "-minmapqual") == 0) {arg_i++; sscanf(argv[arg_i], "%lu", &min_mqual); }
      else if (strcmp(argv[arg_i], "-nolowfreq") == 0) {arg_i++; sscanf(argv[arg_i], "%u", &m_bar); }
      else if (strcmp(argv[arg_i], "-outgroup") == 0) {outgroup_available=1; arg_i++; outgroup_file_name=argv[arg_i];}
//      else if (strcmp(argv[arg_i], "-fstpop2") == 0) {compute_fst=1; arg_i++; pileup2_file_name=argv[arg_i];}
      //else if (strcmp(argv[arg_i], "-pileup") == 0) {}
      else if (strcmp(argv[arg_i], "-snpfile") == 0) {ext_snps=1; arg_i++; snp_file_name=argv[arg_i];}
      else if (strcmp(argv[arg_i], "-annot") == 0) {if_gff=1; arg_i++; gff_file_name=argv[arg_i];}
/*      else if (strcmp(argv[arg_i], "-snpcolumns") == 0) 
	{
	  arg_i++; sscanf(argv[arg_i], "%u", &col_chrom); 
	  arg_i++; sscanf(argv[arg_i], "%u", &col_pos);
	};
*/    }
  if (arg_i==argc-1) {
    if (strcmp(argv[arg_i], "-") == 0) from_stdin=1; else {from_stdin=0; pileup_file_name=argv[arg_i];}
  };
  //from_stdin=0; pileup_file_name=argv[arg_i-1];
  
  if ((n01==0)||(window_size==0)||(from_stdin==2))
    {
      fprintf(stderr,"Missing values in command line!\n  Command:\n    npstat [options] file.pileup\n  or to read from standard input:\n    npstat [options] -\n  Options:\n   -n samplesize : haploid sample size\n   -l windowlength : window length\n   -mincov minimum_coverage : filter on minimum coverage (default 4)\n   -maxcov maximum_coverage : filter on maximum coverage (default 100)\n   -minqual minimum_base_quality : filter on base quality (default 10)\n   -nolowfreq m : filter on minimum allele count mac>m\n   -outgroup file.fa : outgroup file in FASTA\n   -annot file.gff3 : annotation file in GFF3\n   -snpfile file.snp : consider SNPs only if present in file.snp\n"); 
      return(-1);
    };
//        fprintf(stderr,"Missing values in command line!\n  Command:\n    NPStat [options] file.pileup\n  or to read from standard input:\n    NPStat [options] -\n  Options:\n   -n samplesize : haploid sample size\n   -l windowlength : window length\n   -mapqual : pileup includes mapping quality\n   -mincov minimum_coverage : filter on minimum coverage (default 4)\n   -maxcov maximum_coverage : filter on maximum coverage (default 100)\n   -minqual minimum_base_quality : filter on base quality (default 10)\n   -minmapqual minimum_mapping_quality : filter on mapping quality (default 10)\n   -nolowfreq m : filter on minimum allele count mac>m\n   -outgroup file.fa : outgroup file in FASTA\n   -fstpop2 file2.pileup : computes Fst with a second population \n    contained in file2.pileup\n   -n2 : sample size of the second population\n"); 
  //-snpfile file.snp : consider SNPs only if present in file.snp\n
  //-snpcolumns chrom pos : columns for chromosome and \n
  //position in file.snp (default values 1, 2)\n

  ////////outgroup_available=1;
  
  /*outgroup_available=0;
    if (argc<3)
    {
    fprintf(stderr, "Usage: bam_read_routine file.bam reference.fasta outgroup.fasta") ;
    return 1;
    }
    else if (argc>3)
    { 
    outgroup_available=1; 
    };*/
  
  /* Initialize and open file */
  
  //temp_file_name1=malloc(strlen(pileup_file_name)+strlen(pileup2_file_name)+10);
  //temp_file_name2=malloc(strlen(pileup_file_name)+strlen(pileup2_file_name)+10);
  temp_file_name1=malloc(strlen(pileup_file_name)+strlen(pileup_file_name)+10);
  temp_file_name2=malloc(strlen(pileup_file_name)+strlen(pileup_file_name)+10);
  
  if (from_stdin==0) {
    bam_file1=fopen(pileup_file_name,"r");
    if (bam_file1 == NULL){
      fprintf(stderr,"Error: the pileup file cannot be opened!");
      return(-1);
    };     
  } else {
    bam_file1=stdin;
    pileup_file_name="NPStat_input";
  };
  //bam_file2=fopen(pileup2_file_name,"r");
  
  if(ext_snps==1){
  list_snps=fopen(snp_file_name,"r");
  if (list_snps == NULL){
      fprintf(stderr,"Error: the SNP list file cannot be opened!");
      return(-1);
  };
  };
  
  if(if_gff==1){
  gff=fopen(gff_file_name,"r");
  if (gff == NULL){
      fprintf(stderr,"Error: the GFF3 file cannot be opened!");
      return(-1);
  };
  };
  
  strcpy(temp_file_name1,pileup_file_name);
  strcpy(temp_file_name2,".stats");
  strcat(temp_file_name1,temp_file_name2);
  output_stat1=fopen(temp_file_name1,"w");
  if (output_stat1 == NULL){
      fprintf(stderr,"Error: the output file cannot be written!");
      return(-1);
  };
  /*strcpy(temp_file_name1,pileup2_file_name);
  strcpy(temp_file_name2,",stats");
  strcat(temp_file_name1,temp_file_name2);
  output_stat2=fopen(temp_file_name1,"w"); 
  
  strcpy(temp_file_name1,"fst_");
  strcpy(temp_file_name2,pileup_file_name);
  strcat(temp_file_name1,temp_file_name2);
  strcpy(temp_file_name2,"_");
  strcat(temp_file_name1,temp_file_name2);
  strcpy(temp_file_name2,pileup2_file_name);
  strcat(temp_file_name1,temp_file_name2);
  //output_fst=fopen(temp_file_name1,"w"); */
  
  //SNPS
/*  strcpy(temp_file_name1,"snps_");
  strcpy(temp_file_name2,pileup_file_name);
  strcat(temp_file_name1,temp_file_name2);
  strcpy(temp_file_name2,"_");
  strcat(temp_file_name1,temp_file_name2);
  strcpy(temp_file_name2,pileup2_file_name);
  strcat(temp_file_name1,temp_file_name2);
#if PRINTSNPS == 1
  FILE *output_snps;
  output_snps=fopen(temp_file_name1,"w");
#endif
*/

  /*bam_file=fopen(argv[1],"r");
    fasta_ref=fopen(argv[2],"r");
    if (outgroup_available==1) {
    fasta_out=fopen(argv[3],"r");
    };*/
  /*for (;iscntrl(fgetc(fasta_ref))==0;) {};*/
  fasta_length=0;
  if (outgroup_available==1)
    { fasta_length=0;
      fasta_out=fopen(outgroup_file_name,"r");
      for (;iscntrl(fgetc(fasta_out))==0;) 
	{
	};
      // fasta_pos=ftell(fasta_out);
      for (count_i=0;iscntrl(fgetc(fasta_out))==0;count_i++)
	{
	};
      fasta_length=count_i;
      rewind(fasta_out);   
      for (;iscntrl(fgetc(fasta_out))==0;) 
	{
	};
    };
  //for (;ct!="\n";ct=fgetc(fasta_ref)) {};
  //for (;ct!="\n";ct=fgetc(fasta_out)) {};
  
  /* Load combinatorics */
  printf("Initializing combinatorics...\n");

  vec_rd=(unsigned long *)malloc(max_cov*sizeof(unsigned long));
  vec_s=(double *)malloc(max_cov*(n01-1)*sizeof(double));
  vec_p=(double *)malloc(max_cov*(n01-1)*sizeof(double));
  vec_h=(double *)malloc(max_cov*(n01-1)*sizeof(double));
  vec0_s=(double *)malloc(max_cov*sizeof(double));
  vec0_d=(double *)malloc(max_cov*sizeof(double));
  vec0_h=(double *)malloc(max_cov*sizeof(double));
  
  comb1.d_t=(double *)malloc(max_cov*sizeof(double));
  comb1.d_p=(double *)malloc(max_cov*sizeof(double));
  comb1.d_hl=(double *)malloc(max_cov*sizeof(double));
  comb1.d_hq=(double *)malloc(max_cov*sizeof(double));
  /*comb2.d_t=(double *)malloc(max_cov*sizeof(double));
  comb2.d_p=(double *)malloc(max_cov*sizeof(double));
  comb2.d_hl=(double *)malloc(max_cov*sizeof(double));
  comb2.d_hq=(double *)malloc(max_cov*sizeof(double));
  combfst.c_s=(double **)malloc(max_cov*sizeof(double *));
  for(count_i=0;count_i<max_cov;count_i++)
    {
      combfst.c_s[count_i]=(double *)malloc(max_cov*sizeof(double));
    };*/
  for(count_i=0;count_i<max_cov;count_i++)
    {
      comb1.d_t[count_i]=0;
      comb1.d_p[count_i]=0;
      comb1.d_hl[count_i]=0;
      comb1.d_hq[count_i]=0;
/*      comb2.d_t[count_i]=0;
      comb2.d_p[count_i]=0;
      comb2.d_hl[count_i]=0;
      comb2.d_hq[count_i]=0;
      for(count_j=0;count_j<max_cov;count_j++)
	{
	  combfst.c_s[count_i][count_j]=0;
	};*/
    };
  for(rd=2*m_bar+1;rd<=max_cov;rd++)
    {
      int i,j,k;
      comb1.d_p[rd-1]+=((double)n01-1)/(double)n01;
      //comb2.d_p[rd-1]+=((double)n02-1)/(double)n02;
      comb1.d_hl[rd-1]+=(double)rd*(n01-1)/(double)(n01*(rd-1));
      //comb2.d_hl[rd-1]+=(double)rd*(n02-1)/(double)(n02*(rd-1));
      comb1.d_hq[rd-1]+=(double)rd*((double)((n01-1)*(rd+1))/(double)(2*rd))/(double)(n01*(rd-1));
      //comb2.d_hq[rd-1]+=(double)rd*((double)((n02-1)*(rd+1))/(double)(2*rd))/(double)(n02*(rd-1));  
      for(k=1;k<n01;k++)
	{
	  comb1.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n01,rd)-gsl_pow_int(1-(double)k/(double)n01,rd))/(double)k;
	  for(lt=1;lt<=m_bar;lt++)
	    {
	      comb1.d_t[rd-1]+=-gsl_sf_choose(rd,lt)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt-1)/(double)n01;
	      comb1.d_p[rd-1]+=-2*gsl_sf_choose(rd-2,lt-1)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt-1)/(double)n01;
	      comb1.d_hq[rd-1]+=-1/(double)(n01*(rd-1))*((double)lt*gsl_sf_choose(rd-1,lt-1)*gsl_pow_int((double)k/(double)n01,lt-1)*gsl_pow_int(1-(double)k/(double)n01,rd-lt) +(double)(rd-lt)*gsl_sf_choose(rd-1,rd-lt-1)*gsl_pow_int((double)k/(double)n01,rd-lt-1)*gsl_pow_int(1-(double)k/(double)n01,lt));
	    };
	  comb1.d_hl[rd-1]+=((double)rd/(double)n01)*((1-2/((double)rd-1))*(double)k/(double)n01-1)
	    *gsl_pow_int((double)k/(double)n01,rd-2);
	  comb1.d_hq[rd-1]+=-(double)rd/(double)(n01*(rd-1))*(gsl_pow_int((double)k/(double)n01,rd-1));
	  //(((double)rd-1)/(double)n01)*((1+(double)(rd+1)/gsl_pow_int((double)(rd-1),2))*(double)k/(double)n01-1)*gsl_pow_int((double)k/(double)n01,rd-2);
	};
      /*for(k=1;k<n02;k++)
	{
	  comb2.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n02,rd)-gsl_pow_int(1-(double)k/(double)n02,rd))/(double)k;
	  for(lt=1;lt<=m_bar;lt++)
	    {
	      comb2.d_t[rd-1]+=gsl_sf_choose(rd,lt)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt-1)/(double)n02;
	      comb2.d_p[rd-1]+=2*gsl_sf_choose(rd-2,lt-1)*gsl_pow_int((double)k/(double)n02,lt-1)*gsl_pow_int(1-(double)k/(double)n02,rd-lt-1);
	    };
	  comb2.d_hl[rd-1]+=((double)rd/(double)n02)*((1-2/((double)rd-1))*(double)k/(double)n02-1)
	    *gsl_pow_int((double)k/(double)n02,rd-2);
	  comb2.d_hq[rd-1]+=(((double)rd-1)/(double)n02)*((1+(double)(rd+1)/gsl_pow_int((double)(rd-1),2))*(double)k/(double)n02-1)*gsl_pow_int((double)k/(double)n02,rd-2);
	  /*    comb2.d_t[rd-1]+=(1-gsl_pow_int((double)k/(double)n02,rd-1)*((double)k/(double)n02+(double)rd)-gsl_pow_int(1-(double)k/(double)n02,rd))/(double)k;
		comb2.d_p[rd-1]+=(-2*gsl_pow_int((double)k/(double)n02,rd-2))/(double)n02;
		comb2.d_hl[rd-1]+=((double)rd/(double)n02)*((1-2/((double)rd-1))*(double)k/(double)n02-1)
		*gsl_pow_int((double)k/(double)n02,rd-2);
		comb2.d_hq[rd-1]+=(((double)rd-1)/(double)n02)*((1+(double)(rd+1)/gsl_pow_int((double)(rd-1),2))*(double)k/(double)n02-1)*gsl_pow_int((double)k/(double)n02,rd-2);
	  */
	/*};*/
    };   
  /* NOFST BEGIN
     for(rd1=2;rd1<=max_cov;rd1++)
     {
     for(rd2=2;rd2<=max_cov;rd2++)
     {
     int i,j,k;
     double x1,y1,x2,y2;
     for(k=1;k<=rd1+rd2-1;k++)
     { 
     for(i=0;i<=k;i++)
     {
	x1=(double)(k-i)/(double)n01;
	x2=(double)(i)/(double)n02;
	y1=1-x1;
	y2=1-x2;
	for(m1t=0;m1t<=m_bar;m1t++){
	for(m2t=0;m2t<=m_bar;m2t++){
	combfst.c_s[rd1-1][rd2-1]+=
	gsl_ran_hypergeometric_pdf(k-i,n01,n02,k)*
	((double)(m1t*(rd2-m2t)+m2t*(rd1-m1t))/(double)(rd1*rd2))*gsl_sf_choose(rd1,m1t)*gsl_sf_choose(rd2,m2t)*
	(gsl_pow_int(x2,m1t)*gsl_pow_int(y2,rd1-m1t)+gsl_pow_int(x2,rd1-m1t)*gsl_pow_int(y2,m1t))*
	(gsl_pow_int(x1,m2t)*gsl_pow_int(y1,rd2-m2t)+gsl_pow_int(x1,rd2-m2t)*gsl_pow_int(y1,m2t))
	//( (y2-x2)*x1*y1*(gsl_pow_int(y1,rd1-2)-gsl_pow_int(x1,rd1-2))+ (y1-x1)*x2*y2*(gsl_pow_int(y2,rd2-2)-gsl_pow_int(x2,rd2-2))
	//-(double)(rd1+rd2)*x1*y1*x2*y2*(gsl_pow_int(x1,rd1-2)+gsl_pow_int(y1,rd1-2))*(gsl_pow_int(x2,rd2-2)+gsl_pow_int(y2,rd2-2))+
	// 2*x1*y1*x2*y2*(gsl_pow_int(x1,rd1-2)-gsl_pow_int(y1,rd1-2))*(gsl_pow_int(x2,rd2-2)-gsl_pow_int(y2,rd2-2)) )
	/(double)k;
	};
	};
	};
      };
    };  
  };
  END */ 
  
    for(rd=1;rd<=max_cov*(n01-1); rd++){
	    vec_s[rd-1]=0;
	    vec_p[rd-1]=0;
	    vec_h[rd-1]=0;
    };
    
  for(rd=2*m_bar+1;rd<=max_cov; rd++){
    int k,covt;
	    vec0_s[rd-1]=0;
	    vec0_d[rd-1]=0;
	    vec0_h[rd-1]=0;
   for(k=1;k<n01;k++){
     vec_s[(k-1)+(n01-1)*(rd-1)]=(1-gsl_pow_int((double)k/(double)n01,rd)-	gsl_pow_int(1-(double)k/(double)n01,rd));
     vec_p[(k-1)+(n01-1)*(rd-1)]=(double)(2*k*(n01-k))/(double)(n01*n01);
     vec_h[(k-1)+(n01-1)*(rd-1)]=(double)(k*k)/(double)(n01*n01)+(double)(k)/(double)(n01*(rd-1))-(double)rd/(double)(rd-1)*gsl_pow_int((double)(k)/(double)(n01),rd);
     for(covt=1;covt<=m_bar;covt++){
	      vec_s[(k-1)+(n01-1)*(rd-1)]+=-( gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd)); 
	      vec_p[(k-1)+(n01-1)*(rd-1)]+=-( 2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))*(gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd))); 
	      vec_h[(k-1)+(n01-1)*(rd-1)]+=-( (double)(covt*covt)/(double)(rd*(rd-1))*gsl_ran_binomial_pdf(covt,	(double)k/(double)n01,rd)+(double)(rd-covt)*(double)(rd-covt)/(double)(rd*(rd-1))*gsl_ran_binomial_pdf(rd-covt,(double)k/(double)n01,rd)); 
     };       
     for(covt=m_bar+1;covt<=rd-m_bar-1;covt++){
       vec0_s[rd-1]+=gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
       vec0_d[rd-1]+=gsl_pow_int(2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/comb1.d_p[rd-1]-1/comb1.d_t[rd-1],2)*gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
       vec0_h[rd-1]+=gsl_pow_int((double)(covt*covt)/(double)(rd*(rd-1))/comb1.d_hq[rd-1]-2*(double)(covt*(rd-covt))/(double)(rd*(rd-1))/comb1.d_p[rd-1],2)*gsl_ran_binomial_pdf(covt,(double)k/(double)n01,rd)/k;
     };
   };
  };
      
  covmat=(double *)malloc((n01-1)*(n01-1)*sizeof(double));
  generate_covariance_matrix(covmat,n01);

  
  /* Initialize variables */
  pos=0;
  oldpos=0;
  pos_base1=0;
  pos_base2=0;
  n_window=1;
  cds_start=0;
  cds_end=0;
  
  //read_line_ms(bam_file1, window_size, n01, &n_pos1, int_positions1, array_counts1);
  
  fprintf(output_stat1, "window\tlength\tlength_outgroup\tread_depth\tS\tWatterson\tPi\tTajima_D\tvar_S\tvar_Watterson\tunnorm_FayWu_H\tFayWu_H\tdiv\tnonsyn_pol\tsyn_pol\tnonsyn_div\tsyn_div\talpha\n");
	  
  printf("Computing statistics for the window:");
  /* Run across all bases */
  for(pos=1;(ct1!=EOF); pos++)
    {
      DEB(printf("new line\n")); //debug 

      if (pos==(n_window-1)*window_size+1) 
	{	  
	  // INITIALIZE EVERYTHING HERE!!!!!!!!!!!!!!
	  
	  // printf("inizialize\n"); //debug 
	  
	  printf(" %lu ", n_window); 
	  if ( n_window % 10 == 0 ) printf("\n");
	  
	  test1.cov=0;
	  test1.l=0;
	  test1.l_out=0;
	  test1.s=0;
	  test1.num_t=0;
	  test1.num_p=0;
	  test1.num_hl=0;
	  test1.num_hq=0;
	  test1.den_t=0;
	  test1.den_p=0;
	  test1.den_hl=0;
	  test1.den_hq=0;
	  div=0;
	  for(rd=1;rd<=max_cov; rd++){
	    vec_rd[rd-1]=0;
//	    vec_s[rd-1]=0;
//	    vec_p[rd-1]=0;
//	    vec_h[rd-1]=0;
	  };
	  psyn=0;
	  dsyn=0;
	  pnon=0;
	  dnon=0;
	  
/*	  test2.cov=0;
	  test2.l=0;
	  test2.l_out=0;
	  test2.s=0;
	  test2.num_t=0;
	  test2.num_p=0;
	  test2.num_hl=0;
	  test2.num_hq=0;
	  test2.den_t=0;
	  test2.den_p=0;
	  test2.den_hl=0;
	  test2.den_hq=0;
	  
	  fst.l=0;
	  fst.gen_diff=0;
	  fst.c_s=0;	  */
	};
      /* Read data */
      //if(pos==TESTPOS){
        //           ct2="1";
          //   };
      if(ext_snps==1){
	while (pos_snp<pos) pos_snp=extract_pos_snpinput(list_snps);
      };
      
      if(if_gff==1){
	while (cds_end<pos) {
	  char *line_gff;
	  size_t n_line_gff;
	  n_line_gff=1;
	  line_gff=malloc(sizeof(char));
	  if(getline(&line_gff,&n_line_gff,gff)!=-1){
	  while ((line_gff[0]=='#')||(line_gff[0]=='\0')) {
	    getline(&line_gff,&n_line_gff,gff);
	  }; sscanf(line_gff,"%*s\t%*s\t%s\t%lu\t%lu\t%*s\t%c\t%u\t",feature,&cds_start,&cds_end,&strand,&phase_cds);
	  if((strcmp(feature,"CDS")!=0)&&(strcmp(feature,"cds")!=0)&&(strcmp(feature,"SO:0000316")!=0)) {
	    cds_end=0;
	  };
	  } else {cds_end=-1;};
	};
      };
      
          if (pos_base1<pos) 
	{
	  DEB(printf("reading new base from file 1\n")); //debug
	  read_line_pileup(bam_file1, min_qual, min_mqual, &pos_base1, &n_ref1, &n_alt_allele1, &rd1, n_alt_1, &ref_base1, &alt_base1);
	  //  read_line_ms(bam_file1, window_size, n01, &n_pos1, int_positions1, array_counts1);
	  //base_repetition(window_size, n01, &pos_base1, &n_pos1, int_positions1, &n_variation1, array_counts1, input_coverage1, error_rate1, r, &n_ref1, &n_alt_allele1, &n_tot_allele1, n_alt_1, &ref_base1, &alt_base1); rd1=n_tot_allele1;  
	};
/*      if (pos_base2<pos)
	{
	  DEB(printf("reading new base from file 2\n")); //debug
	  // All outputs set manually to 0 just for now.
	  //read_line_pileup(bam_file2, min_qual, min_mqual, &pos_base2, &n_ref2, &n_alt_allele2, &rd2, n_alt_2, &ref_base2, &alt_base2);
	  pos_base2 = 0;
	  n_ref2 = 0;
	  n_alt_allele2 = 0;
	  rd2 = 0;
	  int n_alt_2[5] = {0};
	  ref_base2 = 0;
	  alt_base2 = 0;
	};*/
      /* Extract statistics */
	/*if(pos>=fasta_length){
	  out_base=dmau_extract_outgroup_base(fasta_out,pos,oldpos,&fasta_length,ref_base1); //1;
	  } else {
	    out_base=ref_base1;
	  };*/
      if ((pos_base1==pos)&&(rd1>=min_cov)&&(rd1<=max_cov))
	{
	  if(outgroup_available==1){
	    out_base=extract_outgroup_base(fasta_out,pos,oldpos,fasta_length); //1;
	  } else {
	    out_base=0;
	  };
	  if((ext_snps==1)&&(pos_snp!=pos)) {
	    if (n_ref1>=n_alt_allele1){
	      n_ref1+=n_alt_allele1;
	      n_alt_allele1=0;
	    } else {
	      n_alt_allele1+=n_ref1;
	      n_ref1=0;
	    };
	  };
	  //printf("out base: %u\n",out_base); //debug
	  //printf("sample 1: "); //debug
	  if (rd1>=max(min_cov,2*m_bar+2))
	    {
	      if (if_gff==1){
	      if (cds_start<=pos){
		if(strand=='+'){
		  frame=((pos-cds_start+3-phase_cds)%3)+1;
		} else {
		  if(strand=='-'){
		    frame=((cds_end+3-phase_cds-pos)%3)+1;
		  } else frame=0;
		};
	      } else frame=0;
	      };
            extract_stats(&test1, &comb1, n01, n_ref1, n_alt_allele1, rd1, n_alt_1, ref_base1, alt_base1, out_base, m_bar);//+rd1*(pos_snp!=pos));
	      if(out_base!=0){
	      if (((n_alt_allele1>=rd1-m_bar)&&(alt_base1!=out_base))||((n_ref1>=rd1-m_bar)&&(ref_base1!=out_base))){
		div++;
		if (if_gff==1){
		  if (frame==3) dsyn++;
		  if ((frame==1)||(frame==2)) dnon++;
		};
	      };
	      if ((if_gff==1)&&(n_alt_allele1>m_bar)&&(n_ref1>m_bar)){
		if (frame==3) psyn++;
		if ((frame==1)||(frame==2)) pnon++;
	      };
	      };
	      vec_rd[rd1-1]++;
	    };
	  
	  
	  /*if ((pos_base2==pos)&&(rd2>=min_cov)&&(rd2<=max_cov))
	    {
	      //printf("sample 2: "); //debug
	      if (rd2>=max(min_cov,2*m_bar+2))
		{
		  extract_stats(&test2, &comb2, n02, n_ref2, n_alt_allele2, rd2, n_alt_2, ref_base2, alt_base2, out_base, m_bar);
		};
	      extract_fst(&fst, &combfst, n01, n_ref1, n_alt_allele1, rd1, n_alt_1, ref_base1, alt_base1, n02, n_ref2, n_alt_allele2, rd2, n_alt_2, ref_base2, alt_base2, out_base, m_bar);
	      //SNPS
#if PRINTSNPS == 1
	      extract_snps(pos, output_snps, n_alt_1, n_alt_2, m_bar);
#endif
	    };*/
	  oldpos=pos;
	/*}
      else 
	{
	  if ((pos_base2==pos)&&(rd2>=max(min_cov,2*m_bar+2))&&(rd2<=max_cov))
	    {
	      out_base=0;//extract_outgroup_base(fasta_out,pos,oldpos,fasta_length);
	      //printf("out base: %u\n",out_base); //debug
	      //printf("sample 2: "); //debug
	      extract_stats(&test2, &comb2, n02, n_ref2, n_alt_allele2, rd2, n_alt_2, ref_base2, alt_base2, out_base, m_bar);
	      oldpos=pos;
	    };*/
	};
      /* Print output */
      ct1=fgetc(bam_file1); ungetc(ct1,bam_file1);
      
      if ((pos==(n_window*window_size))||(ct1==EOF))
	{
	  double theta1_val, pi1_val, d1_val, h1_val, theta2_val, pi2_val, d2_val, h2_val, pia_val, fst_val, cov1_val, cov2_val, div_val, var_h, var_d, var_s, var0_s, var0_d, var0_h, vk_s[n01-1], vk_d[n01-1], vk_h[n01-1];
	  int k;
	  DEB(printf("printing output\n")); //debug
	  
	  var_s=0;
	  var_d=0;
	  var_h=0;
	  var0_s=0;
	  var0_d=0;
	  var0_h=0;
	  if ((test1.den_t>0)&&(test1.den_p>0)) {
	  if (test1.den_hq>0) {
	  for(k=1;k<n01;k++){
	    vk_s[k-1]=0;
	    vk_d[k-1]=0;
	    vk_h[k-1]=0;
	    for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
	      vk_s[k-1]+=vec_rd[rd-1]*vec_s[(k-1)+(n01-1)*(rd-1)];
	      vk_d[k-1]+=vec_rd[rd-1]*(vec_p[(k-1)+(n01-1)*(rd-1)]/test1.den_p-vec_s[(k-1)+(n01-1)*(rd-1)]/test1.den_t);
	      vk_h[k-1]+=vec_rd[rd-1]*(vec_h[(k-1)+(n01-1)*(rd-1)]/test1.den_hq-vec_p[(k-1)+(n01-1)*(rd-1)]/test1.den_p);
	    };
	  };
	  for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
	      var0_s+=vec_rd[rd-1]*vec0_s[rd-1];
	      var0_d+=vec_rd[rd-1]*vec0_d[rd-1];
	      var0_h+=vec_rd[rd-1]*vec0_h[rd-1];
	  };
	  for(k=1;k<n01;k++){
	    for(lt=1;lt<n01;lt++){
	      var_s+=vk_s[k-1]*vk_s[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
	      var_d+=vk_d[k-1]*vk_d[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
	      var_h+=vk_h[k-1]*vk_h[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
	    };
	  };
	  var0_d=var0_d/(double)(test1.l*test1.l);
	  var0_h=var0_h/(double)(test1.l_out*test1.l_out);
	  } else { 
	  for(k=1;k<n01;k++){
	    vk_s[k-1]=0;
	    vk_d[k-1]=0;
	    for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
	      vk_s[k-1]+=vec_rd[rd-1]*vec_s[(k-1)+(n01-1)*(rd-1)];
	      vk_d[k-1]+=vec_rd[rd-1]*(vec_p[(k-1)+(n01-1)*(rd-1)]/test1.den_p-vec_s[(k-1)+(n01-1)*(rd-1)]/test1.den_t);
	    };
	  };
	  for(rd=max(min_cov,2*m_bar+2);rd<=max_cov; rd++){
	      var0_s+=vec_rd[rd-1]*vec0_s[rd-1];
	      var0_d+=vec_rd[rd-1]*vec0_d[rd-1];
	  };
	  for(k=1;k<n01;k++){
	    for(lt=1;lt<n01;lt++){
	      var_s+=vk_s[k-1]*vk_s[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
	      var_d+=vk_d[k-1]*vk_d[lt-1]*covmat[(n01-1)*(lt-1)+k-1];
	    };
	  };
	  var0_d=var0_d/(double)(test1.l*test1.l);
	  };
	  };
	  
	  if(test1.l>0) { cov1_val=(double)(test1.cov)/(double)(test1.l); } else { cov1_val=-1; };
	  if(test1.den_t>0) { theta1_val=test1.num_t/test1.den_t; } else { theta1_val=-1; }; 
	  if(test1.den_p>0) { pi1_val=test1.num_p/test1.den_p; } else { pi1_val=-1; }; 
	  if((test1.den_t>0)&&(test1.den_p>0)) { d1_val=pi1_val-theta1_val; } else { d1_val=-1; }; 
	  if((test1.den_p>0)&&(test1.den_hq>0)) { h1_val=test1.num_hq/test1.den_hq-pi1_val; } else { h1_val=-1; }; 
	  if(test1.l_out>0) { div_val=(double)(div)/(double)(test1.l_out); } else { div_val=-1; };
	  /*if(test2.l>0) { cov2_val=(double)(test2.cov)/(double)(test2.l); } else { cov2_val=0; };
	  if(test2.den_t>0) { theta2_val=test2.num_t/test2.den_t; } else { theta2_val=0; }; 
	  if(test2.den_p>0) { pi2_val=test2.num_p/test2.den_p; } else { pi2_val=0; }; 
	  if((test2.den_t>0)&&(test2.den_p>0)) { d2_val=pi2_val-theta2_val; } else { d2_val=0; }; 
	  if((test2.den_hl>0)&&(test2.den_hq>0)) { h2_val=test2.num_hq/test2.den_hq-test2.num_hl/test2.den_hl; } else { h2_val=0; };
	  pia_val=(fst.gen_diff+fst.c_s*(pi1_val+pi2_val)/2)/(double)(fst.l);
	  fst_val=1-2*(pi1_val+pi2_val)/(pi1_val+pi2_val+2*pia_val);*/
	  var0_s=var0_s*theta1_val;
	  var0_d=var0_d*theta1_val;
	  var0_h=var0_h*theta1_val;
	  var_s=var_s*theta1_val*theta1_val;
	  var_d=var_d*theta1_val*theta1_val;
	  var_h=var_h*theta1_val*theta1_val;
	  
	  ////fprintf(output_stat1, "window %u\n",n_window); //debug COMPLETE!!!!!!!!!!!!!!!!
	  //  DEB(printf(output_stat1, "%u\t%u\t%u\t%f\t%u\t%f\t%f\t%f\t%f\n", n_window, test1.l, test1.l_out, cov1_val, test1.s, theta1_val, pi1_val, d1_val, h1_val)); 
	  //	  fprintf(output_stat1, "%lu\t%lu\t%lu\t%f\t%lu\t%f\t%f\t%f\t%f\t%f\n", n_window, test1.l, test1.l_out, cov1_val, test1.s, theta1_val, pi1_val, d1_val, h1_val,test1.den_t); 
	  fprintf(output_stat1, "%lu\t%lu\t%lu",n_window, test1.l, test1.l_out);
	  if (test1.l>0) {
	  fprintf(output_stat1, "\t%f\t%lu\t%f\t%f\t%f\t%f\t%f", cov1_val, test1.s, theta1_val, pi1_val, d1_val/sqrt(var0_d+var_d), var0_s+var_s, (var0_s+var_s)/(test1.den_t*test1.den_t)); 
	  } else {
	  fprintf(output_stat1, "\tNA\t0\tNA\tNA\tNA\tNA\tNA");  
	  };
	  if (test1.l_out>0) {
	  fprintf(output_stat1, "\t%f\t%f\t%f", h1_val, h1_val/sqrt(var0_h+var_h), div_val);	  
	  } else { 
	  fprintf(output_stat1, "\tNA\tNA\tNA");	    
	  };
	  if (if_gff==1) {
	  fprintf(output_stat1, "\t%lu\t%lu\t%lu\t%lu", pnon, psyn, dnon, dsyn);
	  if(dsyn*pnon==0){
	    if(psyn*dnon==0){
	      fprintf(output_stat1, "\tNA");	    
	    } else {
	      fprintf(output_stat1, "\tInf");	      
	    };
	  } else {	    
	    fprintf(output_stat1, "\t%f", 1-(double)(dsyn*pnon)/(double)(psyn*dnon));
	  };
	  } else { 
	  fprintf(output_stat1, "\tNA\tNA\tNA\tNA\tNA");	    
	  };
	  DEB(fprintf(output_stat1, "\tvars\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",var0_s,var_s,var0_s/(test1.den_t*test1.den_t+0.000001),var_s/(test1.den_t*test1.den_t+0.000001),var0_d/(theta1_val),var_d/(theta1_val*theta1_val),var0_h/(theta1_val),var_h/(theta1_val*theta1_val));)
	  fprintf(output_stat1, "\n");	    
	  
	  ////fprintf(output_stat2, "window %u\n",n_window); //debug COMPLETE!!!!!!!!!!!!!!!!
	  //  DEB(printf(output_stat2, "%u\t%u\t%u\t%f\t%u\t%f\t%f\t%f\t%f\n", n_window, test2.l, test2.l_out, cov2_val, test2.s, theta2_val, pi2_val, d2_val, h2_val));
	  //fprintf(output_stat2, "%u\t%u\t%u\t%f\t%u\t%f\t%f\t%f\t%f\n", n_window, test2.l, test2.l_out, cov2_val, test2.s, theta2_val, pi2_val, d2_val, h2_val);
	  ////fprintf(output_fst,  "window %u\n",n_window); //debug COMPLETE!!!!!!!!!!!!!!!!
	  //  DEB(printf(output_fst, "%u\t%u\t%f\t%f\t%f\t%f\n", n_window, fst.l, pi1_val, pi2_val, pia_val, fst_val));
	  //  fprintf(output_fst, "%u\t%u\t%f\t%f\t%f\t%f\n", n_window, fst.l, pi1_val, pi2_val, pia_val, fst_val);
	  n_window++; 
	};
      
      ////printf("%lu\t%lu\t%lu\t|\t%lu\t%lu\t%lu\t%lu\n",n_ref,n_alt_allele,rd,n_alt_A,n_alt_C,n_alt_G,n_alt_T); //debug

      //ct1=fgetc(bam_file1); ungetc(ct1,bam_file1);
      //ct2=fgetc(bam_file2); ungetc(ct2,bam_file2);
      
      ////ct=fgetc(bam_file); ungetc(ct,bam_file);
      ////printf("%c\n",ct); //debug
    };
  /* Close files, free gsl random number generator */  
  fclose(bam_file1);
  //fclose(bam_file2);
  if (outgroup_available==1) { fclose(fasta_out); };
  if (ext_snps==1) { fclose(list_snps); };
  if (if_gff==1) { fclose(gff); };
  fclose(output_stat1);
  //fclose(output_stat2);
  //fclose(output_fst);
  gsl_rng_free(r);
  //SNPS
#if PRINTSNPS == 1
  fclose(output_snps);
#endif 
  printf("\n");
  printf("Computation of statistics completed.\n");
}
