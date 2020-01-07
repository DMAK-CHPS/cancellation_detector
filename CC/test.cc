#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <random>

/**************************************************
*				Macro et typedef				  *
**************************************************/

#define bit(A,i) ((A >> i) & 1)
#define min(a,b) ((a < b) ? a : b)
#define max(a,b) ((a > b) ? a : b)
#define bit_a_b(X,a,b) (((X << (sizeof(X)*8 - max(a,b) - 1)) >> (sizeof(X)*8 - max(a,b) - 1)) & ((X >> min(a,b) ) << min(a,b)))

typedef unsigned long long int u_64;
typedef unsigned int u_32;
std::default_random_engine generator(time(NULL));
std::uniform_real_distribution<double> distribution(-0.5,0.5);

/**************************************************
*					Affichages					  *
**************************************************/

static void print_double(double d)
{
	u_64 a = *((u_64*) &d);

	for(int i = 63; i >= 0; i--)
		(i == 63 || i == 52) ? printf("%lld ", bit(a,i)) : printf("%lld", bit(a,i));

	printf("\n");
}

static void print_float(float f)
{
	u_32 a = *((u_32*) &f);

	for(int i = 31; i >= 0; i--)
		(i == 31 || i == 23) ? printf("%d ", bit(a,i)) : printf("%d", bit(a,i));

	printf("\n");
}

void cancell_generator_f(float *f1, float *f2, size_t size)
{
	std::uniform_int_distribution<int> b(0,1);

	(*f1) = (*f2) = 1;

	for(unsigned int i = 1; i <= 23; i++){
		int bit = b(generator);
		float f = 1/pow(2,i);
		
		if(i <= size){
			(*f1) += bit * f;
			(*f2) += bit * f;
		}else if(i == size+1){
			(*f1) += 1 * f;
		}else{
			(*f1) += bit * f;
		}
	}
}

void cancell_generator_d(double *d1, double *d2, size_t size)
{
	std::uniform_int_distribution<int> b(0,1);

	(*d1) = (*d2) = 1;

	for(unsigned int i = 1; i <= 52; i++){
		int bit = b(generator);
		float f = 1/pow(2,i);
		
		if(i <= size){
			(*d1) += bit * f;
			(*d2) += bit * f;
		}else if(i == size+1){
			(*d1) += 1 * f;
		}else{
			(*d1) += bit * f;
		}
	}
}

/**************************************************
*				Fonctions projet				  *
**************************************************/

double _mca_rand()
{
  return distribution(generator);	
}

static int _mca_inexactq(double* qa, int size) {
  int ex = 0;
  frexp((*qa), &ex);
  (*qa) = (*qa) + pow(2,ex-size) * _mca_rand();
  return 1;
}

static int _mca_inexactd(float* da, int size) {
  int ex = 0;
  frexp((*da), &ex);
  (*da) = (*da) + pow(2,ex-size) * _mca_rand();
  return 1;
}

int cancell_float(float f1, float f2)
{
	int ea, eb, er;
	frexp(f1, &ea);
	frexp(f2, &eb);
	frexp(f1-f2, &er);

	return max(ea,eb) - er;	
}

int cancell_double(double d1, double d2)
{
	int ea, eb, er;
	frexp(d1, &ea);
	frexp(d2, &eb);
	frexp(d1-d2, &er);

	return max(ea,eb) - er;
}

/**************************************************
*				Tris de tableaux	  	  		  *
**************************************************/

void sortf(float *a, u_64 n)
{
  	for (u_64 i = 0; i < n; i++)
    	for (u_64 j = i + 1; j < n; j++)
      		if (a[j] < a[i])
			{
			  float tmp = a[i];

			  a[i] = a[j];
			  a[j] = tmp;
			}
}

void sortd(double *a, u_64 n)
{
  	for (u_64 i = 0; i < n; i++)
    	for (u_64 j = i + 1; j < n; j++)
      		if (a[j] < a[i])
			{
			  double tmp = a[i];

			  a[i] = a[j];
			  a[j] = tmp;
			}
}

/**************************************************
*				Tests fonctionnels	  	  		  *
**************************************************/

void test_f(void (*generator) (float* , float*, size_t), int (*detect) (float, float), int (*inexact) (float*, int))
{
	float f1 = 0, f2 = 0;

	for(unsigned int i = 0; i <= 23; i++){
		printf("\n\n------- %d --------\n\n", i+1);
		for(unsigned int j = 0; j < 10; j++){
			generator(&f1, &f2, i);
			float t[1000] = {}, true_val = f1-f2;

			for(unsigned int x = 0; x < 1000; x++){
				t[x] = true_val;
				inexact(&t[x],detect(f1,f2));
			}
			sortf(t,1000);
			float avg = (t[499]+t[500])/2.0;

			printf("Exact: %.2a | Approximation: %.2a\n", true_val, avg);
		}
	}
}


void test_d(void (*generator) (double* , double*, size_t), int (*detect) (double, double), int (*inexact) (double*, int))
{
	double d1 = 0, d2 = 0;

	for(unsigned int i = 0; i <= 52; i++){
		printf("\n\n------- %d --------\n\n", i+1);
		for(unsigned int j = 0; j < 10; j++){
			generator(&d1, &d2, i);
			double t[1000] = {}, true_val = d1 - d2;
			for(unsigned int x = 0; x < 1000; x++){
				t[x] = true_val;
				inexact(&t[x],detect(d1,d2));
			}
			sortd(t,1000);
			double avg = (t[499]+t[500])/2.0;

			printf("Exact: %.2a | Approximation: %.2a\n", true_val, avg);
		}
	}
}

void test_inexact_f(int (*inexact) (float*, int))
{
	FILE *f = fopen("dat/inexact_f.dat","wr");
	
	for(unsigned int i = 0; i <= 52; i++){
		for(unsigned int x = 0; x < 1000; x++){
			float f1 = 0;
			inexact(&f1,i);
			fprintf(f,"%d %la\n", i, f1 * pow(2,i)/i);
		}
	}

	fclose(f);

	system("gnuplot dat/plot_inexact_f.sh");
}

void test_inexact_d(int (*inexact) (double*, int))
{
	FILE *f = fopen("dat/inexact_d.dat","wr");
	
	for(unsigned int i = 0; i <= 52; i++){
		for(unsigned int x = 0; x < 1000; x++){
			double d1 = 0;
			inexact(&d1,i);
			fprintf(f,"%d %la\n", i, d1 * pow(2,i)/i);
		}
	}

	fclose(f);

	system("gnuplot dat/plot_inexact_d.sh");
}

void test_detect_f(void (*generator) (float* , float*, size_t), int (*detect) (float, float))
{
	FILE *f = fopen("dat/detect_f.dat","wr");
	
	for(unsigned int i = 0; i <= 23; i++){
		for(unsigned int x = 0; x < 1000; x++){
			float f1 = 0, f2 = 0;
			generator(&f1,&f2,i);
			fprintf(f,"%d %d\n", i, detect(f1,f2));
		}
	}

	fclose(f);

	system("gnuplot dat/plot_detect_f.sh");
}

void test_detect_d(void (*generator) (double* , double*, size_t), int (*detect) (double, double))
{
	FILE *f = fopen("dat/detect_d.dat","wr");
	
	for(unsigned int i = 0; i <= 52; i++){
		for(unsigned int x = 0; x < 1000; x++){
			double d1 = 0, d2 = 0;
			generator(&d1,&d2,i);
			fprintf(f,"%d %d\n", i, detect(d1,d2));
		}
	}

	fclose(f);

	system("gnuplot dat/plot_detect_d.sh");
}

/**************************************************
*				Tests fonctionnels	  	  		  *
**************************************************/

void plot_d(void (*generator) (double* , double*, size_t), int (*detect) (double, double), int (*inexact) (double*, int))
{
	double d1 = 0, d2 = 0;
	FILE *f = fopen("dat/double.dat","wr");

	for(unsigned int i = 0; i < 51; i++){
		for(unsigned int j = 0; j < 10; j++){
			generator(&d1, &d2, i);
			double t[1000] = {}, true_val = d1 - d2;
			for(unsigned int x = 0; x < 1000; x++){
				t[x] = true_val;
				inexact(&t[x],detect(d1,d2));
			}

			for(unsigned int x = 0; x < 1000; x++){
				fprintf(f,"%lf\n",  t[x]);
			}
		}
	}

	fclose(f);

	system("gnuplot dat/plot_double.sh");
}

void plot_f(void (*generator) (float* , float*, size_t), int (*detect) (float, float), int (*inexact) (float*, int))
{
	float f1 = 0, f2 = 0;
	FILE *f = fopen("dat/float.dat","wr");

	for(unsigned int i = 0; i < 22; i++){
		for(unsigned int j = 0; j < 10; j++){
			generator(&f1, &f2, i);
			float t[1000] = {}, true_val = f1-f2;

			for(unsigned int x = 0; x < 1000; x++){
				t[x] = true_val;
				inexact(&t[x],detect(f1,f2));
			}

			for(unsigned int x = 0; x < 1000; x++){
				fprintf(f,"%lf\n",  t[x]);
			}
		}
	}

	fclose(f);

	system("gnuplot dat/plot_float.sh");
}

int main(int argc, char const *argv[])
{
	test_inexact_f(_mca_inexactd);

	test_inexact_d(_mca_inexactq);

	test_detect_f(cancell_generator_f, cancell_float);

	test_detect_d(cancell_generator_d, cancell_double);
/*
	test_f(cancell_generator_f, cancell_float, _mca_inexactd);

	test_d(cancell_generator_d, cancell_double, _mca_inexactq);
*/
	plot_f(cancell_generator_f, cancell_float, _mca_inexactd);

	plot_d(cancell_generator_d, cancell_double, _mca_inexactq);

	return 0;
}