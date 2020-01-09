#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <time.h>
#include <string.h>

/**************************************************
*				Macro et typedef				  *
**************************************************/

#define bit(A,i) ((A >> i) & 1)
#define min(a,b) ((a < b) ? a : b)
#define max(a,b) ((a > b) ? a : b)
#define bit_a_b(X,a,b) (((X << (sizeof(X)*8 - max(a,b) - 1)) >> (sizeof(X)*8 - max(a,b) - 1)) & ((X >> min(a,b) ) << min(a,b)))

typedef unsigned long long int u_64;
typedef unsigned int u_32;

static inline u_64 rdtsc (void)
{
  //64 bit variables
  unsigned long long a, d;

  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  
  return (d << 32) | a;
}

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

static void print_u64(u_64 u)
{
	for(int i = 63; i >= 0; i--)
		(i == 63 || i == 52) ? printf("%lld ", bit(u,i)) : printf("%lld", bit(u,i));

	printf("\n");
}

static void print_u32(u_32 u)
{
	for(int i = 31; i >= 0; i--)
		(i == 31 || i == 23) ? printf("%d ", bit(u,i)) : printf("%d", bit(u,i));

	printf("\n");
}

/**************************************************
*				Detecteurs de cancellation		  *
**************************************************/

int cancell_double_v1(double d1, double d2)
{
	unsigned long long int a = *((unsigned long long int*) &d1);
	unsigned long long int b = *((unsigned long long int*) &d2);

	if(bit_a_b(a,62,52) > bit_a_b(b,62,52))
	{
		b = bit_a_b(b,63,63) + bit_a_b(a,62,52) + ((bit_a_b(b,51,0) + (unsigned long long int)pow(2,52)) >> (bit_a_b(a,62,52) - bit_a_b(b,62,52)));
	}
	else if(bit_a_b(a,62,52) < bit_a_b(b,62,52))
	{
		a = bit_a_b(a,63,63) + bit_a_b(b,62,52) + ((bit_a_b(a,51,0) + (unsigned long long int)pow(2,52)) >> (bit_a_b(b,62,52) - bit_a_b(a,62,52)));
	}

	int i;
	for(i = 51; i >= 0 && bit(a,i) == bit(b,i); i--){}

	return (i == -1) ? 0: 53-(i+1);	
}

int cancell_float_v1(float f1, float f2)
{
	unsigned int a = *((unsigned int*) &f1);
	unsigned int b = *((unsigned int*) &f2);

	if(bit_a_b(a,30,23) > bit_a_b(b,30,23))
	{
		b = bit_a_b(b,31,31) + bit_a_b(a,30,23) + ((bit_a_b(b,22,0) + (unsigned int)pow(2,23)) >> (bit_a_b(a,30,23) - bit_a_b(b,30,23)));
	}
	else if(bit_a_b(a,30,23) < bit_a_b(b,30,23))
	{
		a = bit_a_b(a,31,31) + bit_a_b(b,30,23) + ((bit_a_b(a,22,0) + (unsigned int)pow(2,23)) >> (bit_a_b(b,30,23) - bit_a_b(a,30,23)));
	}

	int i;
	for(i = 22; i >= 0 && bit(a,i) == bit(b,i); i--){}

	return (i == -1) ? 0 : 24-(i+1);
}

int cancell_double_v2(double d1, double d2)
{
	d1 = d1 - d2;

	unsigned long long int a = *((unsigned long long int*) &d1);	

	int i;
	for(i = 0; i <= 51 && bit(a,i) == 0; i++){}

	return (i == 52) ? 0: i;
}

int cancell_float_v2(float f1, float f2)
{
	f1 = f1 - f2;

	unsigned int a = *((unsigned int*) &f1);

	int i;
	for(i = 0; i <= 22 && bit(a,i) == 0; i++){}

	return (i == 23) ? 0: i;
}

int cancell_double_v3(double d1, double d2)
{	
	if(round(log2(d1)) > round(log2(d2)))
		d2 *= pow(2, log2(d1) - log2(d2));
	else if(round(log2(d2)) > round(log2(d1)))
		d1 *= pow(2, log2(d2) - log2(d1));
	if(d1 == 0 || d2 == 0)
		return 0;

	return (d1 > d2) ? -round(log2(1-d2/d1)): -round(log2(1-d1/d2));
}

int cancell_float_v3(float f1, float f2)
{
	if(round(log2(f1)) > round(log2(f2)))
		f2 *= pow(log2(f1) - log2(f2), 2);
	else if(round(log2(f2)) > round(log2(f1)))
		f1 *= pow(log2(f2) - log2(f1), 2);
	else if(f1 == 0 || f2 == 0)
		return 0;

	return (f1 > f2) ? -round(log2(1-f2/f1)): -round(log2(1-f1/f2));
}

int cancell_float_v4(float f1, float f2)
{
	int ea, eb, er;
	frexp(f1, &ea);
	frexp(f2, &eb);
	frexp(f1-f2, &er);

	return max(ea,eb) - er;	
}

int cancell_double_v4(double d1, double d2)
{
	int ea, eb, er;
	frexp(d1, &ea);
	frexp(d2, &eb);
	frexp(d1-d2, &er);

	return max(ea,eb) - er;
}

/**************************************************
*			Generateur de cancellation		  	  *
**************************************************/

int rand_a_b(int a, int b){
    return rand()%(b-a) + a;
}

u_64 lrand_a_b(u_64 a, u_64 b){
	return rand()%(b-a) + a;
}

double frand_a_b(double a, double b){
    return (rand()/(double)RAND_MAX)*(b-a)+a;
}

void cancell_generator_f(float *f1, float *f2, u_32 size)
{
	(*f1) = (*f2) = 1;

	for(unsigned int i = 1; i <= 23; i++){
		int bit = rand_a_b(0,2);
		float f = 1/pow(2,i);
		
		if(i <= size){
			(*f1) += bit * f;
			(*f2) += bit * f;
		}else{
			(*f1) += bit * f;
		}
	}

	int e = rand_a_b(0,5);

	(*f1) *= pow(2,e);
	(*f2) *= pow(2,e);
}

void cancell_generator_d(double *d1, double *d2, u_64 size)
{
	(*d1) = (*d2) = 1;

	for(unsigned int i = 1; i <= 52; i++){
		int bit = rand_a_b(0,2);
		float f = 1/pow(2,i);
		
		if(i <= size){
			(*d1) += bit * f;
			(*d2) += bit * f;
		}else{
			(*d1) += bit * f;
		}
	}

	int e = rand_a_b(0,5);

	(*d1) *= pow(2,e);
	(*d2) *= pow(2,e);
}

/**************************************************
*				Tris de tableaux	  	  		  *
**************************************************/

void sort_u(unsigned int *restrict a, u_64 n)
{
  	for (u_64 i = 0; i < n; i++)
    	for (u_64 j = i + 1; j < n; j++)
      		if (a[j] < a[i])
			{
			  unsigned int tmp = a[i];

			  a[i] = a[j];
			  a[j] = tmp;
			}
}

void sort_ll(u_64 *restrict a, u_64 n)
{
  	for (u_64 i = 0; i < n; i++)
    	for (u_64 j = i + 1; j < n; j++)
      		if (a[j] < a[i])
			{
			  unsigned int tmp = a[i];

			  a[i] = a[j];
			  a[j] = tmp;
			}
}

/**************************************************
*				Test de detections	  	  		  *
**************************************************/

void test_detect_f(const char* title, int (*f) (float f1, float f2), unsigned int nb_iteration)
{
	float f1 = 0.0, f2 = 0.0;
	u_32 size = 0;
	unsigned int t[nb_iteration], cpt = 0, sum = 0;

	for(unsigned int i = 0; i < nb_iteration; i++){
		size = rand_a_b(0,23);
		cancell_generator_f(&f1, &f2, size);
		t[cpt] = f(f1,f2);
		if(t[cpt] <= size){
			/*
			print_float(f1);
			print_float(f2);
			print_float(f1-f2);
			printf("%d %ld\n",t[cpt], size);
			*/
			t[cpt] = abs(t[cpt] - size);
			sum += t[cpt];
			cpt++; 
		}
	}

	if(cpt == 0){
		printf("%s\n", title);
		printf("SUCCES !!\n\n");
	}else{
		sort_u(t,cpt);
		printf("%s\n", title);
		printf("- nombre d'erreur = %d\n- pourcentage d'erreur = %f\n- plus petite erreur = %d\n- plus grosse erreur = %d\n- erreur mediane = %f\n- erreur moyenne = %f\n\n", cpt, (float)cpt/nb_iteration, t[0], t[cpt-1], (t[cpt/2]+t[cpt/2+1])/2.0, (float)sum/cpt);
	}
}

void test_detect_q(const char* title, int (*f) (double d1, double d2), unsigned int nb_iteration)
{
	double d1 = 0.0, d2 = 0.0;
	u_64 size = 0;
	unsigned int t[nb_iteration], cpt = 0, sum = 0;

	for(unsigned int i = 0; i < nb_iteration; i++){
		size = rand_a_b(0,52);
		cancell_generator_d(&d1, &d2, size);
		t[cpt] = f(d1,d2);
		if(t[cpt] <= size){
/*			print_float(f1);
			print_float(f2);
			print_float(f1-f2);
			printf("\n");
*/
			t[cpt] = abs(t[cpt] - size);
			sum += t[cpt];
			cpt++; 
		}
	}

	if(cpt == 0){
		printf("%s\n", title);
		printf("SUCCES !!\n\n");
	}else{
		sort_u(t,cpt);
		printf("%s\n", title);
		printf("- nombre d'erreur = %d\n- pourcentage d'erreur = %f\n- plus petite erreur = %d\n- plus grosse erreur = %d\n- erreur mediane = %f\n- erreur moyenne = %f\n\n", cpt, (float)cpt/nb_iteration, t[0], t[cpt-1], (t[cpt/2]+t[cpt/2+1])/2.0, (float)sum/cpt);
	}
}

/**************************************************
*				Test performances	  	  		  *
**************************************************/

void test_perf_f(const char* title, int (*f) (float f1, float f2), unsigned int nb_iteration)
{
	float f1 = 0.0, f2 = 0.0;
	u_32 size = 0;
	u_64 t[nb_iteration], sum = 0;

	for(unsigned int i = 0; i < nb_iteration; i++){
		f1 = frand_a_b(-8,+8);
		f2 = frand_a_b(-8,+8);
		u_64 start = 0, stop = 0;
		start = rdtsc();
		for(unsigned int j = 0; j < 10; j++){
			size = f(f1,f2);
		}
		stop = rdtsc();
		t[i] = (stop-start)/10.0;
		sum += t[i];
	}

	sort_ll(t,nb_iteration);
	printf("%s\n", title);
	printf("- nombre de cycles moyen = %lf\n- nombre de cycles mediane = %lf\n\n", (float)sum/nb_iteration, (t[nb_iteration/2]+t[nb_iteration/2+1])/2.0);
}

void test_perf_q(const char* title, int (*f) (double d1, double d2), unsigned int nb_iteration)
{
	double d1 = 0.0, d2 = 0.0;
	u_64 size = 0;
	u_64 t[nb_iteration], sum = 0;

	for(unsigned int i = 0; i < nb_iteration; i++){
		d1 = frand_a_b(-8,+8);
		d2 = frand_a_b(-8,+8);
		u_64 start = 0, stop = 0;
		start = rdtsc();
		for(unsigned int j = 0; j < 10; j++){
			size = f(d1,d2);
		}
		stop = rdtsc();
		t[i] = (stop-start)/10.0;
		sum += t[i];
	}

	sort_ll(t,nb_iteration);
	printf("%s\n", title);
	printf("- nombre de cycles moyen = %lf\n- nombre de cycles mediane = %lf\n\n", (float)sum/nb_iteration, (t[nb_iteration/2]+t[nb_iteration/2+1])/2.0);
}


int main(int argc, char const *argv[])
{
	srand(time(NULL));

	test_detect_f("Version original float", cancell_float_v1, 1000);
	test_detect_f("Version 2 float", cancell_float_v2, 1000);
	test_detect_f("Version 3 float", cancell_float_v3, 1000);
	test_detect_f("Version verrou float", cancell_float_v4, 1000);

	test_detect_q("Version original double", cancell_double_v1, 1000);
	test_detect_q("Version 2 double", cancell_double_v2, 1000);
	test_detect_q("Version 3 double", cancell_double_v3, 1000);
	test_detect_q("Version verrou double", cancell_double_v4, 1000);

	test_perf_f("Version original float", cancell_float_v1, 1000);
	test_perf_f("Version 2 float", cancell_float_v2, 1000);
	test_perf_f("Version 3 float", cancell_float_v3, 1000);
	test_perf_f("Version verrou float", cancell_float_v4, 1000);

	test_perf_q("Version original double", cancell_double_v1, 1000);
	test_perf_q("Version 2 double", cancell_double_v2, 1000);
	test_perf_q("Version 3 double", cancell_double_v3, 1000);
	test_perf_q("Version verrou double", cancell_double_v4, 1000);

	return 0;
}
