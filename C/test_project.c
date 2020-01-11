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
#define expd(X) ((X << 1) >> 24)
#define expq(X) ((X << 1) >> 53)

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

	if(expq(a) > expq(b))
	{
		b = bit_a_b(b,63,63) + bit_a_b(a,62,52) + ((bit_a_b(b,51,0) + (unsigned long long int)pow(2,52)) >> (bit_a_b(a,62,52) - bit_a_b(b,62,52)));
	}
	else if(expq(b) < expq(a))
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

	if(expd(a) > expd(b))
	{
		b = bit_a_b(b,31,31) + bit_a_b(a,30,23) + ((bit_a_b(b,22,0) + (unsigned int)pow(2,23)) >> (bit_a_b(a,30,23) - bit_a_b(b,30,23)));
	}
	else if(expd(b) < expd(a))
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
	float res = f1-f2;
	unsigned int a = *((unsigned int*) &f1);
	unsigned int b = *((unsigned int*) &f2);
	unsigned int r = *((unsigned int*) &(res));

	int ea = expd(a)-127, eb = expd(b)-127, er  = expd(r)-127;

	return max(ea,eb)-er;	
}

int cancell_double_v4(double d1, double d2)
{
	double res = d1-d2;
	unsigned long long int a = *((unsigned long long int*) &d1);
	unsigned long long int b = *((unsigned long long int*) &d2);
	unsigned long long int r = *((unsigned long long int*) &(res));

	int ea = expq(a)-1023, eb = expq(b)-1023, er = expq(r)-1023;

	return max(ea,eb)-er;
}


int cancell_float_v5(float f1, float f2)
{
	unsigned int a = *((unsigned int*) &f1);
	unsigned int b = *((unsigned int*) &f2);

	int exp_a = expd(a)-127, exp_b = expd(b)-127;

	if(exp_a > exp_b){
		a = 8388608 + bit_a_b(a,0,22);
		b = (8388608 + bit_a_b(b,0,22)) >> (exp_a - exp_b);
	}else if(exp_b > exp_a){
		a = (8388608 + bit_a_b(a,0,22)) >> (exp_b - exp_a);
		b = 8388608 + bit_a_b(b,0,22);
	}

	int i = 23;
	unsigned int res = a-b;

	while(i >= 0 && bit(res,i) == 0)
		i--;

	return (i == -1) ? 0 : 23 - i;
} 

int cancell_double_v5(double d1, double d2)
{
	unsigned long long int a = *((unsigned long long int*) &d1);
	unsigned long long int b = *((unsigned long long int*) &d2);

	int exp_a = expq(a)-1023, exp_b = expq(b)-1023;

	if(exp_a > exp_b){
		a = 4503599627370496 + bit_a_b(a,0,51);
		b = (4503599627370496 + bit_a_b(b,0,51)) >> (exp_a - exp_b);
	}else if(exp_b > exp_a){
		a = (4503599627370496 + bit_a_b(a,0,51)) >> (exp_b - exp_a);
		b = 4503599627370496 + bit_a_b(b,0,51);
	}

	int i = 52;
	unsigned long long int res = a-b;

	while(i >= 0 && bit(res,i)  == 0)
		i--;

	return (i == -1) ? 0 : 52 - i;
} 

int cancell_float_v6(float f1, float f2)
{
	float max = max(f1,f2);
	float min = min(f1,f2);
	return log2(max) - log2(max-min);
}

int cancell_double_v6(double d1, double d2)
{
	double max = max(d1,d2);
	double min = min(d1,d2);
	return log2(max) - log2(max-min);
}

int cancell_float_v7(float f1, float f2)
{
	unsigned int a = *((unsigned int*) &f1);
	unsigned int b = *((unsigned int*) &f2);
	float res_f = f1-f2;
	unsigned int fres = *((unsigned int*) &res_f);

	int exp_a = expd(a)-127, exp_b = expd(b)-127;

	if(exp_a > exp_b){
		a = 8388608 + bit_a_b(a,0,22);
		b = (8388608 + bit_a_b(b,0,22)) >> (exp_a - exp_b);
	}else if(exp_b > exp_a){
		a = (8388608 + bit_a_b(a,0,22)) >> (exp_b - exp_a);
		b = 8388608 + bit_a_b(b,0,22);
	}

	unsigned int res = a-b;
	int rr = 0, rf = 0;

	if (res){
	  res = (res ^ (res - 1)) >> 1;
	  for (rr = 0; res && rr < 23; rr++){
	    res >>= 1;}
	}

	if (fres){
	  fres = (fres ^ (fres - 1)) >> 1;
	  for (rf = 0; fres && rf < 23; rf++){
	    fres >>= 1;}
	}

	return rf - rr;
}

int cancell_double_v7(double d1, double d2)
{
	unsigned long long int a = *((unsigned long long int*) &d1);
	unsigned long long int b = *((unsigned long long int*) &d2);
	double res_d = d1-d2;
	unsigned long long int dres = *((unsigned long long int*) &res_d);

	int exp_a = expq(a)-1023, exp_b = expq(b)-1023;

	if(exp_a > exp_b){
		a = 4503599627370496 + bit_a_b(a,0,51);
		b = (4503599627370496 + bit_a_b(b,0,51)) >> (exp_a - exp_b);
	}else if(exp_b > exp_a){
		a = (4503599627370496 + bit_a_b(a,0,51)) >> (exp_b - exp_a);
		b = 4503599627370496 + bit_a_b(b,0,51);
	}

	unsigned long long int res = a-b;
	int rr = 0, rd = 0;

	if (res){
	  res = (res ^ (res - 1)) >> 1;
	  for (rr = 0; res && rr < 52; rr++){
	    res >>= 1;}
	}

	if (dres){
	  dres = (dres ^ (dres - 1)) >> 1;
	  for (rd = 0; dres && rd < 52; rd++){
	    dres >>= 1;}
	}

	return rd - rr;
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

static float f_tab[] = {
  /* 0 */
  0x1.804f9ep+3,    // 0 10000010 10000000010011111001111
  0x1.0df39cp+0,    // 0 01111111 00001101111100111001110
  0x1.009f3cp+2,    // 0 10000010 10000000010011111001111
  /* 1 */
  0x1.804f9ep+3,    // 0 10000010 10000000010011111001111
  0x1p+3,           // 0 10000010 00000000000000000000000
  0x1.009f3cp+2,    // 0 10000001 00000000100111110011110
  /* 2 */
  0x1.7000a8p+0,    // 0 01111111 01110000000000001010100
  0x1p+0,           // 0 01111111 00000000000000000000000
  0x1.c002ap-2,     // 0 01111101 11000000000000101010000
  /* 3 */
  0x1.bf3402p+1,    // 0 10000000 10111111001101000000001
  0x1.8p+1,         // 0 10000000 10000000000000000000000
  0x1.f9a01p-2,     // 0 01111101 11111001101000000001000
  /* 4 */
  0x1.955b32p+3,    // 0 10000010 10010101010110110011001
  0x1.8p+3,   // 0 10000010 10000000000000000000000
  0x1.55b32p-1,   // 0 01111110 01010101101100110010000
  /* 5 */
  0x1.feb942p+3,    // 0 10000010 11111110101110010100001
  0x1.fp+3,   // 0 10000010 11110000000000000000000
  0x1.d7284p-2,   // 0 01111101 11010111001010000100000
  /* 6 */
  0x1.0d4018p+0,    // 0 01111111 00001101010000000001100
  0x1.08p+0,    // 0 01111111 00001000000000000000000
  0x1.5006p-6,    // 0 01111001 01010000000001100000000
  /* 7 */
  0x1.3be11p+0,   // 0 01111111 00111011111000010001000
  0x1.38p+0,    // 0 01111111 00111000000000000000000
  0x1.f088p-7,    // 0 01111000 11110000100010000000000
  /* 8 */
  0x1.af5eb6p+4,    // 0 10000011 10101111010111101011011
  0x1.aep+4,    // 0 10000011 10101110000000000000000
  0x1.5eb6p-4,    // 0 01111011 01011110101101100000000
  /* 9 */
  0x1.62bb1ep+2,    // 0 10000001 01100010101110110001111
  0x1.62p+2,    // 0 10000001 01100010000000000000000
  0x1.763cp-7,    // 0 01111000 01110110001111000000000
  /* 10 */
  0x1.be67e6p+3,    // 0 10000010 10111110011001111110011
  0x1.bep+3,    // 0 10000010 10111110000000000000000
  0x1.9f98p-7,    // 0 01111000 10011111100110000000000
  /* 11 */
  0x1.b5b782p+1,    // 0 10000000 10110101101101111000001
  0x1.b58p+1,   // 0 10000000 10110101100000000000000
  0x1.bc1p-10,    // 0 01110101 10111100000100000000000
  /* 12 */
  0x1.7c9d3cp+1,    // 0 10000000 01111100100111010011110
  0x1.7c8p+1,   // 0 10000000 01111100100000000000000
  0x1.d3cp-11,    // 0 01110100 11010011110000000000000
  /* 13 */
  0x1.75de86p+4,    // 0 10000011 01110101110111101000011
  0x1.75dp+4,   // 0 10000011 01110101110100000000000
  0x1.d0cp-9,   // 0 01110110 11010000110000000000000
  /* 14 */
  0x1.7976aap+4,    // 0 10000011 01111001011101101010101
  0x1.797p+4,   // 0 10000011 01111001011100000000000
  0x1.aa8p-10,    // 0 01110101 10101010100000000000000
  /* 15 */
  0x1.c612b6p+1,    // 0 10000000 11000110000100101011011
  0x1.c61p+1,   // 0 10000000 11000110000100000000000
  0x1.5bp-14,   // 0 01110001 01011011000000000000000
  /* 16 */
  0x1.9d239p+4,   // 0 10000011 10011101001000111001000
  0x1.9d22p+4,    // 0 10000011 10011101001000100000000
  0x1.9p-12,    // 0 01110011 10010000000000000000000
  /* 17 */
  0x1.f47e9ep+4,    // 0 10000011 11110100011111101001111
  0x1.f47ep+4,    // 0 10000011 11110100011111100000000
  0x1.3cp-13,   // 0 01110010 00111100000000000000000
  /* 18 */
  0x1.bf375cp+3,    // 0 10000010 10111111001101110101110
  0x1.bf37p+3,    // 0 10000010 10111111001101110000000
  0x1.7p-15,    // 0 01110000 01110000000000000000000
  /* 19 */
  0x1.5389ecp+4,    // 0 10000011 01010011100010011110110
  0x1.5389cp+4,   // 0 10000011 01010011100010011100000
  0x1.6p-15,    // 0 01110000 01100000000000000000000
  /* 20 */
  0x1.777d7cp+2,    // 0 10000001 01110111011111010111110
  0x1.777d6p+2,   // 0 10000001 01110111011111010110000
  0x1.cp-18,    // 0 01101101 11000000000000000000000
  /* 21 */
  0x1.b3c9a8p+2,    // 0 10000001 10110011110010011010100
  0x1.b3c9ap+2,   // 0 10000001 10110011110010011010000
  0x1p-19,    // 0 01101100 00000000000000000000000
  /* 22 */
  0x1.6c9ceep+1,    // 0 10000000 01101100100111001110111
  0x1.6c9ce8p+1,    // 0 10000000 01101100100111001110100
  0x1.8p-21,    // 0 01101010 10000000000000000000000
  /* 23 */
  0x1.0df39ep+0,    // 0 01111111 00001101111100111001111
  0x1.0df39cp+0,    // 0 01111111 00001101111100111001110
  0x1p-23           // 0 01101000 00000000000000000000000
};


static double q_tab[] = {
  /* 0 */
  0x1.e9aad7f278b3p+0,    // 0 01111111111 1110100110101010110101111111001001111000101100110000
  0x1.2761330b2acp+4,	  // 0 10000000011 0010011101100001001100110000101100101010110000000000
  -0x1.08c6858c0334dp+4,  // 1 10000000011 0000100011000110100001011000110000000011001101001101

  /* 1 */
  0x1.e9aad7f278b3p+0,    // 0 01111111111 1110100110101010110101111111001001111000101100110000
  0x1p+0,                 // 0 01111111111 0000000000000000000000000000000000000000000000000000
  0x1.d355afe4f166p-1,    // 0 01111111110 1101001101010101101011111110010011110001011001100000
  /* 2 */
  0x1.74734915e38adp+0,   // 0 01111111111 0111010001110011010010010001010111100011100010101101
  0x1p+0,   // 0 01111111111 0000000000000000000000000000000000000000000000000000
  0x1.d1cd24578e2b4p-2,   // 0 01111111101 1101000111001101001001000101011110001110001010110100
  /* 3 */
  0x1.b26429af2c1e8p+4,   // 0 10000000011 1011001001100100001010011010111100101100000111101000
  0x1.8p+4,   // 0 10000000011 1000000000000000000000000000000000000000000000000000
  0x1.93214d7960f4p+1,    // 0 10000000000 1001001100100001010011010111100101100000111101000000
  /* 4 */
  0x1.da7c951993c34p+1,   // 0 10000000000 1101101001111100100101010001100110010011110000110100
  0x1.cp+1,   // 0 10000000000 1100000000000000000000000000000000000000000000000000
  0x1.a7c951993c34p-3,    // 0 01111111100 1010011111001001010100011001100100111100001101000000
  /* 5 */
  0x1.bbfbe3af1d92dp+0,   // 0 01111111111 1011101111111011111000111010111100011101100100101101
  0x1.bp+0,   // 0 01111111111 1011000000000000000000000000000000000000000000000000
  0x1.7f7c75e3b25ap-5,    // 0 01111111010 0111111101111100011101011110001110110010010110100000
  /* 6 */
  0x1.3f3a10ff4d516p+1,   // 0 10000000000 0011111100111010000100001111111101001101010100010110
  0x1.38p+1,    // 0 10000000000 0011100000000000000000000000000000000000000000000000
  0x1.ce843fd35458p-5,    // 0 01111111010 1100111010000100001111111101001101010100010110000000
  /* 7 */
  0x1.3ad1323f74281p+0,   // 0 01111111111 0011101011010001001100100011111101110100001010000001
  0x1.38p+0,    // 0 01111111111 0011100000000000000000000000000000000000000000000000
  0x1.68991fba1408p-7,    // 0 01111111000 0110100010011001000111111011101000010100000010000000
  /* 8 */
  0x1.7f2d263ccbda2p+0,   // 0 01111111111 0111111100101101001001100011110011001011110110100010
  0x1.7ep+0,    // 0 01111111111 0111111000000000000000000000000000000000000000000000
  0x1.2d263ccbda2p-8,   // 0 01111110111 0010110100100110001111001100101111011010001000000000
  /* 9 */
  0x1.37f25202b2fbep+2,   // 0 10000000001 0011011111110010010100100000001010110010111110111110
  0x1.37p+2,    // 0 10000000001 0011011100000000000000000000000000000000000000000000
  0x1.e4a40565f7cp-7,   // 0 01111111000 1110010010100100000001010110010111110111110000000000
  /* 10 */
  0x1.d85f8091b819bp+2,   // 0 10000000001 1101100001011111100000001001000110111000000110011011
  0x1.d8p+2,    // 0 10000000001 1101100000000000000000000000000000000000000000000000
  0x1.7e0246e066cp-8,   // 0 01111110111 0111111000000010010001101110000001100110110000000000
  /* 11 */
  0x1.41f39e547c011p+1,   // 0 10000000000 0100000111110011100111100101010001111100000000010001
  0x1.41cp+1,   // 0 10000000000 0100000111000000000000000000000000000000000000000000
  0x1.9cf2a3e0088p-10,    // 0 01111110101 1001110011110010101000111110000000001000100000000000
  /* 12 */
  0x1.1a9759db495d2p+4,   // 0 10000000011 0001101010010111010110011101101101001001010111010010
  0x1.1a8p+4,   // 0 10000000011 0001101010000000000000000000000000000000000000000000
  0x1.759db495d2p-8,    // 0 01111110111 0111010110011101101101001001010111010010000000000000
  /* 13 */
  0x1.286d1ff1d39dep+3,   // 0 10000000010 0010100001101101000111111111000111010011100111011110
  0x1.286p+3,   // 0 10000000010 0010100001100000000000000000000000000000000000000000
  0x1.a3fe3a73bcp-10,   // 0 01111110101 1010001111111110001110100111001110111100000000000000
  /* 14 */
  0x1.56d7644507949p+3,   // 0 10000000010 0101011011010111011001000100010100000111100101001001
  0x1.56dp+3,   // 0 10000000010 0101011011010000000000000000000000000000000000000000
  0x1.d91141e524p-11,   // 0 01111110100 1101100100010001010000011110010100100100000000000000
  /* 15 */
  0x1.d42e140ab7f2p+4,    // 0 10000000011 1101010000101110000101000000101010110111111100100000
  0x1.d42cp+4,    // 0 10000000011 1101010000101100000000000000000000000000000000000000
  0x1.0a055bf9p-11,   // 0 01111110100 0000101000000101010110111111100100000000000000000000
  /* 16 */
  0x1.a3d7f187c4752p+0,   // 0 01111111111 1010001111010111111100011000011111000100011101010010
  0x1.a3d6p+0,    // 0 01111111111 1010001111010110000000000000000000000000000000000000
  0x1.f187c4752p-16,    // 0 01111101111 1111000110000111110001000111010100100000000000000000
  /* 17 */
  0x1.3fc7f499a631bp+2,   // 0 10000000001 0011111111000111111101001001100110100110001100011011
  0x1.3fc7p+2,    // 0 10000000001 0011111111000111000000000000000000000000000000000000
  0x1.e9334c636p-15,    // 0 01111110000 1110100100110011010011000110001101100000000000000000
  /* 18 */
  0x1.c96c57a36e13cp+3,   // 0 10000000010 1100100101101100010101111010001101101110000100111100
  0x1.c96cp+3,    // 0 10000000010 1100100101101100000000000000000000000000000000000000
  0x1.5e8db84fp-15,   // 0 01111110000 0101111010001101101110000100111100000000000000000000
  /* 19 */
  0x1.fb4c23d19975bp+1,   // 0 10000000000 1111101101001100001000111101000110011001011101011011
  0x1.fb4cp+1,    // 0 10000000000 1111101101001100000000000000000000000000000000000000
  0x1.1e8ccbad8p-18,    // 0 01111101101 0001111010001100110010111010110110000000000000000000
  /* 20 */
  0x1.cd8f99365e341p+0,   // 0 01111111111 1100110110001111100110010011011001011110001101000001
  0x1.cd8f8p+0,   // 0 01111111111 1100110110001111100000000000000000000000000000000000
  0x1.9365e341p-20,   // 0 01111101011 1001001101100101111000110100000100000000000000000000
  /* 21 */
  0x1.8cf87eac7e8c8p+4,   // 0 10000000011 1000110011111000011111101010110001111110100011001000
  0x1.8cf87p+4,   // 0 10000000011 1000110011111000011100000000000000000000000000000000
  0x1.d58fd19p-17,    // 0 01111101110 1101010110001111110100011001000000000000000000000000
  /* 22 */
  0x1.12280ef45eb8cp+4,   // 0 10000000011 0001001000101000000011101111010001011110101110001100
  0x1.122808p+4,    // 0 10000000011 0001001000101000000010000000000000000000000000000000
  0x1.bd17ae3p-18,    // 0 01111101101 1011110100010111101011100011000000000000000000000000
  /* 23 */
  0x1.0bedb2d848ec2p+2,   // 0 10000000001 0000101111101101101100101101100001001000111011000010
  0x1.0bedbp+2,   // 0 10000000001 0000101111101101101100000000000000000000000000000000
  0x1.6c24761p-21,    // 0 01111101010 0110110000100100011101100001000000000000000000000000
  /* 24 */
  0x1.b2dc13ffb1e38p+0,   // 0 01111111111 1011001011011100000100111111111110110001111000111000
  0x1.b2dc12p+0,    // 0 01111111111 1011001011011100000100100000000000000000000000000000
  0x1.ffb1e38p-24,    // 0 01111100111 1111111110110001111000111000000000000000000000000000
  /* 25 */
  0x1.8ecfa6a5e18bap+4,   // 0 10000000011 1000111011001111101001101010010111100001100010111010
  0x1.8ecfa6p+4,    // 0 10000000011 1000111011001111101001100000000000000000000000000000
  0x1.4bc3174p-21,    // 0 01111101010 0100101111000011000101110100000000000000000000000000
  /* 26 */
  0x1.3dd2decc4ca44p+2,   // 0 10000000001 0011110111010010110111101100110001001100101001000100
  0x1.3dd2de8p+2,   // 0 10000000001 0011110111010010110111101000000000000000000000000000
  0x1.313291p-24,   // 0 01111100111 0011000100110010100100010000000000000000000000000000
  /* 27 */
  0x1.f37f94e2de34cp+4,   // 0 10000000011 1111001101111111100101001110001011011110001101001100
  0x1.f37f94cp+4,   // 0 10000000011 1111001101111111100101001100000000000000000000000000
  0x1.16f1a6p-23,   // 0 01111101000 0001011011110001101001100000000000000000000000000000
  /* 28 */
  0x1.c03d40d56011cp+3,   // 0 10000000010 1100000000111101010000001101010101100000000100011100
  0x1.c03d40cp+3,   // 0 10000000010 1100000000111101010000001100000000000000000000000000
  0x1.56011cp-25,   // 0 01111100110 0101011000000001000111000000000000000000000000000000
  /* 29 */
  0x1.dba6fada0d974p+3,   // 0 10000000010 1101101110100110111110101101101000001101100101110100
  0x1.dba6fadp+3,   // 0 10000000010 1101101110100110111110101101000000000000000000000000
  0x1.41b2e8p-26,   // 0 01111100101 0100000110110010111010000000000000000000000000000000
  /* 30 */
  0x1.c5878dd5ce50ep+0,   // 0 01111111111 1100010110000111100011011101010111001110010100001110
  0x1.c5878ddp+0,   // 0 01111111111 1100010110000111100011011101000000000000000000000000
  0x1.739438p-30,   // 0 01111100001 0111001110010100001110000000000000000000000000000000
  /* 31 */
  0x1.50d8321b996d3p+3,   // 0 10000000010 0101000011011000001100100001101110011001011011010011
  0x1.50d83218p+3,    // 0 10000000010 0101000011011000001100100001100000000000000000000000
  0x1.ccb698p-28,   // 0 01111100011 1100110010110110100110000000000000000000000000000000
  /* 32 */
  0x1.276748bdb4826p+0,   // 0 01111111111 0010011101100111010010001011110110110100100000100110
  0x1.276748bcp+0,    // 0 01111111111 0010011101100111010010001011110000000000000000000000
  0x1.b4826p-32,    // 0 01111011111 1011010010000010011000000000000000000000000000000000
  /* 33 */
  0x1.56c34027d7b55p+0,   // 0 01111111111 0101011011000011010000000010011111010111101101010101
  0x1.56c34027p+0,    // 0 01111111111 0101011011000011010000000010011100000000000000000000
  0x1.af6aap-33,    // 0 01111011110 1010111101101010101000000000000000000000000000000000
  /* 34 */
  0x1.0272b29d7ea94p+2,   // 0 10000000001 0000001001110010101100101001110101111110101010010100
  0x1.0272b29dp+2,    // 0 10000000001 0000001001110010101100101001110100000000000000000000
  0x1.faa5p-32,   // 0 01111011111 1111101010100101000000000000000000000000000000000000
  /* 35 */
  0x1.c2284c0c7622bp+1,   // 0 10000000000 1100001000101000010011000000110001110110001000101011
  0x1.c2284c0c4p+1,   // 0 10000000000 1100001000101000010011000000110001000000000000000000
  0x1.b1158p-34,    // 0 01111011101 1011000100010101100000000000000000000000000000000000
  /* 36 */
  0x1.75b85dae5c7cbp+3,   // 0 10000000010 0111010110111000010111011010111001011100011111001011
  0x1.75b85dae4p+3,   // 0 10000000010 0111010110111000010111011010111001000000000000000000
  0x1.c7cbp-33,   // 0 01111011110 1100011111001011000000000000000000000000000000000000
  /* 37 */
  0x1.ce12e7cade9fep+4,   // 0 10000000011 1100111000010010111001111100101011011110100111111110
  0x1.ce12e7cadp+4,   // 0 10000000011 1100111000010010111001111100101011010000000000000000
  0x1.d3fcp-33,   // 0 01111011110 1101001111111100000000000000000000000000000000000000
  /* 38 */
  0x1.030e23a64fb54p+0,   // 0 01111111111 0000001100001110001000111010011001001111101101010100
  0x1.030e23a648p+0,    // 0 01111111111 0000001100001110001000111010011001001000000000000000
  0x1.ed5p-38,    // 0 01111011001 1110110101010000000000000000000000000000000000000000
  /* 39 */
  0x1.bd377c26768bcp+3,   // 0 10000000010 1011110100110111011111000010011001110110100010111100
  0x1.bd377c2674p+3,    // 0 10000000010 1011110100110111011111000010011001110100000000000000
  0x1.45ep-36,    // 0 01111011011 0100010111100000000000000000000000000000000000000000
  /* 40 */
  0x1.181a2a9d5bec6p+2,   // 0 10000000001 0001100000011010001010101001110101011011111011000110
  0x1.181a2a9d5ap+2,    // 0 10000000001 0001100000011010001010101001110101011010000000000000
  0x1.ec6p-38,    // 0 01111011001 1110110001100000000000000000000000000000000000000000
  /* 41 */
  0x1.a5aa5f0b96df3p+0,   // 0 01111111111 1010010110101010010111110000101110010110110111110011
  0x1.a5aa5f0b96p+0,    // 0 01111111111 1010010110101010010111110000101110010110000000000000
  0x1.be6p-41,    // 0 01111010110 1011111001100000000000000000000000000000000000000000
  /* 42 */
  0x1.c26019263bc11p+3,   // 0 10000000010 1100001001100000000110010010011000111011110000010001
  0x1.c26019263b8p+3,   // 0 10000000010 1100001001100000000110010010011000111011100000000000
  0x1.044p-39,    // 0 01111011000 0000010001000000000000000000000000000000000000000000
  /* 43 */
  0x1.d4261910bdf5bp+0,   		// 0 01111111111 1101010000100110000110010001000010111101111101011011
  0x1.d4261910bdcp+0,   		// 0 01111111111 1101010000100110000110010001000010111101110000000000
  0x1.ad8p-43,    				// 0 01111010100 1010110110000000000000000000000000000000000000000000
  /* 44 */
  0x1.9805cc0b84f51p+0,			// 0 01111111111 1001100000000101110011000000101110000100111101010001
  0x1.9805cc0b84ep+0,			// 0 01111111111 1001100000000101110011000000101110000100111000000000
  0x1.51p-44,					// 0 01111010011 0101000100000000000000000000000000000000000000000000
  /* 45 */
  0x1.d23d12ef3f8acp+2,   		// 0 10000000001 1101001000111101000100101110111100111111100010101100
  0x1.d23d12ef3f8p+2,   		// 0 10000000001 1101001000111101000100101110111100111111100000000000
  0x1.58p-43,   				// 0 01111010100 0101100000000000000000000000000000000000000000000000
  /* 46 */
  0x1.303627ab3fb79p+3,   // 0 10000000010 0011000000110110001001111010101100111111101101111001
  0x1.303627ab3fbp+3,   // 0 10000000010 0011000000110110001001111010101100111111101100000000
  0x1.e4p-43,   // 0 01111010100 1110010000000000000000000000000000000000000000000000
  /* 47 */
  0x1.2e54ef49c3ba1p+1,   // 0 10000000000 0010111001010100111011110100100111000011101110100001
  0x1.2e54ef49c3b8p+1,    // 0 10000000000 0010111001010100111011110100100111000011101110000000
  0x1.08p-46,   // 0 01111010001 0000100000000000000000000000000000000000000000000000
  /* 48 */
  0x1.1b758da877839p+3,   // 0 10000000010 0001101101110101100011011010100001110111100000111001
  0x1.1b758da87782p+3,    // 0 10000000010 0001101101110101100011011010100001110111100000100000
  0x1.9p-45,    // 0 01111010010 1001000000000000000000000000000000000000000000000000
  /* 49 */
  0x1.03409e8c0ef38p+2,   // 0 10000000001 0000001101000000100111101000110000001110111100111000
  0x1.03409e8c0ef3p+2,    // 0 10000000001 0000001101000000100111101000110000001110111100110000
  0x1p-47,    // 0 01111010000 0000000000000000000000000000000000000000000000000000
  /* 50 */
  0x1.a96f176c098bep+4,   // 0 10000000011 1010100101101111000101110110110000001001100010111110
  0x1.a96f176c098b8p+4,   // 0 10000000011 1010100101101111000101110110110000001001100010111000
  0x1.8p-46,    		  // 0 01111010001 1000000000000000000000000000000000000000000000000000
  /* 51 */
  0x1.3c58c4e2a9552p+3,   // 0 10000000010 0011110001011000110001001110001010101001010101010010
  0x1.3c58c4e2a955p+3,    // 0 10000000010 0011110001011000110001001110001010101001010101010000
  0x1p-48,    			  // 0 01111001111 0000000000000000000000000000000000000000000000000000
  /* 52 */
  0x1.74e6f391548f5p+2,   // 0 10000000001 0111010011100110111100111001000101010100100011110101
  0x1.74e6f391548f4p+2,   // 0 10000000001 0111010011100110111100111001000101010100100011110100
  0x1p-50,    			  // 0 01111001101 0000000000000000000000000000000000000000000000000000
};

void test_detect_f(const char* title, int (*f) (float f1, float f2))
{
	unsigned int t[24] = {}, cpt = 0, sum = 0;

    for(unsigned int j = 0; j <= 23; j++){
	    float a = f_tab[j*3], b = f_tab[j*3+1], res = f_tab[j*3+2];
	    if(f(a,b) != j){
	    	t[cpt] = abs(f(a,b)-j);
	    	sum += t[cpt];
	    	cpt++;
			printf("%d %d\n", j, f(a,b));
	    }
    }

	if(cpt == 0){
		printf("%s\n", title);
		printf("SUCCES !!\n\n");
	}else{
		sort_u(t,cpt);
		printf("%s\n", title);
		printf("- nombre d'erreur = %d\n- pourcentage d'erreur = %f\n- plus petite erreur = %d\n- plus grosse erreur = %d\n- erreur mediane = %f\n- erreur moyenne = %f\n\n", cpt, (float)cpt/23, t[0], t[cpt-1], (t[cpt/2]+t[cpt/2+1])/2.0, (float)sum/cpt);
	}
}

void test_detect_q(const char* title, int (*f) (double d1, double d2))
{
	unsigned int t[53] = {}, cpt = 0, sum = 0;

    for(unsigned int j = 0; j <= 52; j++){
    	double a = q_tab[j*3], b = q_tab[j*3+1], res = q_tab[j*3+2];
    	if(f(a,b) != j){
	    	t[cpt] = abs(f(a,b)-j);
	    	sum += t[cpt];
	    	cpt++;
			printf("%d %d\n", j, f(a,b));
    	}
    }

	if(cpt == 0){
		printf("%s\n", title);
		printf("SUCCES !!\n\n");
	}else{
		sort_u(t,cpt);
		printf("%s\n", title);
		printf("- nombre d'erreur = %d\n- pourcentage d'erreur = %f\n- plus petite erreur = %d\n- plus grosse erreur = %d\n- erreur mediane = %f\n- erreur moyenne = %f\n\n", cpt, (float)cpt/52, t[0], t[cpt-1], (t[cpt/2]+t[cpt/2+1])/2.0, (float)sum/cpt);
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
/*
	double d1 =  0x1.e9aad7f278b3p+0;    // 0 01111111111 1110100110101010110101111111001001111000101100110000
  	double d2 = 0x1.2761330b2acp+4;	  // 0 10000000011 0010011101100001001100110000101100101010110000000000
  	printf("\t%la,\t// ", d1-d2); print_double(d1 - d2);
*/
	test_detect_f("Version original float", cancell_float_v1);
	test_detect_f("Version 2 float", cancell_float_v2);
	test_detect_f("Version 3 float", cancell_float_v3);
	test_detect_f("Version verrou float", cancell_float_v4);
	test_detect_f("Version bit float", cancell_float_v5);
	test_detect_f("Version 6", cancell_float_v6);
	test_detect_f("Version 7", cancell_float_v7);

	test_detect_q("Version original double", cancell_double_v1);
	test_detect_q("Version 2 double", cancell_double_v2);
	test_detect_q("Version 3 double", cancell_double_v3);
	test_detect_q("Version verrou double", cancell_double_v4);
	test_detect_q("Version bit double", cancell_double_v5);
	test_detect_q("Version 6", cancell_double_v6);
	test_detect_q("Version 7", cancell_double_v7);

	test_perf_f("Version original float", cancell_float_v1, 1000);
	test_perf_f("Version 2 float", cancell_float_v2, 1000);
	test_perf_f("Version 3 float", cancell_float_v3, 1000);
	test_perf_f("Version verrou float", cancell_float_v4, 1000);
	test_perf_f("Version bit float", cancell_float_v5, 1000);
	test_perf_f("Version 6", cancell_float_v6, 1000);
	test_perf_f("Version 7", cancell_float_v6, 1000);

	test_perf_q("Version original double", cancell_double_v1, 1000);
	test_perf_q("Version 2 double", cancell_double_v2, 1000);
	test_perf_q("Version 3 double", cancell_double_v3, 1000);
	test_perf_q("Version verrou double", cancell_double_v4, 1000);
	test_perf_q("Version bit double", cancell_double_v5, 1000);
	test_perf_q("Version 7", cancell_double_v6, 1000);

	return 0;
}
