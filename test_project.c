#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define bit(A,i) ((A >> i) & 1)

typedef unsigned long long int u_64;
typedef unsigned int u_32;
typedef struct{
	char sign;
	u_32 exp;
	u_64 mant;
}double_struct;

typedef struct{
	char sign;
	char exp;
	u_32 mant;
}float_struct;

u_64 double_sign = pow(2,63);
u_64 double_exp = (pow(2,63)-1)-(pow(2,52)-1);
u_64 double_mant = (pow(2,52)-1);

u_32 float_sign = pow(2,31);
u_32 float_exp = (pow(2,31)-1)-(pow(2,23)-1);
u_32 float_mant = (pow(2,23)-1);


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

static void double_to_struct(double* d, double_struct* d_str)
{
	u_64 a = *((u_64*) d);

	d_str->sign = (a & double_sign) >> 63;
	d_str->exp = (a & double_exp) >> 52;
	d_str->mant = (a & double_mant);
}

static void float_to_struct(float* f, float_struct* f_str)
{
	u_32 a = *((u_32*) f);

	f_str->sign = (a & float_sign) >> 31;
	f_str->exp = (a & float_exp) >> 23;
	f_str->mant = (a & float_mant);
}

static void struct_to_double(double_struct* d_str, double* d)
{
	u_64 a = d_str->sign;
	a = (a << 11) + d_str->exp;
	a = (a << 52) + d_str->mant;
}

static void struct_to_float(float_struct* f_str, float* f)
{
	u_32 a = f_str->sign;
	a = (a << 8) + f_str->exp;
	a = (a << 23) + f_str->mant;
}

static u_64 struct_to_u64(double_struct* d_str)
{
	u_64 a = 0;
	a = d_str->sign;
	a = (a << 11) + d_str->exp;
	a = (a << 52) + d_str->mant;

	return a;
}

static u_32 struct_to_u32(float_struct* f_str)
{
	u_32 a = 0;
	a = f_str->sign;
	a = (a << 8) + f_str->exp;
	a = (a << 23) + f_str->mant;

	return a;
}



static int cancell_double(double d1, double d2)
{
	double_struct s_a, s_b;
	double_to_struct(&d1, &s_a);
	double_to_struct(&d2, &s_b);

	if(s_a.exp > s_b.exp)
	{
		s_b.mant = s_b.mant + pow(2,52);
		s_b.mant = s_b.mant >> (s_a.exp - s_b.exp);
		s_b.exp = s_a.exp;
	}
	else if(s_a.exp < s_b.exp)
	{
		s_a.mant = s_a.mant + pow(2,52);
		s_a.mant = s_a.mant >> (s_b.exp - s_a.exp);
		s_a.exp = s_b.exp;
	}

	u_64 a = struct_to_u64(&s_a) ^ struct_to_u64(&s_b);

	int i;
	for(i = 51; i >= 0 && bit(a,i) == 0; i--){}

	return 52-(i+1);	
}


static int cancell_float(float f1, float f2)
{
	float_struct s_a, s_b;
	float_to_struct(&f1, &s_a);
	float_to_struct(&f2, &s_b);

	if(s_a.exp > s_b.exp)
	{
		s_b.mant = s_b.mant + pow(2,23);
		s_b.mant = s_b.mant >> (s_a.exp - s_b.exp);
		s_b.exp = s_a.exp;
	}
	else if(s_a.exp < s_b.exp)
	{
		s_a.mant = s_a.mant + pow(2,23);
		s_a.mant = s_a.mant >> (s_b.exp - s_a.exp);
		s_a.exp = s_b.exp;
	}

	u_32 a = struct_to_u32(&s_a) ^ struct_to_u32(&s_b);

	int i;
	for(i = 22; i >= 0 && bit(a,i) == 0; i--){}

	return 23-(i+1);
}

int main(int argc, char const *argv[])
{
	float d1 = 6.53564;
	float d2 = 6.53565;

	printf("a: ");
	print_float(d1);
	printf("b: ");
	print_float(d2);
	printf("nombre de bits cancellation: %d\n", cancell_float(d1,d2));
	printf("a-b: ");
	print_float(d1 - d2);

	return 0;
}

/*
static int cancell(double d1, double d2)
{
	unsigned long long int a, b, sign_a, sign_b, mantisse_a, mantisse_b, exp_a, exp_b;
	memcpy(&a,&d1,8);
	memcpy(&b,&d2,8);

	sign_a = (a & double_sign) >> 63;
	sign_b = (b & double_sign) >> 63;
	exp_a = (a & double_exp) >> 52;
	exp_b = (b & double_exp) >> 52;
	mantisse_a = (a & double_mant);
	mantisse_b = (b & double_mant);

	if(exp_a > exp_b)
	{
		mantisse_b = mantisse_b + pow(2,52);
		mantisse_b = mantisse_b >> (exp_a - exp_b);
		b = sign_b;
		b = (b << 11) + exp_a;
		b = (b << 52) + mantisse_b;		
	}
	else if(exp_a < exp_b)
	{
		mantisse_a = mantisse_a + pow(2,52);
		mantisse_a = mantisse_a >> (exp_b - exp_a);
		a = sign_a;
		a = (a << 11) & exp_b;
		a = (a << 52) & mantisse_a;
	}

	printf("a: ");
	print_u64(a);
	printf("b: ");
	print_u64(b);

	a = (a ^ b);

	int i;
	for(i = 51; i >= 0 && bit(a,i) == 0; i--){}

	return 52-(i+1);	
}
*/