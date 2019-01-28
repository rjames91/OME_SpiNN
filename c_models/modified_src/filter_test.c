// DO NOT EDIT! THIS FILE WAS GENERATED FROM src/filter_test.c

#include <math.h>
#include <stdio.h>
int main (void)
{
        double filter_1;
	double sub;
	double result;
	double stapes_a;
	double stapes_b;
	stapes_a = -0.99937168146928201384326939660240896046161651611328125;
	stapes_b = 1.0 + stapes_a;
	
	filter_1 = stapes_b * -6.8853628760871392730007883632840612560987092125031061584650160511955618858337402343750000000000000000e-16;
	sub = stapes_a * -5.2230809704382550137617236885060008499323186480049752145049524187925271689891815185546875000000000000e-17;

	result = filter_1 - sub;

	printf("result = %.100e\n",result);
	printf("result = %.100e\n",(float)result);


        return 0;
}

