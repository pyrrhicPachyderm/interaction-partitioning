#ifndef UTILS_MATH_HPP
#define UTILS_MATH_HPP

#include <bit>

//Integer square root.
//Implementing Heron's method.
//Adapting https://en.wikipedia.org/wiki/Integer_square_root#Example_implementation_in_C
template<typename T> T isqrt(T n) {
	//isqrt(0) = 0, isqrt(1) = 1;
	if(n <= 1) return n;
	
	//Initial guess.
	//We want 2^(floor(log2(n)/2) + 1).
	//bit_width gives floor(log2(n))+1.
	T x0 = 1 << ((std::bit_width(n) - 1)/2 + 1);
	
	//Update and loop until equal.
	T x1 = (x0 + n / x0) / 2;
	while(x1 < x0) {
		x0 = x1;
		x1 = (x0 + n / x0) / 2;
	}
	return x0;
}

#endif
