#include <iostream>
#ifdef assert
    #undef assert
#endif
#include <TError.h> // R__ASSERT
#define assert R__ASSERT // Fixing assert for interpreted ROOT macros 

void Test() 
{
	std::cout << "Hello world!" << std::endl; 
	assert(false); 
}
