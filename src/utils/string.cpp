#include "string.hpp"

void strToLowerCase(std::string &s) {
	std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){return std::tolower(c);});
}
