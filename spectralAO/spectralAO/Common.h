// For tokenizer()
#include <vector>

// The operation we are performing
enum RENDERING_TYPE
{
	USUAL = 0x00000001,
	MH_BASES = 0x00000010,
	SPECTRAL_AO = 0x00000100
};

// EXTERN FUNCTIONS

// Takes the line and return each token
extern void tokenizer(const char* delimeters, const char* line, std::vector<std::string>& tokens);