#include "Common.h"

// Takes the line and return each token
extern void tokenizer(const char* delimeters, const char* line, std::vector<std::string>& tokens)
{
	char token[255] = "\0";
	int tokenIndex = 0;
	int index = 0;
	bool isToken = false;

	// Remove delimeter above of the line
	for (int deli = 0; deli < strlen(delimeters); deli++)
	{
		while (line[index] == delimeters[deli])
		{
			index++;
		}
	}


	// Check the each character of the line
	for (index; index < strlen(line); index++)
	{
		for (int deli = 0; deli < strlen(delimeters); deli++)
		{
			if (line[index] != delimeters[deli])	// Add characters to the token until find the delimeter
			{
				token[tokenIndex++] = line[index];

				isToken = true;
			}
			else if (isToken)	// Token ends, add into the std::vector
			{
				token[tokenIndex + 1] = '\0';
				tokens.push_back(std::string(token));

				tokenIndex = 0;
				memset(token, 0, 255);

				isToken = false;
			}
		}	// End of for strlen(delimeters)

	}	// End of for strlen(line)

	if (isToken)
	{
		token[tokenIndex + 1] = '\0';
		tokens.push_back(token);
	}
}