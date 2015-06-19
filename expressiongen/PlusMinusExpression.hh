#ifndef _PLUS_MINUS_EXPRESSION_HH_
#define _PLUS_MINUS_EXPRESSION_HH_

#include <string>
#include <vector>
#include "Expression.hh"

class PlusMinusExpression: public Expression
{
private:
	std::vector <Expression *> v;
	char sign;

public:
	PlusMinusExpression();
	PlusMinusExpression(char sign, const std::vector <Expression *> &v);

	virtual std::string getValue() const;
};

#endif
