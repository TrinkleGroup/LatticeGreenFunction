#ifndef _INVERT_EXPRESSION_HH_
#define _INVERT_EXPRESSION_HH_

#include <string>
#include "Expression.hh"

class InvertExpression: public Expression
{
private:
	Expression *exp;

public:
	InvertExpression();
	InvertExpression(Expression *exp);

	virtual std::string getValue() const;
};

#endif
