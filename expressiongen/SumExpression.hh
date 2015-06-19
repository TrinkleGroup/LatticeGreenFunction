#ifndef _SUM_EXPRESSION_HH_
#define _SUM_EXPRESSION_HH_

#include <string>
#include <vector>
#include "Expression.hh"

class SumExpression: public Expression
{
private:
	std::vector <Expression *> v;

public:
	SumExpression();
	SumExpression(const std::vector <Expression *> &v);

	virtual std::string getValue() const;
};

#endif
