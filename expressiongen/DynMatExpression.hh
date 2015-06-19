#ifndef _DYN_MAT_EXPRESSION_HH_
#define _DYN_MAT_EXPRESSION_HH_

#include <string>
#include "Expression.hh"

class DynMatExpression: public Expression
{
private:
	int n;

public:
	DynMatExpression();
	DynMatExpression(int n);

	virtual std::string getValue() const;
};

#endif
