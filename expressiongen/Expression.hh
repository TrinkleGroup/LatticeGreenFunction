#ifndef _EXPRESSION_HH_
#define _EXPRESSION_HH_

#include <string>

class Expression
{
private:
	std::string expr;

public:
	Expression();
	Expression(std::string &s);
	Expression(const char *s);
	virtual const std::string &getExpr() const;
	virtual void setExpr(const std::string &s);

	virtual std::string getValue() const = 0;
};

#endif
