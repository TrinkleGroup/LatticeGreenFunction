#include <string>
#include <vector>
#include <sstream>
#include "Expression.hh"
#include "SumExpression.hh"

SumExpression::SumExpression() : v(0) {
	setExpr("sum expr");
}

SumExpression::SumExpression(const std::vector <Expression *> &v)
	: v(v)
{
	std::stringstream s;
	s << "Sum expr. Args: " << v.size();
	setExpr(s.str());
}

std::string SumExpression::getValue() const {
	if(v.size() <= 0) {
		return "(0.0)";
	}

	std::stringstream s;
	s << "(";
	s << v[0]->getValue();
	for(int i=1; i<v.size(); ++i)
	{
		s << "+" << v[i]->getValue();
	}
	s << ")";
	return s.str();
}
