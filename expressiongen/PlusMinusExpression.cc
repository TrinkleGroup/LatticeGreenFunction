#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include "Expression.hh"
#include "PlusMinusExpression.hh"

PlusMinusExpression::PlusMinusExpression() : sign('\0'), v(0) {
	setExpr("plus minus product");
}

PlusMinusExpression::PlusMinusExpression(char sign, const std::vector <Expression *> &v) : sign(sign), v(v)
{
	std::stringstream s;
	s << "plus minus. sign: " << sign << " args: " << v.size();
	setExpr(s.str());
}

std::string PlusMinusExpression::getValue() const {
	if(v.size() <= 0) {
		return "(0.0)";
	}

	std::stringstream s;
	s << "(";
	switch(sign) {
	case '+':
		s << v[0]->getValue();
		break;
	case '-':
		s << "-" << v[0]->getValue();
		break;
	default:
		std::cerr << "ERROR: Invalid PlusMinus!" << std::endl;
		exit(-1);
	}

	for(int i=1; i<v.size(); ++i)
	{
		s << "*" << v[i]->getValue();
	}
	s << ")";
	return s.str();
}
