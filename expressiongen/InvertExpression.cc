#include <string>
#include "InvertExpression.hh"

InvertExpression::InvertExpression() : exp(NULL) {
	setExpr("Invert expr");
}

InvertExpression::InvertExpression(Expression *exp) : exp(exp) {
	setExpr("Invert: \"" + exp->getExpr() + "\"");
}

std::string InvertExpression::getValue() const {
	return "(" + exp->getValue() + ".inverse())";
}
