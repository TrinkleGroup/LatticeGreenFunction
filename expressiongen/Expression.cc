#include "Expression.hh"

Expression::Expression()
{
}

Expression::Expression(std::string &s) : expr(s) { }

Expression::Expression(const char *s) : expr(s) { }

const std::string &Expression::getExpr() const {
	return expr;
}

void Expression::setExpr(const std::string &s) {
	expr = s;
}
