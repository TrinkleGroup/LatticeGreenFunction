#include <string>
#include <sstream>
#include "DynMatExpression.hh"

DynMatExpression::DynMatExpression() : n(-1) {
	setExpr("Dyn mat expr");
}

DynMatExpression::DynMatExpression(int n) : n(n) {
	std::stringstream s;
	s << "Dyn mat expr. Arg: " << n;
	setExpr(s.str());
}

std::string DynMatExpression::getValue() const {
	std::stringstream s;
	s << "(dynmat.getOrder(" << n << ",k))";
	return s.str();
}
