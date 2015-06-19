#ifndef _PROJECTOR_EXPRESSION_HH_
#define _PROJECTOR_EXPRESSION_HH_

#include <string>
#include "Expression.hh"

class ProjectorExpression: public Expression
{
private:
	char proj;

public:
	ProjectorExpression();
	ProjectorExpression(char proj);

	virtual std::string getValue() const;
};

#endif
