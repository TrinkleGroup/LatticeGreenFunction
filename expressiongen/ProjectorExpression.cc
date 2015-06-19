#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "ProjectorExpression.hh"

ProjectorExpression::ProjectorExpression() : proj('\0') {
	setExpr("Projector expr");
}

ProjectorExpression::ProjectorExpression(char proj) : proj(proj) {
	std::stringstream s;
	s << "Projector: " << proj;
	setExpr(s.str());
}

std::string ProjectorExpression::getValue() const {
	switch(proj) {
	case 'A':
		return "(dynmat.getAcousticProjector())";
	case 'a':
		return "(dynmat.getAcousticProjectorT())";
	case 'O':
		return "(dynmat.getOpticalProjector())";
	case 'o':
		return "(dynmat.getOpticalProjectorT())";
	case 'U':
		return "(dynmat.getAcousticFullP())";
	case 'u':
		return "(dynmat.getAcousticFullPT())";
	case 'L':
		return "(dynmat.getOpticalFullP())";
	case 'l':
		return "(dynmat.getOpticalFullPT())";
	default:
		std::cerr << "ERROR: Invalid Projector!" << std::endl;
		exit(-1);
	}
}
