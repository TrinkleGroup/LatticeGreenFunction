#include "SystemDimension.hh"
#include <cstdlib>

SystemDimension * SystemDimension::sd = NULL;

SystemDimension::SystemDimension(unsigned int u) : sysdim(u) { }

void SystemDimension::initialize(unsigned int u) {
	if(SystemDimension::sd != NULL) { 
		delete SystemDimension::sd;
	}
	SystemDimension::sd = new SystemDimension(u);
}

SystemDimension& SystemDimension::singleton() {
	if(SystemDimension::sd == NULL) { 
		SystemDimension::initialize(DEFAULT_CARTDIM);
	}
	return *SystemDimension::sd;
}

unsigned int SystemDimension::getSysdim() const {
	return sysdim;
}
