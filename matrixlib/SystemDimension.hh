#ifndef _SYSTEM_DIMENSION_H_
#define _SYSTEM_DIMENSION_H_

#define DEFAULT_CARTDIM 3u
#define CARTDIM SystemDimension::singleton().getSysdim()

class SystemDimension {
private:
	unsigned int sysdim;
	SystemDimension(unsigned int u);
	static SystemDimension *sd;
public:
	static void initialize(unsigned int u);
	static SystemDimension& singleton();
	virtual unsigned int getSysdim() const;
};

#endif
