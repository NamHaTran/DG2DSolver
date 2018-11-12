#ifndef DGPOSTPROCESSLIB_H_INCLUDED
#define DGPOSTPROCESSLIB_H_INCLUDED
namespace debugTool
{
	void checkElemSurPt(int ipoin);
	void checkPtsSurPt(int ipoin);
	void checkElemsSurElem(int ielem);
	void checkElemInfor(int elem);
	void checkPointValue(int element);
}

namespace DG2Matlab
{
	void createMatlabCode();
	void exportData();
}
#endif // DGPOSTPROCESSLIB_H_INCLUDED