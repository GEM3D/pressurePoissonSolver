#ifndef ZEROSUM_H
#define ZEROSUM_H
class ZeroSum
{
	private:
	static bool flag;

	public:
	ZeroSum() {}
	static void setTrue();
	static bool isSet();
};
bool ZeroSum::flag = false;
void ZeroSum::setTrue() { flag = true; }
bool ZeroSum::isSet() { return flag; }
#endif
