#ifndef Serializable_H
#define Serializable_H
#include <memory>
class Serializable
{
	public:
	virtual ~Serializable()                   = default;
	virtual int serialize(char *buffer) const = 0;
	virtual int deserialize(char *buffer)     = 0;
};
#endif
