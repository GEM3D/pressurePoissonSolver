#ifndef BUFFERWRITER_H
#define BUFFERWRITER_H
#include "Serializable.h"
#include <cstddef>
#include <iostream>
#include <type_traits>
template <typename T> constexpr bool isSerializable()
{
	return std::is_base_of<Serializable, T>::value;
}
class BufferWriter
{
	private:
	char *buffer = nullptr;
	int   pos    = 0;

	public:
    BufferWriter()=default;
	BufferWriter(char *buffer) { this->buffer = buffer; }
	int           getPos() { return pos; }
	BufferWriter &operator<<(const Serializable &obj)
	{
		pos += obj.serialize(buffer == nullptr ? nullptr : (buffer + pos));
		return *this;
	}
	template <typename T>
	typename std::enable_if<!isSerializable<T>(), BufferWriter>::type &operator<<(const T &obj)
	{
		if (buffer != nullptr) *(T *) (buffer + pos) = obj;
		pos += sizeof(T);
		return *this;
	}
};
class BufferReader
{
	private:
	char *buffer = nullptr;
	int   pos    = 0;

	public:
    BufferReader()=default;
	BufferReader(char *buffer) { this->buffer = buffer; }
	int           getPos() { return pos; }
	BufferReader &operator>>(Serializable &obj)
	{
		pos += obj.deserialize(buffer + pos);
		return *this;
	}
	template <typename T>
	typename std::enable_if<!isSerializable<T>(), BufferReader>::type &operator>>(T &obj)
	{
		obj = *(T *) (buffer + pos);
		pos += sizeof(T);
		return *this;
	}
};
#endif
