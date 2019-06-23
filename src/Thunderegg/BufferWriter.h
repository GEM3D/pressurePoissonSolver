/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef THUNDEREGG_BUFFERIO_H
#define THUNDEREGG_BUFFERIO_H
#include "Serializable.h"
#include <cstddef>
#include <iostream>
#include <type_traits>
template <typename T> constexpr bool isSerializable()
{
	return std::is_base_of<Serializable, T>::value;
}
/**
 * @brief Class that is used to help serialize objects into a buffer.
 */
class BufferWriter
{
	private:
	char *buffer = nullptr;
	int   pos    = 0;

	public:
	/**
	 * @brief Create a new BufferWriter with the buffer set to nullptr. This is helpful for
	 * determining the size needed for the buffer.
	 */
	BufferWriter() = default;
	/**
	 * @brief Create a new BufferWriter with given buffer.
	 *
	 * @param buffer the pointer to the beginning of the buffer.
	 */
	BufferWriter(char *buffer)
	{
		this->buffer = buffer;
	}
	/**
	 * @brief get the current position in the buffer
	 *
	 * @return  the current position
	 */
	int getPos()
	{
		return pos;
	}
	/**
	 * @brief  Add object to the buffer.
	 *
	 * @param obj the Serializable object.
	 *
	 * @return  this BufferWriter
	 */
	BufferWriter &operator<<(const Serializable &obj)
	{
		pos += obj.serialize(buffer == nullptr ? nullptr : (buffer + pos));
		return *this;
	}
	/**
	 * @brief Add an object to the buffer.
	 *
	 * @tparam T the type of the object.
	 * @param obj the object. This object must be in serialized form.
	 *
	 * @return  this BufferWriter
	 */
	template <typename T>
	typename std::enable_if<!isSerializable<T>(), BufferWriter>::type &operator<<(const T &obj)
	{
		if (buffer != nullptr) *(T *) (buffer + pos) = obj;
		pos += sizeof(T);
		return *this;
	}
};
/**
 * @brief Class that is used to help read serialized objects from a buffer.
 */
class BufferReader
{
	private:
	char *buffer = nullptr;
	int   pos    = 0;

	public:
	/**
	 * @brief Create a new BufferReader with given buffer.
	 *
	 * @param buffer the pointer to the beginning of the buffer.
	 */
	BufferReader(char *buffer)
	{
		this->buffer = buffer;
	}
	/**
	 * @brief get the current position in the buffer
	 *
	 * @return the current position
	 */
	int getPos()
	{
		return pos;
	}
	/**
	 * @brief Get an object of the buffer.
	 *
	 * @param obj the Serializable object.
	 *
	 * @return  this BufferReader
	 */
	BufferReader &operator>>(Serializable &obj)
	{
		pos += obj.deserialize(buffer + pos);
		return *this;
	}
	/**
	 * @brief Get an object from the buffer.
	 *
	 * @tparam T the type of the object.
	 * @param obj the object. This object must be in serialized form.
	 *
	 * @return  this BufferReader
	 */
	template <typename T>
	typename std::enable_if<!isSerializable<T>(), BufferReader>::type &operator>>(T &obj)
	{
		obj = *(T *) (buffer + pos);
		pos += sizeof(T);
		return *this;
	}
};
#endif
