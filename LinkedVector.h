#ifndef LINKEDVECTOR_H_SEEN
#define LINKEDVECTOR_H_SEEN
#include <stdexcept>
#include <vector>
/**
 * @brief This class is is a simple extension of the vector class in the
 * standard library. It can be linked to vectors that are to it's left and
 * right, and it has an at() method that can access the values from those
 * vectors.
 *
 * @tparam T the type that the vector will hold.
 */
template <typename T> class LinkedVector : public std::vector<T>
{
	public:
	typedef int size_type;
	using std::vector<T>::vector;

	/**
	 * @brief the pointer to it's left neighbor
	 */
	LinkedVector<T> *left_nbr_ptr = nullptr;
	/**
	 * @brief the pointer to it's right neighbor
	 */
	LinkedVector<T> *right_nbr_ptr = nullptr;

	/**
	 * @brief Retrieves the value at a given index. If the given index is out of
	 * bounds, it will try to get the value from one of it's neighbors.
	 *
	 * @param i the index to access
	 *
	 * @return a reference to the value at that index.
	 */
	T &at(int i)
	{
		if (i >= 0 && i < (int) this->size()) {
			return (*this)[i];
		} else if (i < 0 && left_nbr_ptr != nullptr) {
			int new_i = i + left_nbr_ptr->size();
			return left_nbr_ptr->at(new_i);
		} else if (i >= (int) this->size() && right_nbr_ptr != nullptr) {
			int new_i = i - (int) this->size();
			return right_nbr_ptr->at(new_i);
		} else {
			throw std::out_of_range("Index is out of bounds");
		}
	}
};

#endif
