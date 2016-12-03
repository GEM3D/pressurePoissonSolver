#ifndef LINKEDVECTOR_H_SEEN 
#define LINKEDVECTOR_H_SEEN
#include <vector>
#include <stdexcept>
template <typename T> class LinkedVector : public std::vector<T> {
    public:
    typedef int size_type;
    using std::vector<T>::vector;
    LinkedVector<T>* left_nbr_ptr = nullptr;
    LinkedVector<T>* right_nbr_ptr = nullptr;
    T& at(int i) {
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
