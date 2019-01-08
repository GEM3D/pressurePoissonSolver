#ifndef IFACETYPE_H
#define IFACETYPE_H
/*
enum class IfaceType {
	normal,
	coarse_to_coarse,
	fine_to_coarse_0,
	fine_to_coarse_1,
	fine_to_coarse_2,
	fine_to_coarse_3,
	fine_to_fine_0,
	fine_to_fine_1,
	fine_to_fine_2,
	fine_to_fine_3,
	coarse_to_fine_0,
	coarse_to_fine_1,
	coarse_to_fine_2,
	coarse_to_fine_3
};
inline IfaceType operator+(const IfaceType &a, const int &b)
{
	return static_cast<IfaceType>(static_cast<int>(a) + b);
}
*/
class IfaceType
{
	private:
	/**
	 * @brief the value of the enum
	 */
	char val = -1;
	char orthant = -1;

	public:
	// enum definitions
	static constexpr char normal   = 0b000;
	static constexpr char coarse_to_coarse   = 0b001;
	static constexpr char fine_to_coarse  = 0b010;
	static constexpr char fine_to_fine  = 0b011;
	static constexpr char coarse_to_fine = 0b100;
	IfaceType()=default;
	IfaceType(const char val)
	{
		this->val = val;
	}
	IfaceType(const char val, const char orthant)
	{
		this->val = val;
		this->orthant = orthant;
	}
    void set(int blah){}
    char toInt(){return val;}
    char getOrthant(){return orthant;}
    void setOrthant(const char orthant){this->orthant=orthant;}
    bool operator <(const IfaceType&b)const{
        return std::tie(val,orthant)<std::tie(b.val,b.orthant);
    }
};

#endif
