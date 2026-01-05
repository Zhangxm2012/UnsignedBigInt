#include "UnsignedBigInt.cpp"

class SignedBigInt{
private:
	UnsignedBigInt abs;
	bool sign;
public:
	SignedBigInt():abs(),sign(0){}
	~SignedBigInt()=default;
	SignedBigInt(const std::string&value){
		abs=(value.data()+(sign=value.front()=='-')),sign=sign&&(abs.True());
	}
	friend std::istream& operator>>(std::istream&in,SignedBigInt&x){
		std::string s;in>>s;
		if(in) x=SignedBigInt(s);
		return in;
	}
	friend std::ostream& operator<<(std::ostream&out,const SignedBigInt&x){
		if(x.sign) out<<'-';
		out<<x.abs;
		return out;
	}
	SignedBigInt& operator*=(const SignedBigInt&b){
		sign^=b.sign;abs*=b.abs;sign&=abs.True();
		return *this;
	}
};