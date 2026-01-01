#pragma GCC target("fma")
#include<vector>
#include<stdexcept>
#include<cmath>
#include<memory>
#include<cstring>
#include<iostream>
#include<string>
#include<climits>
#include<immintrin.h>
#include<complex>
#include<random>
#include<algorithm>
#define lf double
#define ull unsigned long long
#define ll long long
#define u32 unsigned
#define int128 __int128_t
#define __AVX2__ 1
#ifdef SIZE
#define LENGTH SIZE
#else
#define LENGTH 2000004
#endif
namespace ERROR{
	class Exception:public std::exception{
		std::string message;
	public:
		Exception(const std::string&c):message(c){}
		const char* what()const noexcept override{return message.c_str();}
	};
	class Div_by_zero:public Exception{
	public:
		Div_by_zero():Exception("Error:Division by zero!"){}
		Div_by_zero(const std::string&s):Exception("Error:"+s+" Division by zero!"){}
	};
	class MLE:public Exception{
	public:
		MLE():Exception("Error:Memory Limit Exceeded!"){}
		MLE(const std::string&s):Exception("Error:"+s+" is Memory Limit Exceeded!"){}
	};
	class Negative:public Exception{
	public:
		Negative(const std::string&s):Exception("Error:"+s+" is negative!"){}
	};
	class Number:public Exception{
	public:
		Number():Exception("Error:Not a valid numeric string!"){}
	};
	class Out_of_range:public Exception{
	public:
		Out_of_range():Exception("Error:Out of range!"){}
	};
}

namespace Transform{
#ifndef __AVX2__
	template<typename T>
	struct fComplex{
		T rez,imz;
		fComplex(){rez=0.0,imz=0.0;}
		fComplex(T x,T y){rez=x,imz=y;}
		template<class C>fComplex(std::complex<C>x):rez(x.real()),imz(x.imag()){}
		fComplex operator+(const fComplex&b)const{return {rez+b.rez,imz+b.imz};}
		fComplex operator-(const fComplex&b)const{return {rez-b.rez,imz-b.imz};}
		fComplex operator*(const fComplex&b)const{return {rez*b.rez-imz*b.imz,rez*b.imz+imz*b.rez};}
		fComplex& operator*=(const fComplex&b){*this=*this*b;return *this;}
		fComplex& operator*(const T&b){return {rez*b,imz*b};}
		fComplex& operator*=(const T&b){rez*=b,imz*=b;return *this;}
		fComplex operator+(const T&x)const{return rez+=x;}
		fComplex conj()const{return {rez,-imz};}
		T real()const{return rez;}
		T imag()const{return imz;}
	};
	using Complex=fComplex<lf>;
	struct FFT{
		const lf pi=3.141592653589793;
		const lf pi2=6.283185307179586;
		std::vector<Complex>fft_a;
		void init(int len){fft_a.resize(len);}
		std::vector<int>rev;std::vector<Complex>omega;
		void fft(int flag,int len){
			for(int i=0;i<len;i++){
				if(i<rev[i]) std::swap(fft_a[i],fft_a[rev[i]]);
			}
			for(int i=2;i<=len;i<<=1){
				int t=len/i;
				for(int j=0;j<len;j+=i){
					for(int k=0;k<i/2;k++){
						int idx=k*t;
						if(flag<0) idx=len-idx;
						if(idx>=len) idx-=len;
						Complex w=omega[idx];
						Complex x=fft_a[j+k];
						Complex y=w*fft_a[j+k+i/2];
						fft_a[j+k]=x+y;
						fft_a[j+k+i/2]=x-y;
					}
				}
			}
			lf inv=1/len;
			if(flag==-1) for(int i=0;i<len;i++) fft_a[i]*=inv;
		}
		void Init(int k){
			rev.resize(1<<k,0);
			int l=1<<k;
			for(int i=0;i<l;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(k-1));
			omega.clear(),omega.resize(l,{0.0,0.0});
			omega[0]={1.0,0.0};
			for(int i=1;i<l;i<<=1) omega[i]={cos(pi2*i/l),sin(pi2*i/l)};
			for(int i=0;i<l;i++) if(i&(i-1)) omega[i]=omega[i&(-i)]*omega[i&(i-1)];
		}
	};
#else
	struct Complex{
		__m128d val;
		Complex()=default;
		Complex(const __m128d&x):val(x){}
		Complex(lf x,lf y):val(_mm_set_pd(y,x)){}
		template<class C>Complex(std::complex<C>x):val(_mm_set_pd(x.imag(),x.real())){}
		Complex operator+(const Complex&x)const{return _mm_add_pd(val,x.val);}
		Complex operator-(const Complex&x)const{return _mm_sub_pd(val,x.val);}
		Complex operator*(const Complex&x)const{return _mm_fmaddsub_pd(_mm_unpacklo_pd(val,val),x.val,_mm_unpackhi_pd(val,val)*_mm_permute_pd(x.val,1));}
		Complex& operator*=(const Complex&b){val=_mm_fmaddsub_pd(_mm_unpacklo_pd(val,val),b.val,_mm_unpackhi_pd(val,val)*_mm_permute_pd(b.val,1));return *this;}
		Complex operator*(lf x)const{return _mm_mul_pd(val,_mm_set1_pd(x));}
		Complex& operator*=(lf x){val=_mm_mul_pd(val,_mm_set1_pd(x));return *this;}
		Complex operator+(lf x)const{return _mm_add_pd(val,_mm_set_pd(0.0,x));}
		Complex conj()const{return Complex(_mm_xor_pd(val,_mm_set_pd(-0.0,0.0)));}
		Complex operator-()const{return _mm_mul_pd(val,_mm_set1_pd(-1.0));}
		lf real()const{return _mm_cvtsd_f64(val);}
		lf imag()const{return _mm_cvtsd_f64(_mm_unpackhi_pd(val,val));}
	};
	struct FFT{
		const lf pi=3.141592653589793;
		const lf pi2=6.283185307179586;
		std::vector<Complex>omega;
		Complex calc(const Complex&a,const Complex&b){return _mm_fmadd_pd(_mm_unpacklo_pd(a.val,a.val),b.val,_mm_unpackhi_pd(a.val,a.val)*_mm_permute_pd(b.val,1));}
		void init(int Len){
			if(Len<=(int)omega.size()<<1) return;
			int k=std::__lg(Len-1);
			omega.resize(1<<k),omega[0]={1.0,0.0};Len=1<<k;
			for(int i=1;i<Len;i<<=1) omega[i]=std::polar(1.0,pi/(i<<1));
			for(int i=1;i<Len;i++) if(i&(i-1)) omega[i]=omega[i&(-i)]*omega[i&(i-1)];
		}
		void dif(std::vector<Complex>&a){
			int len=a.size();
			for(int Len=len>>1,sp=len;Len;sp=Len,Len>>=1){
				for(int i=0;i<Len;i++){auto temp=a[i];a[i]=temp+a[i+Len],a[i+Len]=temp-a[i+Len];}
				for(int blk=sp,o=1;blk<len;blk+=sp,o++){
					for(int i=blk;i<blk+Len;i++){auto t1=a[i],t2=a[Len+i]*omega[o];a[i]=t1+t2,a[Len+i]=t1-t2;}
				}
			}
		}
		void dit(std::vector<Complex>&a){
			int len=a.size();
			for(int Len=1,sp=2;Len!=len;Len=sp,sp<<=1){
				for(int i=0;i<Len;i++){auto temp=a[i];a[i]=temp+a[i+Len],a[i+Len]=temp-a[i+Len];}
				for(int blk=sp,o=1;blk<len;blk+=sp,o++){
					for(int i=blk;i<blk+Len;i++){auto t1=a[i],t2=a[Len+i];a[i]=t1+t2,a[Len+i]=(t1-t2)*omega[o].conj();}
				}
			}
		}
		void mul(std::vector<Complex>&F,std::vector<Complex>&G){//F,G must be resized
			int len=F.size();
			lf inv=1.0/len,_2=inv*0.25;
			F[0]=calc(F[0],G[0])*inv;
			F[1]=F[1]*G[1]*inv;
			for(int st=2,ed=3;st<len;st<<=1,ed<<=1){
				for(int i=st,j=i+st-1;i<ed;i++,j--){
					Complex oi=(F[i]+F[j].conj()),hi=(F[i]-F[j].conj());
					Complex Oi=(G[i]+G[j].conj()),Hi=(G[i]-G[j].conj());
					Complex A=oi*Oi-hi*Hi*((i&1)?-omega[i>>1]:omega[i>>1]),B=Oi*hi+oi*Hi;
					F[i]=(A+B)*_2,F[j]=(A-B).conj()*_2;
				}
			}
		}
	};
#endif
}

namespace IO{
	static const int BUFSIZE=1<<20;
	char ibuf[BUFSIZE],*is=ibuf,*it=ibuf,obuf[BUFSIZE];int cnt=0;
	inline void flush(){
		fwrite(obuf,1,cnt,stdout),cnt=0;
		return;
	}
	inline char get(){
		if(is==it) it=(is=ibuf)+fread(ibuf,1,BUFSIZE,stdin);
		return is==it? EOF:*is++;
	}
	inline void put(char c){
		obuf[cnt++]=c;
		if(cnt==BUFSIZE) flush();
		return;
	}
	struct AutoFlush{~AutoFlush(){flush();}}flusher;
	inline int read(){
		char c=get();int x=0,neg=0;
		while(!isdigit(c)){
			if(c=='-') neg^=1;
			c=get();
		}
		do x=(x<<1)+(x<<3)+(c&15); while(isdigit(c=get()));
		return neg? -x:x;
	}
	inline void c_write(const char*s,const char&c='\0'){
		while(*s) put(*s++);
		if(c) put(c);
	}
	inline void s_write(const std::string&s,const char&c='\0'){
		for(char cc:s) put(cc);
		if(c) put(c);
	}
	template<typename Tp>
	inline void write(Tp x,const char& c='\0'){
		if(x<0) put('-'),x=-x;
		static int top=0,wr[50];
		do wr[++top]=x%10; while(x/=10);
		while(top) put(wr[top--]|48);
		if(c) put(c);
		return;
	}
}

#if __cplusplus<201703L
static constexpr ull p10[]={1,10,100,1000,10000,100000,1000000,10000000,100000000};
#endif
class UnsignedBigInt{
private:
	static const u32 LEN=8,BASE=100000000,FFT_BASE=10000,T=255,INV=64,DEFAULT=1000;
	static const ull MAX=ULLONG_MAX/(BASE+1),LOG=std::__lg(MAX)-1;
	std::unique_ptr<ull[]>num;int len,Max;
#if __cplusplus>=201703L
	static constexpr ull p10[]={1,10,100,1000,10000,100000,1000000,10000000,100000000};
#endif
	ull Pow_10(int p){return p10[p];}
	bool is_zero()const{return len==1&&num[0]==0;}
	void Expand(int nMax){
		if(nMax<=Max) return;
		int _Max=std::max(Max<<1,nMax);
		if(_Max>LENGTH) throw ERROR::MLE("This number");
		auto _num=std::make_unique<ull[]>(_Max);
		if(num) std::copy(num.get(),num.get()+len,_num.get());
		num=std::move(_num);
		Max=_Max;
	}
	void Mul(const UnsignedBigInt&b){
		ull* temp=new ull[len+b.len+5]();
		for(int i=0;i<len;i++){
			ull carry=0,t=num[i];
			for(int j=0;j<b.len;j++){carry+=temp[i+j]+t*b.num[j];temp[i+j]=carry%BASE,carry/=BASE;}
			temp[i+b.len]+=carry;
		}
		Expand(len+b.len+5);len=len+b.len+5;
		std::fill(num.get(),num.get()+len,0);
		for(int i=0;i<len;i++) num[i]=temp[i];
		for(int i=0;i<len;i++) num[i+1]+=num[i]/BASE,num[i]%=BASE;
		while(len>1&&num[len-1]==0) len--;
		delete[] temp;
	}
	ull Get(const UnsignedBigInt&a,int pos)const{
		return 10ull*BASE*(pos+1>=a.len?0:a.num[pos+1])+10ull*a.num[pos]+(pos?a.num[pos-1]:0)/(BASE/10);
	}
	std::pair<UnsignedBigInt,UnsignedBigInt> Simple_Mod(const UnsignedBigInt&b)const{
		if(b.is_zero()) throw ERROR::Div_by_zero();
		if(*this<b) return std::make_pair(0,*this);
		if(*this==b) return std::make_pair(1,0);
		if(b.len<=2){ull q=b.num[0]+b.num[1]*BASE;return std::make_pair(*this/q,*this%q);}
		UnsignedBigInt Q,R(*this);Q.Expand(len-b.len+1),Q.len=len-b.len+1;
		for(int i=len-b.len;i>=0;i--){
			ull q=0;
			auto Sub=[&](){
				ll t=0;
				for(int j=0;j<b.len;j++){
					t=t-q*b.num[j]+R.num[i+j];
					R.num[i+j]=(ull)(t%BASE),t/=BASE;
					if(R.num[i+j]>=BASE) R.num[i+j]+=BASE,t--;
				}
				if(t) R.num[i+b.len]+=(ull)(t);
				Q.num[i]+=q;
			};
			while((q=Get(R,i+b.len-1)/(Get(b,b.len-1)+1))) Sub();
			q=1;
			for(int j=b.len-1;j>=0;j--){if(R.num[j+i]!=b.num[j]&&(q=b.num[j]<R.num[i+j],true)) break;}
			if(q) Sub();
		}
		while(Q.len>1&&Q.num[Q.len-1]==0) Q.len--;
		while(R.len>1&&R.num[R.len-1]==0) R.len--;
		return std::make_pair(Q,R);
	}
	UnsignedBigInt Left(int cnt)const{
		if(cnt<0) throw ERROR::Negative("Left Shift count");
		if(len+cnt>LENGTH) throw ERROR::MLE("Left Shift");
		if(cnt==0||is_zero()) return *this;
		UnsignedBigInt res;
		res.Expand(len+cnt),res.len=len+cnt;
		std::copy(num.get(),num.get()+len,res.num.get()+cnt);
		while(res.len>1&&res.num[res.len-1]==0) res.len--;
		return res;
	}
	UnsignedBigInt Right(int cnt)const{
		if(cnt<0) throw ERROR::Negative("Right Shift count");
		if(cnt>=len) return UnsignedBigInt();
		UnsignedBigInt res;
		res.Expand(len-cnt);res.len=len-cnt;
		std::copy(num.get()+cnt,num.get()+len,res.num.get());
		while(res.len>1&&res.num[res.len-1]==0) res.len--;
		return res;
	}
	UnsignedBigInt Inv(int n)const{
		if(is_zero()) throw ERROR::Div_by_zero();
		if(len<=(int)INV||n<=(int)INV+len){
			UnsignedBigInt a;a.Expand(n+1),a.len=n+1;
			std::fill(a.num.get(),a.num.get()+a.len,0),a.num[n]=1;
			return a.Simple_Mod(*this).first;
		}
		int k=(n-len+5)>>1,kk=k>len?0:len-k;
		UnsignedBigInt t=Right(kk);
		int n1=k+t.len;
		UnsignedBigInt t1=t.Inv(n1);
		UnsignedBigInt res=(t1+t1).Left(n-n1-kk)-(*this*t1*t1).Right(2*(n1+kk)-n);
		return --res;
	}
	std::pair<UnsignedBigInt,UnsignedBigInt> Mod(const UnsignedBigInt&b)const{
		if(*this<b) return std::make_pair(0,*this);
		if(len<=(int)T||b.len<=(int)T) return Simple_Mod(b);
		int Len=len-b.len+5,cnt=Len>b.len?0:b.len-Len;
		UnsignedBigInt tem=b.Right(cnt);
		if(cnt) tem++;
		int inv=Len+tem.len;
		UnsignedBigInt Q=(*this*tem.Inv(inv)).Right(inv+cnt);
		while(Q*b>*this) Q--;
		UnsignedBigInt R=*this-Q*b;
		while(R>=b) Q++,R-=b;
		return std::make_pair(Q,R);
	}
	ull Get_sqrt()const{
		if(len<2) throw ERROR::Out_of_range();
		int top=len-1;ull res=num[top]*BASE+num[top-1];
		return std::sqrt(res)+1;
	}
public:
	UnsignedBigInt():len(1),Max(DEFAULT){
		num=std::make_unique<ull[]>(DEFAULT);
	}
	~UnsignedBigInt()=default;
	UnsignedBigInt(const UnsignedBigInt&b):len(b.len),Max(b.Max){
		num=std::make_unique<ull[]>(Max);
		std::copy(b.num.get(),b.num.get()+len,num.get());
	}
	UnsignedBigInt(UnsignedBigInt&&b)noexcept:num(std::move(b.num)),len(b.len),Max(b.Max){
		b.len=1;
		b.Max=DEFAULT;
		b.num=std::make_unique<ull[]>(DEFAULT);
	}
	UnsignedBigInt(const ull&b):len(0),Max(DEFAULT){
		num=std::make_unique<ull[]>(DEFAULT);
		ull x=b;do{num[len]=x%BASE,len++;}while(x/=BASE);
		while(len>1&&!num[len-1]) len--;
	}
	UnsignedBigInt& operator=(const UnsignedBigInt&b){
		if(this!=&b){
			if(Max<b.len) num.reset(),num=std::make_unique<ull[]>(b.len),Max=b.len;
			std::copy(b.num.get(),b.num.get()+b.len,num.get());len=b.len;
		}
		return *this;
	}
	UnsignedBigInt& operator=(UnsignedBigInt&&b)noexcept{
		if(this!=&b){
			num=std::move(b.num),len=b.len,Max=b.Max;
			b.len=1,b.Max=DEFAULT,b.num=std::make_unique<ull[]>(DEFAULT);
		}
		return *this;
	}
	void init(const u32&size=DEFAULT){
		Expand(size),len=1;
		std::fill(num.get(),num.get()+size,0);
	}
	UnsignedBigInt(const std::string&s){
		if(s.empty()) throw ERROR::Number();
		if(s[0]=='-') throw ERROR::Negative("This number");
		for(int i=0;i<(int)s.size();i++) if(!isdigit(s[i])) throw ERROR::Number();
		len=(s.size()+LEN-1)/LEN;
		if(len>=LENGTH) throw ERROR::MLE("This number");
		Max=std::max((int)DEFAULT,len+5),num=std::make_unique<ull[]>(Max);
		for(int i=(int)s.size()-1,j=0;i>=0;i--,j++){
			num[j/LEN]+=(s[i]-'0')*Pow_10(j%LEN);
		}
		while(len>1&&num[len-1]==0) len--;
	}
	UnsignedBigInt(const char*s){
		u32 Len=std::strlen(s);
		if(Len==0) throw ERROR::Number();
		if(s[0]=='-') throw ERROR::Negative("This number");
		for(int i=0;i<(int)Len;i++) if(!isdigit(s[i])) throw ERROR::Number();
		len=(Len+LEN-1)/LEN;
		if(len>=LENGTH) throw ERROR::MLE("This number");
		Max=std::max((int)DEFAULT,len+5),num=std::make_unique<ull[]>(Max);
		for(int i=(int)Len-1,j=0;i>=0;i--,j++){
			num[j/LEN]+=(s[i]-'0')*Pow_10(j%LEN);
		}
		while(len>1&&num[len-1]==0) len--;
	}
#if __cplusplus>=201703L
	UnsignedBigInt(std::string_view s){
		if(s.empty()) throw ERROR::Number();
		if(s[0]=='-') throw ERROR::Negative("This number");
		for(int i=0;i<(int)s.size();i++) if(!isdigit(s[i])) throw ERROR::Number();
		len=(s.size()+LEN-1)/LEN;
		if(len>=LENGTH) throw ERROR::MLE("This number");
		Max=std::max((int)DEFAULT,len+5),num=std::make_unique<ull[]>(Max);
		for(int i=(int)s.size()-1,j=0;i>=0;i--,j++){
			num[j/LEN]+=(s[i]-'0')*Pow_10(j%LEN);
		}
		while(len>1&&num[len-1]==0) len--;
	}
#endif
	void fread(){
		std::string s;
		char c=IO::get();
		while(!isgraph(c)) c=IO::get();
		while(isgraph(c)&&c!=EOF) s+=c,c=IO::get();
		*this=UnsignedBigInt(s);
	}
	void fwrite(const char&c='\0'){
		IO::write(num[len-1]);
		for(int i=len-2;i>=0;i--){
			char buf[10];
			sprintf(buf,"%08llu",num[i]);
			IO::s_write(buf);
		}
		IO::put(c);
		IO::flush();
	}
	friend std::istream& operator>>(std::istream&in,UnsignedBigInt&x){
		std::string s;in>>s;
		if(in) x=UnsignedBigInt(s);
		return in;
	}
	friend std::ostream& operator<<(std::ostream&out,const UnsignedBigInt&x){
		out<<x.num[x.len-1];
		for(int i=x.len-2;i>=0;i--){
			char buf[10];
			sprintf(buf,"%08llu",x.num[i]);
			out<<buf;
		}
		return out;
	}
	int Two(){
		if(num[0]%2){return 0;}
		*this/=2;
		return Two()+1;
	}
	int Cmp(const UnsignedBigInt&b)const{
		if(len!=b.len) return len>b.len?1:-1;
		for(int i=len-1;i>=0;i--){
			if(num[i]!=b.num[i]) return num[i]>b.num[i]?1:-1;
		}
		return 0;
	}
#if __cplusplus>=202002L
	std::strong_ordering operator<=>(const UnsignedBigInt&b)const{return Cmp(b)<=>0;}
#endif
	bool operator<(const UnsignedBigInt&b)const{return Cmp(b)<0;}
	bool operator>=(const UnsignedBigInt&b)const{return Cmp(b)>=0;}
	bool operator<=(const UnsignedBigInt&b)const{return Cmp(b)<=0;}
	bool operator>(const UnsignedBigInt&b)const{return Cmp(b)>0;}
	bool operator==(const UnsignedBigInt&b)const{return Cmp(b)==0;}
	bool operator!=(const UnsignedBigInt&b)const{return !(*this==b);}
	ull operator[](const int&b)const{return num[b];}
	ull at(const int&b)const{if(b>=Max||b<0) throw ERROR::Out_of_range();return num[b];}
	UnsignedBigInt operator+(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);c+=b;
		return c;
	}
	UnsignedBigInt& operator+=(const UnsignedBigInt&b){
		int n=std::max(len,b.len)+1;Expand(n),len=n;
		for(int i=0;i<n;i++) num[i]+=((i<b.len)?b.num[i]:0);
		for(int i=0;i<n;i++){
			if(num[i]>=BASE){
				num[i+1]+=num[i]/BASE;num[i]%=BASE;
			}
		}
		while(len>1&&num[len-1]==0) len--;
		return *this;
	}
	UnsignedBigInt& operator++(){return *this+=1;}
	UnsignedBigInt operator++(int){UnsignedBigInt res(*this);return *this+=1,res;}
	UnsignedBigInt operator-(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);
		c-=b;return c;
	}
	UnsignedBigInt& operator-=(const UnsignedBigInt&b){
		if(*this<b) throw ERROR::Negative("Result");
		ull t=0;
		for(int i=0;i<len;i++){
			ull sub=((i<b.len)?b.num[i]:0)+t;
			if(num[i]<sub) num[i]=BASE+num[i]-sub,t=1;
			else num[i]-=sub,t=0;
			if(t==0&&i>=b.len) break;
		}
		while(len>1&&num[len-1]==0) len--;
		return *this;
	}
	UnsignedBigInt& operator--(){return *this-=1;}
	UnsignedBigInt operator--(int){UnsignedBigInt res(*this);return *this-=1,res;}
	UnsignedBigInt operator*(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);c*=b;
		return c;
	}
	UnsignedBigInt& operator*=(const UnsignedBigInt&b){
		//https://judge.yosupo.jp/submission/341166
		if(len<=(int)T||b.len<=(int)T){Mul(b);return *this;}
		int k=1,Len=2,n=len,m=b.len;
		while((1<<k)<n+m) k++,Len<<=1;
		if(Len>LENGTH) throw ERROR::MLE("FFT Length");
#ifndef __AVX2__
		Len<<=1,k++;
		if(Len>LENGTH) throw ERROR::MLE("FFT Length");
		Transform::FFT FFT_a,FFT_b;FFT_a.init(Len+1),FFT_b.init(Len+1);
		for(int i=0;i<len;i++) FFT_a.fft_a[i<<1]={(lf)(num[i]%FFT_BASE),0.0},FFT_a.fft_a[i<<1|1]={(lf)(num[i]/FFT_BASE),0.0};
		for(int i=0;i<b.len;i++) FFT_b.fft_a[i<<1]={(lf)(b.num[i]%FFT_BASE),0.0},FFT_b.fft_a[i<<1|1]={(lf)(b.num[i]/FFT_BASE),0.0};
		FFT_a.Init(k),FFT_b.Init(k);
		FFT_a.fft(1,Len);
		FFT_b.fft(1,Len);
		for(int i=0;i<Len;i++) FFT_a.fft_a[i]*=FFT_b.fft_a[i];
		FFT_a.fft(-1,Len);
		Expand(Len);
		std::fill(num.get(),num.get()+len,0);
		for(int i=0;i<Len;i+=2){
			int idx=i>>1;
			__uint128_t t=(ull)(FFT_a.fft_a[i|1].real()+0.5)*FFT_BASE+(ull)(FFT_a.fft_a[i].real()+0.5);
			num[idx]+=t%BASE;
			num[idx+1]+=num[idx]/BASE;
			num[idx]%=BASE;
			num[idx+1]+=t/BASE;
		}
		len=(Len>>1)+1;
		while(len>1&&!num[len-1]) len--;
		return *this;
#else
		if(Len>LENGTH) throw ERROR::MLE("FFT Length");
		Transform::FFT H;H.init(Len);
		std::vector<Transform::Complex>F(len),G(b.len);
		for(int i=0;i<len;i++) F[i]={(lf)(num[i]%FFT_BASE),(lf)(num[i]/FFT_BASE)};
		for(int i=0;i<b.len;i++) G[i]={(lf)(b.num[i]%FFT_BASE),(lf)(b.num[i]/FFT_BASE)};
		F.resize(Len),G.resize(Len);
		H.dif(F),H.dif(G);H.mul(F,G);H.dit(F);
		Expand(Len);
		std::fill(num.get(),num.get()+len,0);
		for(int i=0;i<Len;i++){
			__uint128_t t=(ull)(F[i].imag()+0.5)*FFT_BASE+(ull)(F[i].real()+0.5);
			num[i]+=t%BASE;
			num[i+1]+=num[i]/BASE;
			num[i]%=BASE;
			num[i+1]+=t/BASE;
		}
		len=Len+1;
		while(len>1&&!num[len-1]) len--;
		return *this;
#endif
	}
	UnsignedBigInt& operator/=(const UnsignedBigInt&b){*this=Mod(b).first;return *this;}
	UnsignedBigInt operator/(const UnsignedBigInt&b)const{return Mod(b).first;}
	UnsignedBigInt& operator%=(const UnsignedBigInt&b){*this=Mod(b).second;return *this;}
	UnsignedBigInt operator%(const UnsignedBigInt&b)const{return Mod(b).second;}
	UnsignedBigInt& operator*=(const ull&b){
		if(b<=MAX){
			Expand(len+3),len+=3;
			for(int i=0;i<len;i++){num[i]*=b;}
			for(int i=0;i<len;i++){num[i+1]+=num[i]/BASE,num[i]%=BASE;}
			while(len>1&&!num[len-1]) len--;
		}
		else Mul(UnsignedBigInt(b));
		return *this;
	}
	UnsignedBigInt operator*(const ull&b){
		UnsignedBigInt c(*this);c*=b;
		return c;
	}
	UnsignedBigInt& operator/=(const ull&b){
		if(b==0) throw ERROR::Div_by_zero();
		int128 d=0;
		for(int i=len-1;i>=0;i--){
			d=d*BASE+num[i];num[i]=d/b;d%=b;
		}
		while(len>1&&!num[len-1]) len--;
		return *this;
	}
	UnsignedBigInt operator/(const ull&b)const{
		UnsignedBigInt c(*this);c/=b;
		return c;
	}
	UnsignedBigInt& operator%=(const ull&b){
		if(b==0) throw ERROR::Div_by_zero();
		int128 d=0;
		for(int i=len-1;i>=0;i--){d=d*BASE+num[i];d%=b;}
		return *this=d;
	}
	UnsignedBigInt operator%(const ull&b)const{UnsignedBigInt res(*this);res%=b;return res;}
	UnsignedBigInt& operator<<=(const ull&b){
		UnsignedBigInt base("2");ull p=b;
		for(;p;p>>=1,base.Square()) if(p&1) *this*=base;
		return *this;
	}
	UnsignedBigInt& operator>>=(const ull&b){
		if(b<=10000){
			ull x=b;
			auto Div=[&](int cnt){
				ull d=0;
				for(int i=len-1;i>=0;i--){d=d*BASE+num[i];num[i]=(d>>cnt);d&=((1ull<<cnt)-1);}
				while(len>1&&!num[len-1]) len--;
				x-=cnt;
			};
			while(x>=LOG) Div(LOG);
			Div(x);
		}
		else *this/=UnsignedBigInt("2").pow(b);
		return *this;
	}
	UnsignedBigInt operator>>(const ull&b)const{UnsignedBigInt res(*this);res>>=b;return res;}
	UnsignedBigInt operator<<(const ull&b)const{UnsignedBigInt res(*this);res<<=b;return res;}
	UnsignedBigInt pow(const ull&b)const{
		if(b==0) return UnsignedBigInt("1");
		if(b==1) return *this;
		UnsignedBigInt res("1"),t(*this);ull p=b;
		for(;p;p>>=1){
			if(p&1) res*=t;
			if(p>1) t.Square();
		}
		return res;
	}
	UnsignedBigInt root(int m)const{
		if(m<0) throw ERROR::Negative("Index");
		if(is_zero()) return *this;
		if(m==1) return *this;
		UnsignedBigInt x(std::min(*this,UnsignedBigInt(BASE-1).Left((len+m-1)/m-1))),xx;
		int top=x.len-1;
		int l=0,r=BASE-1;
		while(l<r){
			int mid=(l+r)>>1;
			x.num[top]=mid;
			if(x.pow(m)<=*this) l=mid+1;
			else r=mid;
		}
		x.num[top]=l;
		while(x.len>1&&!x.num[x.len-1]) x.len--;
		xx=(x*(m-1)+*this/x.pow(m-1))/m;
		while(xx<x){
			std::swap(x,xx);
			xx=(x*(m-1)+*this/x.pow(m-1))/m;
		}
		return x;
	}
	UnsignedBigInt sqrt()const{
		if(is_zero()) return *this;
		if(len==1) return (ull)std::sqrt(num[0]);
		if(len==2) return (ull)std::sqrt(num[1]*BASE+num[0]);
		UnsignedBigInt x,xx;x.Expand((len+1)>>1),x.len=(len+1)>>1;
		int top=x.len-1;ull res=Get_sqrt();
		x.num[top]=res;
		while(x.len>1&&!x.num[x.len-1]) x.len--;
		xx=(x+*this/x)/2;
		while(xx<x) std::swap(x,xx),xx=(x+*this/x)/2;
		return x;
	}
	operator std::string()const{
		std::string s="";s+=std::to_string(num[len-1]);
		for(int i=len-2;i>=0;i--){char buf[10];sprintf(buf,"%08llu",num[i]);s+=buf;}
		return s;
	}
	bool True()const{return !is_zero();}
	void test()const{std::cout<<LENGTH<<'\n';}
	void square(){
		if(len<=(int)T){Mul(*this);return;}
		int k=1,Len=2,n=len;
		while((1<<k)<(n<<1)) k++,Len<<=1;
		Len<<=1,k++;
		if(Len>LENGTH) throw ERROR::MLE("FFT Length");
#ifndef __AVX2__
		Transform::FFT FFT_a;FFT_a.init(Len+1);
		for(int i=0;i<len;i++) FFT_a.fft_a[i<<1]={(lf)(num[i]%FFT_BASE),0.0},FFT_a.fft_a[i<<1|1]={(lf)(num[i]/FFT_BASE),0.0};
		FFT_a.Init(k);
		FFT_a.fft(1,Len);
		for(int i=0;i<Len;i++) FFT_a.fft_a[i]*=FFT_a.fft_a[i];
		FFT_a.fft(-1,Len);
		Expand(Len);
		std::fill(num.get(),num.get()+len,0);
		for(int i=0;i<Len;i+=2){
			int idx=i>>1;
			__uint128_t t=(ull)(FFT_a.fft_a[i|1].real()+0.5)*FFT_BASE+(ull)(FFT_a.fft_a[i].real()+0.5);
			num[idx]+=t%BASE;
			num[idx+1]+=num[idx]/BASE;
			num[idx]%=BASE;
			num[idx+1]+=t/BASE;
		}
		len=(Len>>1)+1;
		while(len>1&&!num[len-1]) len--;
#else
		*this*=*this;
#endif
	}
	UnsignedBigInt Square()const{
		UnsignedBigInt res(*this);
		res.square();return res;
	}
};
UnsignedBigInt operator""_UI(const char*literal,size_t len){return UnsignedBigInt(literal);}

namespace Operation{
	UnsignedBigInt Pow(const UnsignedBigInt&a,int p){
		if(p<0) throw ERROR::Negative("Exponent");
		if(p==0) return UnsignedBigInt("1");
		if(p==1) return a;
		UnsignedBigInt res("1"),t(a);
		for(;p;p>>=1){
			if(p&1) res*=t;
			if(p>1) t.Square();
		}
		return res;
	}
	UnsignedBigInt Pow(const UnsignedBigInt&a,int p,const UnsignedBigInt&Mod){
		if(p<0) throw ERROR::Negative("Exponent");
		if(p==0) return UnsignedBigInt("1");
		UnsignedBigInt t(a%Mod);
		if(p==1) return t;
		UnsignedBigInt res("1");
		for(;p;p>>=1){
			if(p&1) res*=t,res%=Mod;
			if(p>1) t.Square(),t%=Mod;
		}
		return res%Mod;
	}
	UnsignedBigInt Fact(int st,int n){
		if(n<=16){
			UnsignedBigInt res=1;
			for(int i=st;i<st+n;i++) res*=i;
			return res;
		}
		int mid=(n+1)/2;
		return Fact(st,mid)*Fact(st+mid,n-mid);
	}
	UnsignedBigInt Fact(int n){return Fact(1,n);}
	UnsignedBigInt Gcd(const UnsignedBigInt&a,const UnsignedBigInt&b){
		UnsignedBigInt c(a),d(b);
		int p=std::min(c.Two(),d.Two());
		while(true){
			int res=c.Cmp(d);
			if(res>0) c-=d,c.Two();
			else if(res<0) d-=c,d.Two();
			else break;
		}
		c<<=p;
		return c;
	}
	UnsignedBigInt Lcm(const UnsignedBigInt&a,const UnsignedBigInt&b){return a*b/Gcd(a,b);}
	UnsignedBigInt Root(const UnsignedBigInt&a,ull p){return a.root(p);}
	UnsignedBigInt Random(int len){
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_int_distribution<int>digit(0,9);
		std::string s="";
		s.push_back('1'+digit(gen)%9);
		for(int i=1;i<len;i++){
			s.push_back('0'+digit(gen));
		}
		return UnsignedBigInt(s);
	}
	UnsignedBigInt Sqrt(const UnsignedBigInt&a){return a.sqrt();}
}

#undef ull
#undef u32
#undef lf
#undef ll
#undef LENGTH
#undef int128