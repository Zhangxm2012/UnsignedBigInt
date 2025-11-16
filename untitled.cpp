#include<bits/stdc++.h>
#define QWQ "QwQ,Cannot handle negative numbers."
#define ll long long
#define ull unsigned long long
#define Exit(str,val) {cout<<str;exit(val);}
#define Re_val 3221225477
#define LENGTH 700000
#define MLE "Memory Limit Exceeded"
#define lf double
#define llf long double
using namespace std;
const lf pi=acos(-1);
template<typename T>
struct Complex{
	T rez,imz;
	Complex(){rez=0.0,imz=0.0;}
	Complex(lf x,lf y){rez=x,imz=y;}
	Complex operator+(const Complex&b)const{return {rez+b.rez,imz+b.imz};}
	Complex operator-(const Complex&b)const{return {rez-b.rez,imz-b.imz};}
	Complex operator*(const Complex&b)const{return {rez*b.rez-imz*b.imz,rez*b.imz+imz*b.rez};}
	void operator/=(const int&b){rez/=b,imz/=b;}
	void operator*=(const Complex&b){*this=*this*b;}
};
ull Max;
template<typename T>
struct FFT{
	vector<Complex<T>>fft_a;
	void init(int len){fft_a.resize(len);}
	int rev[LENGTH];
	void fft(int flag,int len){
		for(int i=0;i<len;i++){
//			cout<<i<<':'<<rev[i]<<'\n';
			if(i<rev[i]) swap(fft_a[i],fft_a[rev[i]]);
		}
		for(int i=1;i<len;i<<=1){
			Complex<T> w_n={cos(flag*pi/i),sin(flag*pi/i)};
			for(int j=0;j<len;j+=(i<<1)){
				Complex<T>w={1,0};
				for(int k=j;k<j+i;k++){
					Complex<T>x=fft_a[k];
					Complex<T>y=w*fft_a[i+k];
					fft_a[k]=x+y;
					fft_a[i+k]=x-y;
					w*=w_n;
				}
			}
		}
		if(flag==-1) for(int i=0;i<len;i++) fft_a[i]/=len;
	}
	void Init(int k){
		memset(rev,0,sizeof rev);
		int l=1<<k;
		for(int i=0;i<l;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(k-1));//cout<<i<<':'<<rev[i]<<'\n';
	}
};
class UnsignedBigInt{
private:
	static const int LEN=8,BASE=100000000,FFT_BASE=10000;
	ull num[LENGTH];int len;
	int Pow(int a,int p){
		int res=1;
		for(;p;p>>=1,a=a*a) if(p&1) res=res*a;
		return res;
	}
public:
	UnsignedBigInt(){
		memset(num,0,sizeof num);
		len=1;
	}
	UnsignedBigInt(const string&s){
		init();
		int st=0;
		if(s[0]=='-') Exit(QWQ,Re_val);
		len=(s.size()-st+LEN-1)/LEN;
		for(int i=s.size()-1,j=0;i>=st;i--,j++){
			num[j/LEN]+=(s[i]-'0')*Pow(10,j%LEN);
		}
		while(len>1&&!num[len-1]) len--;
	}
	UnsignedBigInt(const UnsignedBigInt&b){
		memset(num,0,sizeof num);
		memcpy(num,b.num,b.len*sizeof(ull));
		len=b.len;
	}
	void init(){
		memset(num,0,sizeof num);
		len=1;
	}
	friend istream& operator>>(istream&in,UnsignedBigInt&x){
		string s;in>>s;
		if(in) x=UnsignedBigInt(s);
		return in;
	}
	friend ostream& operator<<(ostream&out,const UnsignedBigInt&x){
		out<<x.num[x.len-1];
		for(int i=x.len-2;i>=0;i--){
			char buf[10];
			sprintf(buf,"%08lld",x.num[i]);
			out<<buf;
		}
		return out;
	}
	bool operator<(const UnsignedBigInt&b)const{
		if(len!=b.len) return len<b.len;
		for(int i=len-1;i>=0;i--){
			if(num[i]!=b.num[i]) return num[i]<b.num[i];
		}
		return 0;
	}
	bool operator>=(const UnsignedBigInt&b)const{return !((*this)<b);}
	bool operator<=(const UnsignedBigInt&b)const{
		if(len!=b.len) return len<b.len;
		for(int i=len-1;i>=0;i--){
			if(num[i]!=b.num[i]) return num[i]<b.num[i];
		}
		return 1;
	}
	bool operator>(const UnsignedBigInt&b)const{return !((*this)<=b);}
	bool operator==(const UnsignedBigInt&b)const{
		if(len!=b.len) return 0;
		for(int i=0;i<len;i++){
			if(num[i]!=b.num[i]) return 0;
		}
		return 1;
	}
	bool operator!=(const UnsignedBigInt&b)const{return !(*this==b);}
	UnsignedBigInt& operator=(const UnsignedBigInt&b){
		if(b.len>=LENGTH) Exit(MLE,Re_val);
		if(this!=&b){memcpy(num,b.num,b.len*sizeof(ull)),len=b.len;}
		return *this;
	}
	ull operator[](const int&b)const{return num[b];}
	ull at(const int&b)const{if(b>=LENGTH) Exit("RE",Re_val);return num[b];}
	UnsignedBigInt operator-()const{UnsignedBigInt c;return c;}
	UnsignedBigInt operator+(const UnsignedBigInt&b)const{
		UnsignedBigInt c;
		int n=max(len,b.len)+1;c.len=n;
		for(int i=0;i<n;i++) c.num[i]+=num[i]+b.num[i];
		for(int i=0;i<n;i++){
			if(c.num[i]>=BASE){
				c.num[i+1]+=c.num[i]/BASE;c.num[i]%=BASE;
			}
		}
		for(int i=n-1;i>=1;i--){
			if(c.num[i]==0) c.len--;
			else break;
		}
		return c;
	}
	void operator+=(const UnsignedBigInt&b){
		int n=max(len,b.len)+1;len=n;
		for(int i=0;i<n;i++) num[i]+=b.num[i];
		for(int i=0;i<n;i++){
			if(num[i]>=BASE){
				num[i+1]+=num[i]/BASE;num[i]%=BASE;
			}
		}
		for(int i=n-1;i>=1;i--){
			if(num[i]==0) len--;
			else break;
		}
	}
	UnsignedBigInt operator-(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);
		c-=b;return c;
	}
	void operator-=(const UnsignedBigInt&b){
		if(*this<b) Exit(QWQ,Re_val);
		ull t=0;
		for(int i=0;i<len;i++){
			ull sub=((i<b.len)?b.num[i]:0)+t;
			if(num[i]<sub) num[i]=BASE+num[i]-sub,t=1;
			else num[i]-=sub,t=0;
			if(t==0&&i>=b.len) break;
		}
		while(len>1&&num[len-1]==0) len--;
	}
	UnsignedBigInt operator*(const UnsignedBigInt&b)const{
		int k=1,Len=2,n=len,m=b.len;
		while((1<<k)<n+m) k++,Len<<=1;
		Len<<=1,k++;
		FFT<lf>FFT_a,FFT_b;FFT_a.init(Len+1),FFT_b.init(Len+1);
		for(int i=0;i<len;i++) FFT_a.fft_a[i<<1].rez=(lf)(num[i]%FFT_BASE),FFT_a.fft_a[i<<1|1].rez=(lf)(num[i]/FFT_BASE);
		for(int i=0;i<b.len;i++) FFT_b.fft_a[i<<1].rez=(lf)(b.num[i]%FFT_BASE),FFT_b.fft_a[i<<1|1].rez=(lf)(b.num[i]/FFT_BASE);
		if(Len>LENGTH) Exit(MLE,Re_val);
		FFT_a.Init(k),FFT_b.Init(k);
		FFT_a.fft(1,Len);
		FFT_b.fft(1,Len);
		for(int i=0;i<Len;i++) FFT_a.fft_a[i]*=FFT_b.fft_a[i];
		FFT_a.fft(-1,Len);
		UnsignedBigInt c;
		for(int i=0;i<Len;i+=2){
			int idx=i>>1;
			__uint128_t t=(ull)(FFT_a.fft_a[i|1].rez+0.5)*FFT_BASE+(ull)(FFT_a.fft_a[i].rez+0.5);
			c.num[idx]+=t%BASE;
			c.num[idx+1]+=c.num[idx]/BASE;
			c.num[idx]%=BASE;
			c.num[idx+1]+=t/BASE;
		}
		c.len=(Len>>1)+1;
		while(c.len>1&&!c.num[c.len-1]) c.len--;
		return c;
	}
}a,b;
int main(){
//	freopen("P1919_1(1).in","r",stdin);
//	freopen("P1919_1(1).ans","w",stdout);
	cin>>a>>b;
	a=a*b;
	cout<<a<<'\n';
	cout<<Max;
	return 0;
}
/*
99999999
99999999
*/

