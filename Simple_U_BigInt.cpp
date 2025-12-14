#include<iostream>
#include<vector>
#include<cmath>
#include<cstring>
#include<climits>
#define ll long long
#define ull unsigned long long
#define Exit(str,val) {cout<<str;exit(val);}
#define MLE "Memory Limit Exceeded"
#define lf double
#define llf long double
#ifdef ONLINE_JUDGE
#define Re_val 0
#define Div_val 0
#else 
#define Re_val 3221225477
#define Div_val 3221225620
#endif
#ifdef SIZE
#define LENGTH SIZE
#else
#define LENGTH 10000
#endif
using namespace std;
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
	inline void s_write(const string&s,const char&c='\0'){
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
class UnsignedBigInt{
private:
	static const int LEN=8,BASE=100000000,FFT_BASE=10000,T=255,INV=64;
	static const ull MAX=ULLONG_MAX/(BASE+1),LOG=__lg(MAX)-1;
	ull num[LENGTH];int len;
	int Pow(int a,int p){
		int res=1;
		for(;p;p>>=1,a=a*a) if(p&1) res=res*a;
		return res;
	}
	void Mul(const UnsignedBigInt&b){
		ull* temp=new ull[len+b.len+5]();
		for(int i=0;i<len;i++){
			ull carry=0,t=num[i];
			for(int j=0;j<b.len;j++){carry+=temp[i+j]+t*b.num[j];temp[i+j]=carry%BASE,carry/=BASE;}
			temp[i+b.len]+=carry;
		}
		memset(num,0,sizeof(ull)*(len));len=len+b.len+5;
		for(int i=0;i<len;i++) num[i]=temp[i];
		for(int i=0;i<len;i++) num[i+1]+=num[i]/BASE,num[i]%=BASE;
		while(len>1&&num[len-1]==0) len--;
		delete[] temp;
	}
	ull Get(const UnsignedBigInt&a,int pos)const{
		return 10ull*BASE*(pos+1>=a.len?0:a.num[pos+1])+10ull*a.num[pos]+(pos?a.num[pos-1]:0)/(BASE/10);
	}
	pair<UnsignedBigInt,UnsignedBigInt> Simple_Mod(const UnsignedBigInt&b)const{
		if(b==0) Exit("Error:Division by zero!",Re_val);
		if(*this<b) return make_pair(0,*this);
		if(*this==b) return make_pair(1,0);
		if(b.len<=2){ull q=b.num[0]+b.num[1]*BASE;return make_pair(*this/q,*this%q);}
		UnsignedBigInt Q,R(*this);Q.len=len-b.len+1;
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
		return make_pair(Q,R);
	}
public:
	UnsignedBigInt(){
		memset(num,0,sizeof num);
		len=1;
	}
	UnsignedBigInt(const string&s){
		if(s[0]=='-') Exit("QwQ,Cannot handle negative numbers.",Re_val);
		for(int i=0;i<(int)s.size();i++){if(!isdigit(s[i])) Exit("Not a valid numeric string",Re_val);}
		init();
		int st=0;
		len=(s.size()-st+LEN-1)/LEN;
		if(len>=LENGTH) Exit(MLE,Re_val);
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
	UnsignedBigInt(const ull&b){
		memset(num,0,sizeof num);len=0;
		ull x=b;do{num[len]=x%BASE,len++;}while(x/=BASE);
		while(len>1&&!num[len-1]) len--;
	}
	void init(){
		memset(num,0,sizeof num);
		len=1;
	}
	void fread(){
		string s;
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
	int Two(){
		if(num[0]%2){return 0;}
		*this/=2;
		return Two()+1;
	}
	int Cmp(const UnsignedBigInt&b){
		if(len!=b.len) return len>b.len?1:-1;
		for(int i=len-1;i>=0;i--){
			if(num[i]!=b.num[i]) return num[i]>b.num[i]?1:-1;
		}
		return 0;
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
			sprintf(buf,"%08llu",x.num[i]);
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
	UnsignedBigInt operator+(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);c+=b;
		return c;
	}
	UnsignedBigInt operator+=(const UnsignedBigInt&b){
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
		return *this;
	}
	UnsignedBigInt operator++(){return *this+=1;}
	UnsignedBigInt operator++(int){UnsignedBigInt res(*this);return *this+=1,res;}
	UnsignedBigInt operator-(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);
		c-=b;return c;
	}
	UnsignedBigInt operator-=(const UnsignedBigInt&b){
		if(*this<b) Exit("QwQ,Cannot handle negative numbers.",Re_val);
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
	UnsignedBigInt operator--(){return *this-=1;}
	UnsignedBigInt operator--(int){UnsignedBigInt res(*this);return *this-=1,res;}
	UnsignedBigInt operator*(const UnsignedBigInt&b)const{
		UnsignedBigInt c(*this);c*=b;
		return c;
	}
	UnsignedBigInt operator*=(const UnsignedBigInt&b){Mul(b);return *this;}
	UnsignedBigInt operator/=(const UnsignedBigInt&b){*this=Simple_Mod(b).first;return *this;}
	UnsignedBigInt operator/(const UnsignedBigInt&b)const{return Simple_Mod(b).first;}
	UnsignedBigInt operator%=(const UnsignedBigInt&b){*this=Simple_Mod(b).second;return *this;}
	UnsignedBigInt operator%(const UnsignedBigInt&b)const{return Simple_Mod(b).second;}
	UnsignedBigInt operator*=(const ull&b){
		if(b<=MAX){
			len+=3;
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
	UnsignedBigInt operator/=(const ull&b){
		if(b==0) Exit("Error:Division by zero!",Div_val);
		__int128_t d=0;
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
	UnsignedBigInt operator%=(const ull&b){
		if(b==0) Exit("Error:Division by zero!",Div_val);
		__int128_t d=0;
		for(int i=len-1;i>=0;i--){d=d*BASE+num[i];d%=b;}
		return *this=d;
	}
	UnsignedBigInt operator%(const ull&b)const{UnsignedBigInt res(*this);res%=b;return res;}
	UnsignedBigInt operator<<=(const ull&b){
		UnsignedBigInt base("1");ull p=b;
		for(;p;p>>=1,base*=2) if(p&1) *this*=base;
		return *this;
	}
	UnsignedBigInt operator>>=(const ull&b){
		ull x=b;
		auto Div=[&](int cnt){
			ull d=0;
			for(int i=len-1;i>=0;i--){d=d*BASE+num[i];num[i]=(d>>cnt);d&=((1ull<<cnt)-1);}
			while(len>1&&!num[len-1]) len--;
			x-=cnt;
		};
		while(x>=LOG) Div(LOG);
		Div(x);
		return *this;
	}
	UnsignedBigInt operator>>(const ull&b)const{UnsignedBigInt res(*this);res>>=b;return res;}
	UnsignedBigInt operator<<(const ull&b)const{UnsignedBigInt res(*this);res<<=b;return res;}
	operator string()const{
		string s="";
		for(int i=len-1;i>=0;i--){char buf[10];sprintf(buf,"%08llu",num[i]);s+=buf;}
		return s;
	}
	bool True()const{return !(len==1&&num[0]==0);}
};
namespace Operation{
	UnsignedBigInt Pow(const UnsignedBigInt&a,int p){
		if(p==0) return UnsignedBigInt("1");
		if(p==1) return a;
		UnsignedBigInt res("1"),t(a);
		for(;p;p>>=1){
			if(p&1) res*=t;
			if(p>1) t*=t;
		}
		return res;
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
	UnsignedBigInt Gcd(const UnsignedBigInt&a,const UnsignedBigInt&b){
		UnsignedBigInt c(a),d(b);
		int p=min(c.Two(),d.Two());
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
}
using namespace Operation;
#undef ll
#undef ull
#undef Exit
#undef MLE
#undef lf
#undef llf
#undef Re_val
#undef Div_val
#undef LENGTH