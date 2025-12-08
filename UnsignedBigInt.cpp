#include<bits/stdc++.h>
#define debug cout<<"Here\n"
#define QWQ "QwQ,Cannot handle negative numbers."
#define ll long long
#define ull unsigned long long
#define Exit(str,val) {cout<<str;exit(val);}
#define Re_val 3221225477
#define LENGTH 4000
#define MLE "Memory Limit Exceeded"
#define lf double
#define llf long double
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
	Complex conj(){return {rez,-imz};}
};
template<typename T>
struct FFT{
	const lf pi=3.141592653589793;
	const lf pi2=6.283185307179586;
	vector<Complex<T>>fft_a;
	void init(int len){fft_a.resize(len);}
	vector<int>rev;vector<Complex<lf>>omega;
	void fft(int flag,int len){
		for(int i=0;i<len;i++){
			if(i<rev[i]) swap(fft_a[i],fft_a[rev[i]]);
		}
		for(int i=2;i<=len;i<<=1){
			int t=len/i;
			for(int j=0;j<len;j+=i){
				for(int k=0;k<i/2;k++){
					int idx=k*t;
					if(flag<0) idx=len-idx;
					if(idx>=len) idx-=len;
					Complex<T>w=omega[idx];
					Complex<T>x=fft_a[j+k];
					Complex<T>y=w*fft_a[j+k+i/2];
					fft_a[j+k]=x+y;
					fft_a[j+k+i/2]=x-y;
				}
			}
		}
		if(flag==-1) for(int i=0;i<len;i++) fft_a[i]/=len;
	}
	void Init(int k){
		rev.resize(1<<k,0);
		int l=1<<k;
		for(int i=0;i<l;i++) rev[i]=(rev[i>>1]>>1)|((i&1)<<(k-1));
		omega.clear(),omega.resize(l);
		omega[0]={1.0,0.0};
		for(int i=1;i<l;i<<=1) omega[i]={cos(pi2*i/l),sin(pi2*i/l)};
		for(int i=0;i<l;i++) omega[i]=omega[i&(-i)]*omega[i&(i-1)];
	}
};
class UnsignedBigInt{
private:
	static const int LEN=8,BASE=100000000,FFT_BASE=10000,T=255;
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
		if(s[0]=='-') Exit(QWQ,Re_val);
		for(int i=0;i<(int)s.size();i++){if(!isdigit(s[i])) Exit("Not a valid numeric string",Re_val);}
		init();
		int st=0;
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
	UnsignedBigInt Left(int cnt)const{
		if(cnt<=0) Exit("Left Shift count is negative!",Re_val);
		if(len+cnt>LENGTH) Exit("Left Shift:Memory Limit Exceeded",Re_val);
		UnsignedBigInt res;
		memmove(res.num+cnt,num,len*sizeof(ull)),memset(res.num,0,cnt*sizeof(ull));
		res.len=len+cnt;return res;
	}
	UnsignedBigInt Right(int cnt)const{
		if(cnt<=0) Exit("Right Shift count is negative!",Re_val);
		if(cnt>=len) return UnsignedBigInt();
		UnsignedBigInt res;
		memmove(res.num,num+cnt,(len-cnt)*sizeof(ull));
		memset(res.num+(len-cnt),0,cnt*sizeof(ull));
		res.len=len-cnt;
		while(res.len>1&&res.num[res.len-1]==0) res.len--;
		return res;
	}
	UnsignedBigInt Inv(int n)const{
		if(*this==0) Exit("Error:Division by zero!",Re_val);
		if(min(len,n-len)<=16){
			UnsignedBigInt a;a.len=n+1;a.num[len]=1;
		}
	}
	bool Cmp_eq(const UnsignedBigInt&a,const UnsignedBigInt&b,int last)const{
		if(last+b.len>a.len) return 0;
		if(a.num[last+b.len]!=0) return 1;
		for(int i=b.len-1;i>=0;i--){
			if(a.num[last+i]!=b.num[i]){return a.num[last+i]>b.num[i];}
		}
		return 1;
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
		if(*this<b) Exit(QWQ,Re_val);
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
	UnsignedBigInt operator*=(const UnsignedBigInt&b){
		if(len<=T||b.len<=T){Mul(b);return *this;}
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
		memset(num,0,sizeof(int)*(Len));
		for(int i=0;i<Len;i+=2){
			int idx=i>>1;
			__uint128_t t=(ull)(FFT_a.fft_a[i|1].rez+0.5)*FFT_BASE+(ull)(FFT_a.fft_a[i].rez+0.5);
			num[idx]+=t%BASE;
			num[idx+1]+=num[idx]/BASE;
			num[idx]%=BASE;
			num[idx+1]+=t/BASE;
		}
		len=(Len>>1)+1;
		while(len>1&&!num[len-1]) len--;
		return *this;
	}
	UnsignedBigInt operator*=(const ull&b){
		if(b<=184467438892){
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
		if(b==0) Exit("Error:Division by zero!",Re_val);
		ull d=0;
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
	UnsignedBigInt operator%=(const ull&b){return *this-=(*this/b*b);}
	UnsignedBigInt operator%(const ull&b)const{return *this-(*this/b*b);}
	UnsignedBigInt operator<<=(const ull&b){
		ull x=b;
		while(x>=37){
			len+=5;
			for(int i=0;i<len;i++) num[i]<<=37;
			for(int i=0;i<len;i++) num[i+1]+=num[i]/BASE,num[i]%=BASE;
			while(len>1&&!num[len-1]) len--;
			x-=37;
		}
		len+=5;
		for(int i=0;i<len;i++) num[i]<<=x;
		for(int i=0;i<len;i++) num[i+1]+=num[i]/BASE,num[i]%=BASE;
		while(len>1&&!num[len-1]) len--;
		return *this;
	}
	UnsignedBigInt operator>>=(const ull&b){
		ull x=b;
		while(x>=37){
			ull d=0;
			for(int i=len-1;i>=0;i--){d=d*BASE+num[i];num[i]=(d>>37);d&=((1ull<<37)-1);}
			x-=37;
		}
		ull d=0;
		for(int i=len-1;i>=0;i--){d=d*BASE+num[i];num[i]=(d>>x);d&=((1ull<<x)-1);}
		return *this;
	}
	
};
namespace Operation{
	UnsignedBigInt pow(const UnsignedBigInt&a,int p){
		if(p==0) return UnsignedBigInt("1");
		if(p==1) return a;
		UnsignedBigInt res("1"),t(a);
		for(;p;p>>=1){
			if(p&1) res*=t;
			if(p>1) t*=t;
		}
		return res;
	}
	UnsignedBigInt fact(int st,int n){
		if(n<=16){
			UnsignedBigInt res=1;
			for(int i=st;i<st+n;i++) res*=i;
			return res;
		}
		int mid=(n+1)/2;
		return fact(st,mid)*fact(st+mid,n-mid);
	}
	UnsignedBigInt gcd(const UnsignedBigInt&a,const UnsignedBigInt&b){
		UnsignedBigInt c(a),d(b);
		int p=min(c.Two(),d.Two());
		while(true){
			int res=c.Cmp(d);
			if(res>0) c-=d,c.Two();
			else if(res<0) d-=c,d.Two();
			else break;
		}
		for(int i=1;i<=p;i++) c*=2;
		return c;
	}
}
//using namespace Operation;
int main(){
	UnsignedBigInt a,b;
	cin>>a>>b;
	cout<<(a+b)<<'\n';
	if(a<b) cout<<'-'<<(b-a);
	else cout<<a-b;
	cout<<'\n';
	cout<<(a*b)<<'\n';
	auto res=a.Simple_Mod(b);
	cout<<res.first<<'\n'<<res.second;
	return 0;
}
/*
1145141919810
114514
*/
