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

namespace Transform{
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
		void mul(std::vector<Complex>&F,std::vector<Complex>&G){
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
}

using Transform::Complex;
using Transform::FFT;
using namespace std;
int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
	cout.tie(0);
    int n,m;cin>>n>>m;
    int len=1;
    while(len<=n+m) len<<=1;
    vector<Complex>F(len);
    for(int i=0;i<=n;i++){int x;cin>>x;F[i]=Complex(x,0);}
    for(int i=0;i<=m;i++){int x;cin>>x;F[i]=Complex(F[i].real(),x);}
    FFT fft;fft.init(len);
    fft.dif(F);
    fft.mul(F,F);
    fft.dit(F);
	// __bit_ceil
    for(int i=0;i<=n+m;i++) cout<<(int)(F[i].imag()/2+0.5)<<' ';
    return 0;
}