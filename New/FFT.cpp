#pragma GCC target("fma")
#include<bits/stdc++.h>
#include<immintrin.h>
#include<memory>
#define lf double
#define ll long long
#define __AVX2__ 1
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