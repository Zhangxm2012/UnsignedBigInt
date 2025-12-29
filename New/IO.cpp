#include<stdio.h>
#include<string>
namespace IO{
	static const int BUFSIZE=1<<20;
	char ibuf[BUFSIZE],*is=ibuf,*it=ibuf,obuf[BUFSIZE];int cnt=0;
	inline void flush(){
		fwrite(obuf,1,cnt,stdout),cnt=0;
		return;
	}
	inline char get(){
		if(is==it) it=(is=ibuf)+fread(ibuf,1,BUFSIZE,stdin);
		return is==it?EOF:*is++;
	}
	inline void put(char c){
		obuf[cnt++]=c;
		if(cnt==BUFSIZE) flush();
		return;
	}
	struct AutoFlush{~AutoFlush(){flush();}}flusher;
	template<typename Tp>
	inline Tp read(){
		char c=get();Tp x=0,neg=0;
		while(!isdigit(c)){
			if(c=='-') neg^=1;
			c=get();
		}
		do x=(x<<1)+(x<<3)+(c&15);while(isdigit(c=get()));
		return neg?-x:x;
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
		do wr[++top]=x%10;while(x/=10);
		while(top) put(wr[top--]|48);
		if(c) put(c);
		return;
	}
}