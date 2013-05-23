#ifndef RSTREAM_H_
#define RSTREAM_H_

std::streamsize 
rstream::xsputn(const char *s, std::streamsize n)
{
	R::Rprintf( "%.*s", n, s);

    return n;
}

int 
rstream::overflow(int c)
{
	 if (c != EOF) {
		R::Rprintf( "%.1s", &c);
	 }

	 return c;
}

int
rstream::sync(){
	R::R_FlushConsole();
    return 0;
}

void write_msg(std::string const& msg) {
	R::Rprintf(msg.c_str());
}

#endif /* RSTREAM_H_ */
