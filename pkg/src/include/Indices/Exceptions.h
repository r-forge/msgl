/* 
 * File:   Exception.h
 * Author: Martin Vincent
 *
 * Created on 28. december 2010, 22:45
 */

#ifndef EXCEPTION_H
#define	EXCEPTION_H

#include <iostream>
using namespace std;

#include <string>

class Exception {

	friend ostream & operator<< (ostream & theStream, Exception & exception) {
			theStream << exception.msg << endl;
			return theStream;
		}

public:
	const string msg;

	Exception(string msg) : msg(msg) {};
};

class GenralException : public Exception {
public:
	GenralException(string msg) : Exception("Genral Exception: "+msg) {};
};

#define ASSERT(condition, msg) if(!(condition)) throw GenralException(msg);

#endif	/* EXCEPTION_H */

