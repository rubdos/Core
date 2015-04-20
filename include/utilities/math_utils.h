/****************************************************************************
 *      This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef Y_MATHUTIL_H
#define Y_MATHUTIL_H

#include <cmath>

//#if ( defined(__i386__) || defined(_M_IX86) || defined(_X86_) )
//#define FAST_INT 1
#define _doublemagicroundeps	      (.5-1.4e-11)
	//almost .5f = .5f - 1e^(number of exp bit)
#define _doublemagic			double (6755399441055744.0)
	//2^52 * 1.5,  uses limited precision to floor
//#endif


inline int Round2Int(double val) {
	#ifdef FAST_INT
		val		= val + _doublemagic;
		return ((long*)&val)[0];
	#else
//	#warning "using slow rounding"
		return int (val+_doublemagicroundeps);
	#endif
}

inline int Float2Int(double val) {
	#ifdef FAST_INT
		return (val<0) ?  Round2Int(val+_doublemagicroundeps) :
		   Round2Int(val-_doublemagicroundeps);
	#else
//	#warning "using slow rounding"
		return (int)val;
	#endif
}

inline int Floor2Int(double val) {
	#ifdef FAST_INT
		return Round2Int(val - _doublemagicroundeps);
	#else
//	#warning "using slow rounding"
		return (int)std::floor(val);
	#endif
}

inline int Ceil2Int(double val) {
	#ifdef FAST_INT
		return Round2Int(val + _doublemagicroundeps);
	#else
//	#warning "using slow rounding"	
		return (int)std::ceil(val);
	#endif
}
inline float inRange(float up, float down, float val){
	// verify only values out of the range (up/down)
	// for equal values, no action is necessary
	if (val > up) val = up;
	if (val < down) val = down;
	return val;
}

#endif // Y_MATHUTIL_H
