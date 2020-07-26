#pragma once

#include "D4D config.h"

class D4D_vec4
{
	using value_t = real;
public:
	real x, y, z, w;

	//ctor and dtor
	D4D_vec4();
	D4D_vec4(const D4D_vec4& src);
	D4D_vec4(D4D_vec4&& src);
	D4D_vec4(real src_x, real src_y, real src_z, real src_w);

	D4D_vec4& operator=(const D4D_vec4& src);
	D4D_vec4& operator=(D4D_vec4&& src);

	virtual ~D4D_vec4();

	//arithmatic
	D4D_vec4 operator+(const D4D_vec4& src)const;
	D4D_vec4 operator-(const D4D_vec4& src)const;
	D4D_vec4 operator*(const real& src)const ;
	D4D_vec4 operator/(const real& src)const;
	//util
	real sqr_length()const;
	// Length needs to invoke sqrt which may cause perfomance issue. Generally libraries don't provide methods which may cause performance issue. If you need length, do it like #include<math.h> real result = sqrtf( my_vec4.sqr_length()); or sqrt for double.
	//normalization
	D4D_vec4 norm(real epsilon= D4D_default_epsilon())const;
	D4D_vec4& norm_self(real epsilon = D4D_default_epsilon());
	//comparison
	bool operator==(const D4D_vec4& rhs);
	bool operator<(const D4D_vec4& rhs);
	bool is_zero() const;

	//vector products
	real dot(const D4D_vec4& rhs)const;
	static real dot(const D4D_vec4& lhs, const D4D_vec4& rhs);
	
	//Due to I specify i*j*k!=l, and i*j!=k*l, so you should not use the cross product to do anything with septernion.
	//This function only helps calculate the second normal for a 4d plane.
	D4D_vec4 cross_prod(const D4D_vec4& v1, const D4D_vec4& v2);
	static D4D_vec4 cross_prod(const D4D_vec4& v1, const D4D_vec4& v2, const D4D_vec4& v3);

	//spatial relationship
	bool is_parallel(const D4D_vec4& v1, real epsilon = D4D_default_epsilon() * 1000)const;
	bool is_perpendicular(const D4D_vec4& v1, real epsilon = D4D_default_epsilon() * 10)const;
	
};