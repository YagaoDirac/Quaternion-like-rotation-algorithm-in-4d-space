#include <math.h>

#include "../include/D4D vec4.h"

//ctor dtor
D4D_vec4::D4D_vec4() { x = y = z = w = 0.f; };
D4D_vec4::D4D_vec4(const D4D_vec4& src) { x = src.x; y = src.y; z = src.z; w = src.w; };
D4D_vec4::D4D_vec4(D4D_vec4&& src) { x = src.x; y = src.y; z = src.z; w = src.w; };
D4D_vec4::D4D_vec4(real src_x, real src_y, real src_z, real src_w) { x = src_x; y = src_y; z = src_z; w = src_w; }
D4D_vec4& D4D_vec4::operator=(const D4D_vec4& src) { x = src.x; y = src.y; z = src.z; w = src.w; return *this; };
D4D_vec4& D4D_vec4::operator=(D4D_vec4&& src) { x = src.x; y = src.y; z = src.z; w = src.w; return *this; };
D4D_vec4::~D4D_vec4() {};//not very sure whether this one has to be virtual.


//operator for arithmatic.
D4D_vec4 D4D_vec4::operator+(const D4D_vec4 & src)const
{
	return D4D_vec4( this->x + src.x, this->y + src.y, this->z + src.z, this->w + src.w );
}
D4D_vec4 D4D_vec4::operator-(const D4D_vec4 & src)const
{
	return D4D_vec4(this->x - src.x, this->y - src.y, this->z - src.z, this->w - src.w);
}
D4D_vec4 D4D_vec4::operator*(const real & src)const
{
	return D4D_vec4(x * src, y * src, z * src, w * src);
}
D4D_vec4 D4D_vec4::operator/(const real & src)const
{
	return D4D_vec4(x / src, y / src, z / src, w / src);
}
//util
real D4D_vec4::sqr_length() const
{
	return x * x + y * y + z * z + w * w;
}
//normalization
D4D_vec4 D4D_vec4::norm(real epsilon)const
{
	real length_sqr = sqr_length();
	if (length_sqr < epsilon)
		return D4D_vec4();
	
#if __D4D_REAL_FLOAT__ == 1
	real length = sqrtf(length_sqr);
#else
	real length = sqrt(length_sqr);
#endif // __D4D_REAL__
	return *this / length;
}
D4D_vec4 & D4D_vec4::norm_self(real epsilon)
{
	real length_sqr = sqr_length();
	if (length_sqr < epsilon)
	{
		x = y = z = w = 0.f;
		return *this;
	}
	
#if __D4D_REAL_FLOAT__ == 1
	real length = sqrtf(length_sqr);
#else
	real length = sqrt(length_sqr);
#endif // __D4D_REAL__
	x /= length; y /= length; z /= length; w /= length;
	return *this;
}

//comparison
bool D4D_vec4::operator==(const D4D_vec4 & rhs) { return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w; }
bool D4D_vec4::operator<(const D4D_vec4 & rhs) { return x < rhs.x && y < rhs.y && z < rhs.z && w < rhs.w; }
bool D4D_vec4::is_zero() const { return 0.f == x && 0.f == y && 0.f == z && 0.f == w; }

//vector products
real D4D_vec4::dot(const D4D_vec4 & rhs) const { return (x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w); }
real D4D_vec4::dot(const D4D_vec4 & lhs, const D4D_vec4 & rhs) { return lhs.dot(rhs); }

//only for calculating the second normal for plane.
D4D_vec4 D4D_vec4::cross_prod(const D4D_vec4 & v1, const D4D_vec4 & v2)
{
	return D4D_vec4::cross_prod(*this, v1, v2);
}
D4D_vec4 D4D_vec4::cross_prod(const D4D_vec4 & v1, const D4D_vec4 & v2, const D4D_vec4 & v3)
{
	 D4D_vec4 r;
	 real xy = v2.x*v3.y;
	 real xz = v2.x*v3.z;
	 real xw = v2.x*v3.w;
	 
	 real yx = v2.y*v3.x;
	 real yz = v2.y*v3.z;
	 real yw = v2.y*v3.w;

	 real zx = v2.z*v3.x;
	 real zy = v2.z*v3.y;
	 real zw = v2.z*v3.w;

	 real wx = v2.w*v3.x;
	 real wy = v2.w*v3.y;
	 real wz = v2.w*v3.z;

	 r.w = v1.x*(yz - zy) + v1.y*(zx - xz) + v1.z*(xy - yx);
	 r.z = v1.x*(-yw + wy) + v1.y*(xw - wx) + v1.w*(-xy + yx);
	 r.y = v1.x*(zw - wz) + v1.z*(-xw + wx) + v1.w*(xz - zx);
	 r.x = v1.y*(-zw + wz) + v1.z*(yw - wy) + v1.w*(-yz + zy);
	 
	 return r;
}

//spatial relationship
bool D4D_vec4::is_parallel(const D4D_vec4& v1, real epsilon) const
{
	real product_of_2_sqr_length = this->sqr_length()*v1.sqr_length();
	real cross_dot = this->dot(v1);
	real cross_dot_and_then_sqr = cross_dot * cross_dot;

	//the cmath version of abs refuses to work for me..
	real abs_value;
	if (product_of_2_sqr_length > cross_dot_and_then_sqr)
		abs_value = product_of_2_sqr_length - cross_dot_and_then_sqr;
	else 
		abs_value = cross_dot_and_then_sqr - product_of_2_sqr_length;

	return (abs_value < epsilon);
}
bool D4D_vec4::is_perpendicular(const D4D_vec4 & v1, real epsilon) const
{
	real dot_result = this->dot(v1);
	return abs(dot_result < epsilon);
}




