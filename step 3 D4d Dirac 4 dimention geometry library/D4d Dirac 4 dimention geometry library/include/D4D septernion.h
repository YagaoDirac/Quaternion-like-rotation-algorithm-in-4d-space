/*
This class is the core of this lib. I named it Septermion, which means a group of 7 numbers, 1 real, 6 imaginary.
In fact, the septernion is organized like (real part, second imaginary part).
In this scope, 4 imaginary units are involved, I represent them with i, j, k and l.
Generally, an imaginary unit multiply itself is -1, means:
i*i==-1, j*j==-1, k*k==-1, l*l==-1, 
Next, something looks like the commutative, but not exactly. When tow different imaginary units multiply along a spcified sequency, the result is opposite to the different sequency:
i*j == -j*i and so on.
But different from the Quaternion, in Septernion, the 4 imaginary units don't stand in a cycle, means:
i*j*k!=l, and i*j!=k*l.
This causes the cross production for vector4 is not part of this system.

Septernion consist of 1 real scalar and 6 second imaginary parts, and it should look like:
S == r + aij + bik + cil + djk + ejl + fkl
Now you see, if ij==kl, the septernion collapses.
If vector4 == xi + yj + zk + wl
In this way, the result of S multiply vector4 consists of 2 kinds of items.
One is 1-imaginary-unit items, from r*vec4 and when ij*i or ij*j, because iji == -iij == j and ijj == -i.
Other items has 3 different imaginary units, such as when ij*k, it's ijk.
According to the way quaternion pair clamping vector3 makes vector rotated,
(Quaternion*vector3*(conjugate of Quaternion)) >>> rotated vector3
The septernion also does the save thing in the same way.
(Septernion*vector4*(conjugate of Septernion)) >>> rotated vector4

Everything looks so beautiful, but would they work?


[[[delete this part
The 2 imaginary parts are all vector 4. The all 4 elements in each imaginary part are comsidered with its own imaginary unit, something like i or j.
So, an unnormalized Novemion should look somgthing like
Nov == am + bn + c (m and n are imaginary units, c represents some abitrary scalar)
and for am:
am == a*(xi+yj+zk+wl)*m (where i,j,k,l are all another series of imaginary units)
So, the whole Nov should look like:
Nov == ax1im + ay1jm + az1km + aw1lm + bx2in + by2jn + bz2kn + bw2ln + c
a,b,x1,y1,z1,w1,x2,y2,z2,w2,c are scalars, m,n are top level imaginary units, i,j,k,l are bottom level imaginary units.
According to the Quaternion,
i*j*k==l; j*k*l==i; k*l*i==j; l*i*j==k; i*i==j*j==k*k==l*l==1;
i*j==-j*i; and so on.
For the m and n:
m*m==n*n==-1;
But this units system doesn't appear in any part of code.

In this file, the Novemion is organized like (vec4, vec4, real).
I hope it to function just the same way how quaternion works in 3d rotation.
delete this part]]]

Would it work for me?

The only info you should keep in mind is that, the default empty constructor function assigned the real part to 1, and all the imaginary parts to 0.
This looks like the empty Quaternion, which is (1,0,0,0). The 1 is at the real part.
*/


#pragma once

#include <math.h>

#include "D4D config.h"

#include "D4D vec4.h"
#include "D4D plane4d.h"


class D4D_septernion
{
	using real_type = real;
public:
	real ij, ik, il, jk, jl, kl;
	real real_part;

	//It's expensive to calculate this. Notice, when this bool is false, it means the septernion is not sure about if itself is normalized, which means, when it's false, it might have a chance to be normalized. Invoke check_normal() to get the certain result.
	bool normalized;
	
	//ctor
	D4D_septernion();
	D4D_septernion(const D4D_septernion& src);
	D4D_septernion(D4D_septernion&& src);
	D4D_septernion(const D4D_plane &plane, const real degree, real epsilon = D4D_default_epsilon());
	D4D_septernion(real in_r, real in_ij, real in_ik, real in_il, real in_jk, real in_jl, real in_kl);

	D4D_septernion& operator=(const D4D_septernion& rhs);
	D4D_septernion& operator=(D4D_septernion&& rhs);
	
	virtual ~D4D_septernion();

	//arithmatic
	D4D_septernion operator*(const real)const;
	D4D_septernion operator/(const real)const;

	//util
	enum class error_type
	{
		OK,
		septernion_is_too_short,
		//*bad_param,
		//zero_vector_is_bad_here,
		//parallel_vectors_are_bad_here,

		//call_func_set_value_first,*/
	};
	void set_value_Rad(const D4D_plane &plane, const real rad, real epsilon = D4D_default_epsilon());
	void set_value_degree(const D4D_plane &plane, const real degree, real epsilon = D4D_default_epsilon());
	real sqr_length()const;

	bool check_normal(real epsilon= D4D_default_epsilon());//this function calculate the length of the whole novemion. It's expensive.

	//this when
	D4D_septernion norm(error_type* return_error_info = nullptr, real epsilon = D4D_default_epsilon())const;
	error_type norm_self(real epsilon = D4D_default_epsilon());
	D4D_septernion Conjugation()const;
	bool Conjugation_self();





private:
	const char* inventor = "Dirac";
	//It's reported some compilers refuse to compile this line. 
	//If your compiler does refuse, uninstall it.
};

//Notice: if a novem corresponds to a rotation 0f 90 degrees is pass to the function, the result vector rotated 180 degrees. 
D4D_vec4 D4D_rotate(const D4D_vec4& vec, const D4D_septernion& novem);






//cpp

//ctor dtor
/*
D4D_septernion::D4D_septernion() 
{
	real_part = 1;//Notice that, mathmatically, the "empty"septernion, which contains no info, should be (1+0imaginary part).
	ij = ik = il = jk = jl = kl;
	normalized = false; 
};
D4D_septernion::D4D_septernion(const D4D_septernion& src) 
{
	ij = src.ij; ik = src.ik; il = src.il; jk = src.jk; jl = src.jl; kl = src.kl;
	real_part = src.real_part; normalized = src.normalized;
};
D4D_septernion::D4D_septernion(D4D_septernion&& src)
{
	ij = src.ij; ik = src.ik; il = src.il; jk = src.jk; jl = src.jl; kl = src.kl;
	real_part = src.real_part; normalized = src.normalized;
};
inline D4D_septernion::D4D_septernion(const D4D_plane & plane, const real degree, real epsilon)
{
	set_value_Rad(plane, degree / 180.f*D4D_pi(), epsilon);
	return;
}
*/
/**/
/*
inline D4D_septernion::D4D_septernion(real in_r, real in_ij, real in_ik, real in_il, real in_jk, real in_jl, real in_kl)
{
	real_part = in_r;
	ij = in_ij; ik = in_ik; il = in_il; jk = in_jk; jl = in_jl; kl = in_kl;
	normalized = false;
}

D4D_septernion& D4D_septernion::operator=(const D4D_septernion& rhs) 
{
	ij = rhs.ij; ik = rhs.ik; il = rhs.il; jk = rhs.jk; jl = rhs.jl; kl = rhs.kl;
	real_part = rhs.real_part; normalized = rhs.normalized;
	return *this;
};
D4D_septernion& D4D_septernion::operator=(D4D_septernion&& rhs) 
{
	ij = rhs.ij; ik = rhs.ik; il = rhs.il; jk = rhs.jk; jl = rhs.jl; kl = rhs.kl;
	real_part = rhs.real_part; normalized = rhs.normalized;
	return *this;
};

D4D_septernion::~D4D_septernion() {};//not very sure whether this one has to be virtual.

//arithmatic
inline D4D_septernion D4D_septernion::operator*(const real src) const
{
	return D4D_septernion(real_part*src, ij*src, ik*src, il*src, jk*src, jl*src, kl*src);
}
inline D4D_septernion D4D_septernion::operator/(const real src) const
{
	return D4D_septernion(real_part/src, ij/src, ik/src, il/src, jk/src, jl/src, kl/src);

}

//util
inline void D4D_septernion::set_value_Rad(const D4D_plane & plane, const real rad, real epsilon)
{
	if (!plane.avail)
	{
		real_part = 1.f; ij = ik = il = jk = jl = kl = 0.f;
		return;
	}

	real_part = -(plane.normal_1.dot(plane.normal_2));//and this part should always near equal to 0.
	if (real_part > -epsilon && real_part < epsilon) real_part = 0;

	ij = plane.normal_1.x *plane.normal_2.y + plane.normal_1.y*plane.normal_2.x;
	ik = plane.normal_1.x *plane.normal_2.z + plane.normal_1.z*plane.normal_2.x;
	il = plane.normal_1.x *plane.normal_2.w + plane.normal_1.w*plane.normal_2.x;
	jk = plane.normal_1.y *plane.normal_2.z + plane.normal_1.z*plane.normal_2.y;
	jl = plane.normal_1.y *plane.normal_2.w + plane.normal_1.w*plane.normal_2.y;
	kl = plane.normal_1.z *plane.normal_2.w + plane.normal_1.w*plane.normal_2.z;

	real_part = real_part + (1 - real_part)*cosf(rad);
	real sin_value = sinf(rad);
	ij = ij * sin_value;
	ik = ik * sin_value;
	il = il * sin_value;
	jk = jk * sin_value;
	jl = jl * sin_value;
	kl = kl * sin_value;
}
inline void D4D_septernion::set_value_degree(const D4D_plane & plane, const real degree, real epsilon)
{
	set_value_Rad(plane, degree / 180.f*D4D_pi(), epsilon);
	return;
}




real D4D_septernion::sqr_length()const 
{
	return real_part * real_part + ij * ij + ik * ik + il * il + jk * jk + jl * jl + kl * kl;
}
inline bool D4D_septernion::check_normal(real epsilon) 
{
	real diff = sqr_length() - epsilon;
	normalized = ((diff > -epsilon) && (diff < epsilon));
	return normalized;
}
inline D4D_septernion D4D_septernion::norm(error_type* return_error_info, real epsilon) const
{
	real length_sqr = sqr_length();
	if (length_sqr < epsilon)
	{
		if (return_error_info)
			*return_error_info = error_type::septernion_is_too_short;
		return D4D_septernion();
	}

#if __D4D_REAL_FLOAT__ == 1
	real length = sqrtf(length_sqr);
#else
	real length = sqrt(length_sqr);
#endif // __D4D_REAL__
	return *this / length;
}
inline D4D_septernion::error_type D4D_septernion::norm_self(real epsilon)
{
	real length_sqr = sqr_length();
	if (length_sqr < epsilon)
	{
		return error_type::septernion_is_too_short;
	}

#if __D4D_REAL_FLOAT__ == 1
	real length = sqrtf(length_sqr);
#else
	real length = sqrt(length_sqr);
#endif // __D4D_REAL__
	return error_type::OK;
}

inline D4D_septernion D4D_septernion::Conjugation() const
{
#if __D4D_REAL_FLOAT__ == 1
	return D4D_septernion(real_part, ij*-1.f, ik*-1.f, il*-1.f, jk*-1.f, jl*-1.f, kl*-1.f);
#else
	return D4D_septernion(real_part, ij*-1., ik*-1., il*-1., jk*-1., jl*-1., kl*-1.);
#endif // __D4D_REAL__
}
inline bool D4D_septernion::Conjugation_self()
{
#if __D4D_REAL_FLOAT__ == 1
	ij = ij * -1.f; ik = ik * -1.f; il = il * -1.f; jk = jk * -1.f; jl = jl * -1.f; kl = kl * -1.f;
#else
	ij = ij * -1.; ik = ik * -1.; il = il * -1.; jk = jk * -1.; jl = jl * -1.; kl = kl * -1.;
#endif // __D4D_REAL__
}

*/



//out of the class

//because ijk!=l, ij!=kl, so, these kind of septernion*vector4 doesn't give any valid result inside what could be represent in this lib.
/*
D4D_vec4 D4Df_mul(const D4D_septernion& sep, D4D_vec4 vec);
D4D_vec4 operator*(const D4D_septernion& sep, D4D_vec4 vec);
*/



//The rotation carried out by this function is double the degree of sep. If the Septernion is set to 90 degrees, the rotation is 180 degrees.
D4D_vec4 D4D_rotate(const D4D_vec4& vec, const D4D_septernion& sep);
/*

D4D_vec4 D4Df_rotate(const D4D_vec4& vec, const D4D_septernion& sep)
{
	//This function calculates set * vec * (conjugate of sep)
	/* This is readme.txt   lmao..
		S == r+v2, in this formula, r means the real part. It's a scalar. v2 means the imaginary part. Notice that, v2 mean, all the items inside are followed by 2 different imaginary units.
		so, conj(S) == r-v2
		V is the vector. But in here, V == v1. The number of 1 means all the 4 items in it are followed by only 1 imaginary units.
		result == S*V*conj(S)
		== (r+v2)*v1*(r-v2)
		== r*r*v1 + r*v2*v1 - r*v1*v2 - v2*v1*v2
		
		r*r*v1 is a D4D_vec4. It's OK.
		
		But what is v2*v1.
		The 3 cases are, aij in v2 multiplied by bi in v1. This produces a*b*i*j*i. a and b are real number, scalar, and they are welcome. What is iji? Due to ij==-ji, iji==- jii. Then ii ==-1, -jii ==-j*-1, it's j.
		The same for ij*j, it's -i.
		The 2 cases upon result in something with only one imaginary unit.
		The third case is somgthing like, aij in v2 multiplied by bk(or bl) in v1. It's ab*ijk. Because I specified the ijk!=l, so nothing could get simplified here. The results s*ck. 
		But they don't s*ck that hard as imagined. Notice, ijk==-ikj==kij. And for the 2 cases before, ij*i==iji==-iij==j and also, i*ij==-j. Means for v2*v1 and v1*v2, the single unit parts are opposite( left == right*-1 ), and the triple units part are the same. (For more info, check the folder of "About the math" out)
		So, the main fomular should derived to 
		result == S*V*conj(S)
		== r*r*v1 + r*v2*v1 - r*v1*v2 - v2*v1*v2
		== rr*v2 + 2r*( single unit part of(v2*v1)) - v2*v1*v2
		
		Now let's move on to the last part. v2*v1*v2
		According to the derivation, all the triple unit parts turned out to be 0. This is great.
		So, v2*v1*v2 == single unit part of (v2*v1*v2)

		Now the conclusion is very clear, S*V*conj(S) has only single unit parts. Means it's a vector4, a D4D_vec4.
		The job left is only completing all the details.
	*/
/*
	//temp
	real r(sep.real_part), rr(r*r);
	real a(sep.ij), b(sep.ik), c(sep.il), d(sep.jk), e(sep.jl), f(sep.kl), aa(a*a), bb(b*b), cc(c*c), dd(d*d), ee(e*e), ff(f*f);
	real x(vec.x), y(vec.y), z(vec.z), w(vec.w);
	real ax(a*x), ay(a*y), az(a*z), aw(a*w),
		bx(b*x), by(b*y), bz(b*z), bw(b*w),
		cx(c*x), cy(c*y), cz(c*z), cw(c*w),
		dx(d*x), dy(d*y), dz(d*z), dw(d*w),
		ex(e*x), ey(e*y), ez(e*z), ew(e*w),
		fx(f*x), fy(f*y), fz(f*z), fw(f*w);


	//r*r*v1
	D4D_vec4 rrv1 = vec * rr; 

	//2r*(single unit part of(v2*v1))
	D4D_vec4 v2v1;
	v2v1.x = -ay - bz - cw;
	v2v1.y = ax - dz - ew;
	v2v1.z = bx + dy - fw;
	v2v1.w = cx + ey + fz;
	D4D_vec4 part2 = v2v1 * (2.f*r);

	//v2*v1*v2
	D4D_vec4 v2v1v2;
	v2v1v2.x = x * (aa + bb + cc - dd - ee - ff)
		+ 2.f * a * (-dz - ew)
		+ 2.f * b * (dy - fw)
		+ 2.f * c * (fz + ey);
	v2v1v2.y = y * (aa + dd + ee - bb - cc - ff)
		+ 2.f * a * (bz + cw)
		+ 2.f * d * (bx - fw)
		+ 2.f * e * (cx + fz);
	v2v1v2.z = z * (bb + dd + ff - aa - cc - ee)
		+ 2.f * b * (ay + cw)
		+ 2.f * d * (-ax + ew)
		+ 2.f * f * (cx + ey); 
	v2v1v2.w = w * (cc + ee + ff - aa - bb - dd)
		+ 2.f * c * (ay + bz)
		+ 2.f * e * (-ax + dz)
		+ 2.f * f * (-bx - dy);

	return rrv1 + part2 + v2v1v2;

	/*You know what, I'm absolutely not gonna do anything that has anything to do with the Quaternion style rotation in 5d or anything similar to this.*/
	/*
}
*/