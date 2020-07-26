#include "../include/D4D plane4d.h"




//ctor dtor
D4D_plane::D4D_plane() { avail = false; normals_need_update = false; };
D4D_plane::D4D_plane(const D4D_plane& src)
{
	main = src.main; minor = src.minor;
	avail = src.avail;
	normals_need_update = src.normals_need_update;
};
D4D_plane::D4D_plane(D4D_plane&& src)
{
	main = src.main; minor = src.minor;
	avail = src.avail;
	normals_need_update = src.normals_need_update;
};

D4D_plane& D4D_plane::operator=(const D4D_plane& rhs)
{
	main = rhs.main; minor = rhs.minor;
	avail = rhs.avail;
	normals_need_update = rhs.normals_need_update;
	return *this;
};
D4D_plane& D4D_plane::operator=(D4D_plane&& rhs)
{
	main = rhs.main; minor = rhs.minor;
	avail = rhs.avail;
	normals_need_update = rhs.normals_need_update;
	return *this;
};

D4D_plane::~D4D_plane() {};//not very sure whether this one has to be virtual.

//init

D4D_plane::error_type D4D_plane::set_value(const D4D_vec4 & v1, const D4D_vec4 & v2, real epsilon)
{
	if (v1.is_zero() || v2.is_zero())return error_type::zero_vector_is_bad_here;

	if (v1.is_parallel(v2))return error_type::parallel_vectors_are_bad_here;

	main = v1.norm();
	minor = v2 - main * (v2.dot(main));//get temp_minor to the perpendicular
	minor.norm_self();
	avail = true;
	normals_need_update = true; //not very sure about this line...
	return error_type::OK;
}
//check
inline bool D4D_plane::is_avail() const { return avail; }

inline D4D_plane::error_type D4D_plane::calc_normals()
{
	if (!is_avail())return error_type::call_func_set_value_first;
	if (!normals_need_update) return error_type::OK;

	D4D_vec4 temp(1.f, 0.f, 0.f, 0.f);
	//func tries (1,0,0,0) first. If this is parallel to any of the 2 vectors, tests(0,1,0,0).
	if (temp.is_parallel(main) || temp.is_parallel(minor)) temp = D4D_vec4(0.f, 1.f, 0.f, 0.f);
	//If the attempt fails again, then tests(0,0,1,0).
	if (temp.is_parallel(main) || temp.is_parallel(minor)) temp = D4D_vec4(0.f, 0.f, 1.f, 0.f);

	//Now, the temp vector is good to go. Bases on it, func calculates the first normal.
	//Notice that the main and minor vectors are already normalized.
	temp = temp - main * (temp.dot(main));
	temp = temp - minor * (temp.dot(minor));
	normal_1 = temp.norm();
	normal_2 = D4D_vec4::cross_prod(main, minor, normal_1);//I'm not sure if this is normalized already.
	normal_2.norm_self();
}
inline D4D_plane::error_type D4D_plane::get_normals_from_ref(D4D_vec4 & get_1, D4D_vec4 & get_2)
{
	if (!is_avail())return error_type::call_func_set_value_first;

	if (normals_need_update) calc_normals();
	get_1 = normal_1; get_2 = normal_2;
	return error_type::OK;
}