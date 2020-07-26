#pragma once

#include <utility>

#include "D4D config.h"
#include "D4D vec4.h"


class D4D_plane
{
	using value_t = real;
	enum class error_type
	{
		OK,
		bad_param,
		zero_vector_is_bad_here,
		parallel_vectors_are_bad_here,

		call_func_set_value_first,
	};

public:
	D4D_vec4 main;
	D4D_vec4 minor;
	D4D_vec4 normal_1;
	D4D_vec4 normal_2;
	bool avail;

	//true when 2 vectors in the plane are ready, but normals are not ready.
	bool normals_need_update;

	//ctor dtor
	D4D_plane();
	D4D_plane(const D4D_plane& src);
	D4D_plane(D4D_plane&& src);

	D4D_plane& operator=(const D4D_plane& rhs);
	D4D_plane& operator=(D4D_plane&& rhs);
	
	virtual ~D4D_plane();

	//init
	error_type set_value(const D4D_vec4& v1, const D4D_vec4& v2, real epsilon = D4D_default_epsilon());
	//check
	bool is_avail()const;

	//about the normals
	error_type calc_normals();
	error_type get_normals_from_ref(D4D_vec4& get_1, D4D_vec4& get_2);

};

