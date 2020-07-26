#include "include\D4D septernion.h"




int main(int argc, char *argv[])
{

	D4D_vec4 v1(0, 0, 1, 0), v2(0, 0, 0, 1);

	D4D_plane p;
	p.set_value(v1, v2);
	
	float angle = 30;
	D4D_septernion sep1;
	sep1.set_value_degree(p, angle/2);
	//float len = sep1.sqr_length();  //Only for verification. len should always equals to 1.

	D4D_vec4 v_ori(1, 0, 0, 0);




	D4D_vec4 result = D4D_rotate(v_ori, sep1);

	return 0;
}
