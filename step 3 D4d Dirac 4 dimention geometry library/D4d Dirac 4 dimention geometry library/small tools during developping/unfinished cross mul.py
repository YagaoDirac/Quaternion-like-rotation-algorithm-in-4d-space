print("symbol doesn't work correctly. DON'T use this. Not finished yet.")
exit()

##################
origin_indices = [0,1,2,3]
members = ["x","y","z","w"]
result = []
for i in range(len(origin_indices)):
	result.append("")
indices = []
temp_str = []
for i in range(len(origin_indices)-1):
	indices.append([])
	temp_str.append("")

###################

count = 0# to calculate simbol
layer_count = 0

for i in origin_indices:
	layer_count = 0
	temp_str[layer_count] = "(v1."+members[i]
	indices[layer_count] = origin_indices.copy()
	indices[layer_count].remove(i)

	for ii in indices[layer_count]:
		layer_count+=1
		temp_str[layer_count] ="+v2."+members[ii]
		indices[layer_count] = indices[layer_count-1].copy()
		indices[layer_count].remove(ii)

		for iii in indices[layer_count]:
			layer_count+=1
			temp_str[layer_count] ="+v3."+members[iii]+")"
			indices[layer_count] = indices[layer_count-1].copy()
			indices[layer_count].remove(iii)
			count +=1#加count在这儿

			if count%2==1:
				result[indices[layer_count][0]]+="+"
			else:
				result[indices[layer_count][0]]+="-"
			#end if
			#3tabs
			for s in temp_str:
				result[indices[layer_count][0]]+=s
			#print(temp_str)
			#print(result[indices[layer_count][0]])

			layer_count-=1#3tabs
			
		layer_count-=1#2tabs

		

##############

#save to file
with open("cross code.txt","w") as f:
	f.write("r.x = ")
	f.write(result[0])
	f.write(";\n")
	f.write("r.y = ")
	f.write(result[1])
	f.write(";\n")
	f.write("r.z = ")
	f.write(result[2])
	f.write(";\n")
	f.write("r.w = ")
	f.write(result[3])
	f.write(";\n")

