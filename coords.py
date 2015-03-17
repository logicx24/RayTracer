
import sys

def range_xyz(verts):
	max_x = 0;
	max_y = 0;
	max_z = 0;
	min_x = 0;
	min_y = 0;
	min_z = 0;
	for vert in verts:
		if vert[0] > max_x:
			max_x = vert[0]
		if vert[1] > max_y:
			max_y = vert[1]
		if vert[2] > max_z:
			max_z = vert[2]

		if vert[0] < min_x:
			min_x = vert[0]
		if vert[1] < min_y:
			min_y = vert[1]
		if vert[2] < min_z:
			min_z = vert[2]

	return ((max_x, max_y, max_z), (min_x, min_y, min_z))

def get_verts(file):
	obj_file = open(file);
	lines = obj_file.read()
	lines = lines.split("\n")
	lines = [l for l in lines if l]
	lines = [l for l in lines if l[0] == 'v']
	#lines = [line for line in lines if lines[0] == "v"]
	verts = []
	for l in lines:
		l = l.split()[1:]
		verts.append((float(l[0]), float(l[1]), float(l[2])))
	return verts

def main():
	obj_file = sys.argv[1]
	print(obj_file)
	verts = get_verts(obj_file)
	print(range_xyz(verts))