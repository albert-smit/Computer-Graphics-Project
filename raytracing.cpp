#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glut.h>
#include "raytracing.h"


//temporary variables
Vec3Df testRayOrigin;
Vec3Df testRayDestination;

//bounding box values
Vec3Df bmin;
Vec3Df bmax;

void boundingBox()
{
	//initialise variables to something reasonable
	bmin = Vec3Df(1000, 1000, 1000);
	bmax = Vec3Df(-1000, -1000, -1000);

	//loop through all vertices
	for (int i = 0; i < MyMesh.vertices.size(); i++)
	{
		//each vertex has 3 coords
		for (int j = 0; j < 3; j++)
		{
			//if coord is less than bmin, make it bmin
			if (MyMesh.vertices[i].p[j] < bmin[j])
			{
				bmin[j] = MyMesh.vertices[i].p[j];
			}
			//same story for bmax but other way around
			if (MyMesh.vertices[i].p[j] > bmax[j])
			{
				bmax[j] = MyMesh.vertices[i].p[j];
			}
		}
	}
}

//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file
	//feel free to replace cube by a path to another model
	//please realize that not all OBJ files will successfully load.
	//Nonetheless, if they come from Blender, they should.
    MyMesh.loadMesh("C:/Users/Vlad/Desktop/raytracing/cube.obj", true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);

	//make a bounding box around the entire model
	boundingBox();
}

/* a = b - c */
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];

#define crossProduct0(a,b,c) \
(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define innerProduct(v,q) \
((v)[0] * (q)[0] + \
	(v)[1] * (q)[1] + \
	(v)[2] * (q)[2])

float rayIntersectsTriangle(const float *p, const float *d,
	const float *v0, const float *v1, const float *v2) {

	float e1[3], e2[3], h[3], s[3], q[3];
	float a, f, u, v;
	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct0(h, d, e2);
	a = innerProduct(e1, h);

	if (a > -0.00001 && a < 0.00001)
		return(-1);

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));

	if (u < 0.0 || u > 1.0)
		return(-1);

	crossProduct0(q, s, e1);
	v = f * innerProduct(d, q);

	if (v < 0.0 || u + v > 1.0)
		return(-1);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	float t = f * innerProduct(e2, q);

	if (t > 0.00001) // ray intersection
		return(t);

	else // this means that there is a line intersection
		// but not a ray intersection
		return (-1);

}

//function that gets the colour
Vec3Df getTriangleColour(int i, Vec3Df ray)
{
	//make the normal
	Vec3Df edge01 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p - MyMesh.vertices[MyMesh.triangles[i].v[0]].p;
	Vec3Df edge02 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p - MyMesh.vertices[MyMesh.triangles[i].v[0]].p;
	Vec3Df n = Vec3Df::crossProduct(edge01, edge02);
	n.normalize();

	Vec3Df normal = n; //Vec3Df(MyMesh.vertices[MyMesh.triangles[i].v[j]].n[0], MyMesh.vertices[MyMesh.triangles[i].v[j]].n[1], MyMesh.vertices[MyMesh.triangles[i].v[j]].n[2]);

	Vec3Df vertexPos = ray; //MyMesh.vertices[MyMesh.triangles[i].v[1]].p;

	//ONLY 1 LIGHT FOR NOW
	Vec3Df lightvector = MyLightPositions[(MyLightPositions.size() - 1)] - vertexPos;
	lightvector.normalize();

	//camera vector
	Vec3Df cameravector = MyCameraPosition - vertexPos;
	cameravector.normalize();

	//"halfway" vector, see wiki page for blinn-phong
	Vec3Df H = lightvector + cameravector;
	H.normalize();

	//exponent value
	float s = 1;

	//blinn-phong = N dot H
	float dotprodblinn = Vec3Df::dotProduct(normal, H);

	//clamp to zero
	if (dotprodblinn < 0)
	{
		dotprodblinn = 0;
	}
	else
	{
		//add exponent to equation
		dotprodblinn = pow(dotprodblinn, s);
	}

	//for diffuse
	float dotprod = Vec3Df::dotProduct(normal, lightvector);

	if (dotprod < 0)
	{
		dotprod = abs(dotprod);
	}

	unsigned int triMat = MyMesh.triangleMaterials.at(i);
	Vec3Df col = MyMesh.materials.at(triMat).Kd();

	//diffuse
	Vec3Df diff = col * dotprod;
	//blinn
	Vec3Df blinn = col * dotprodblinn;

	return blinn + diff;
}

//type for bounding box
template<typename T>
class Ray
{
public:
	Ray(Vec3D<T> orig, Vec3D<T> dir) : orig(orig),
		dir(dir),
		tmin(T(0)),
		tmax(std::numeric_limits<T>::max())
	{
		invdir[0] = T(1) / dir[0];
		invdir[1] = T(1) / dir[1];
		invdir[2] = T(1) / dir[2];
		sign[0] = (invdir[0] < 0);
		sign[1] = (invdir[1] < 0);
		sign[2] = (invdir[2] < 0);
	}
	Vec3D<T> orig, dir; /// ray orig and dir 
	mutable T tmin, tmax; /// ray min and max distances 
	Vec3D<T> invdir;
	int sign[3];
};

//another thing for bounding box
template<typename T>
class Box3 
{ 
public: 
	Box3(Vec3D<T> vmin, Vec3D<T> vmax) 
	{ 
		bounds[0] = vmin; 
		bounds[1] = vmax; 
	} 
	Vec3D<T> bounds[2];

	bool intersect(const Ray<T> &r) const {
		T tmin, tmax, tymin, tymax, tzmin, tzmax;
		tmin = (bounds[r.sign[0]][0] - r.orig[0]) * r.invdir[0];
		tmax = (bounds[1 - r.sign[0]][0] - r.orig[0]) * r.invdir[0];
		tymin = (bounds[r.sign[1]][1] - r.orig[1]) * r.invdir[1];
		tymax = (bounds[1 - r.sign[1]][1] - r.orig[1]) * r.invdir[1];
		if ((tmin > tymax) || (tymin > tmax)) 
			return false;
		if (tymin > tmin) 
			tmin = tymin; 
		if (tymax < tmax) 
			tmax = tymax;
		tzmin = (bounds[r.sign[2]][2] - r.orig[2]) * r.invdir[2]; 
		tzmax = (bounds[1 - r.sign[2]][2] - r.orig[2]) * r.invdir[2]; 
		if ((tmin > tzmax) || (tzmin > tmax)) 
			return false; 
		if (tzmin > tmin) 
			tmin = tzmin; 
		if (tzmax < tmax) 
			tmax = tzmax; 
		if (tmin > r.tmin)
			r.tmin = tmin; 
		if (tmax < r.tmax) 
			r.tmax = tmax; 
		return true;
	}
};



//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	//create ray and box
	Ray<float> ray1(origin, dest);
	Box3<float> box(bmin, bmax);

	//if the ray doesn't intersect the box, return colour (default black, change to something else for debugging)
	if (!(box.intersect(ray1)))
	{
		return Vec3Df(0, 0, 0);
	}

	//initialise to far, this float determines closest distance to origin
	float mindistance = 10000;

	int triangleind = -1;

	//used to get vertexPos, a precise location on triangle
	Vec3Df ray;

	const float *p = origin.p;
	const float *d = dest.p;

	for (int i = 0; i < MyMesh.triangles.size(); i++){
		const float *v0 = MyMesh.vertices[MyMesh.triangles[i].v[0]].p.p;
		const float *v1 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p.p;
		const float *v2 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p.p;

		float t = rayIntersectsTriangle( p, d, v0, v1, v2);

		//t < 0 means behind camera
		if (t > 0)
		{
			//closest triangle
			if (mindistance > t)
			{
				ray = dest - origin;
				ray.normalize();
				ray *= t;
				ray += origin;
				mindistance = t;
				triangleind = i;
			}
		}
	}

	//if there was an intersection with a triangle
	if (triangleind >= 0)
	{
		return getTriangleColour(triangleind, ray);
	}

	return Vec3Df(0, 0, 0);

}


void yourDebugDraw()
{
	//draw open gl debug stuff
	//this function is called every frame

	//as an example: 
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glColor3f(0,1,1);
	glBegin(GL_LINES);
	glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();

}

void yourKeyboardFunc(char t, int x, int y)
{
	// do what you want with the keyboard input t.
	// x, y are the screen position

	//here I use it to get the coordinates of a ray, which I then draw in the debug function.
	produceRay(x, y, testRayOrigin, testRayDestination);

	std::cout<<t<<" pressed! The mouse was in location "<<x<<","<<y<<"!"<<std::endl;
}
