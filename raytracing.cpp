#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glut.h>
#include "raytracing.h"


//temporary variables
Vec3Df testRayOrigin;
Vec3Df testRayDestination;

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
}
/* a = b - c */
#define vector(a,b,c) \
	(a)[0] = (b)[0] - (c)[0];	\
	(a)[1] = (b)[1] - (c)[1];	\
	(a)[2] = (b)[2] - (c)[2];

#define crossProduct(a,b,c) \
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

	crossProduct(h, d, e2);
	a = innerProduct(e1, h);

	if (a > -0.00001 && a < 0.00001)
		return(-1);

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));

	if (u < 0.0 || u > 1.0)
		return(-1);

	crossProduct(q, s, e1);
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

Vec3Df getTriangleColour(int i)
{
	unsigned int triMat = MyMesh.triangleMaterials.at(i);
	Vec3Df col = MyMesh.materials.at(triMat).Kd();
	return col;
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	std::vector<float> depthbuffer;
	//std::vector<int> triangleind;
	float mindistance = 10000;
	int triangleind = -1;
	//float distance = (origin - dest).getLength();

	for (int i = 0; i < MyMesh.triangles.size(); i++){

		const float *p = origin.p;
		const float *d = dest.p;
		const float *v0 = MyMesh.vertices[MyMesh.triangles[i].v[0]].p.p;
		const float *v1 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p.p;
		const float *v2 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p.p;

		float t = rayIntersectsTriangle( p, d, v0, v1, v2);

		if (t > 0)
		{
			if (mindistance > t)
			{
				mindistance = t;
				triangleind = i;
			}
			/*depthbuffer.push_back(t);
			return Vec3Df(1, 1, 1);*/
		}
	}
	if (triangleind >= 0)
	{
		return getTriangleColour(triangleind);
	}

	return Vec3Df(0, 0, 0);


	//Vec3Df L = origin;
	//Vec3Df r = dest - origin;
	//float cosangle = Vec3Df::dotProduct(L, r) / (r.getLength() * L.getLength());
	//float angle = acos(cosangle);
	//float d = L.getLength() * sin(angle);
	//float sphereradius = 1.5;

	//if (d > sphereradius)
	//	return Vec3Df(0, 0, 0);

	//float thc = sqrt(sphereradius * sphereradius + d * d);
	//float tca = L.getLength() * cosangle;
	//float t0 = tca - thc;

	////get r vector, make it length t0
	//Vec3Df t0vec = r;
	//t0vec.normalize();
	//t0vec *= t0;

	//Vec3Df P = origin - t0vec;

	//Vec3Df vertexPos = P;

	//Vec3Df normal = P;
	//normal.normalize();

	//Vec3Df lightvector = origin - vertexPos;
	//lightvector.normalize();
	////lambertian shading value
	//float dotprod = Vec3Df::dotProduct(normal, lightvector);

	////clamp to zero when negative
	//if (dotprod < 0)
	//{
	//	dotprod = abs(dotprod);
	//}

	////turn lambertian value into a vector, red colour
	//Vec3Df diffusecolour(1, 0, 0);
	//diffusecolour *= dotprod;
	//return diffusecolour;

	//return Vec3Df(dest[0],dest[1],dest[2]);
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
