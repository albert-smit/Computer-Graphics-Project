#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#endif
#include "raytracing.h"
#include <limits>


//temporary variables
Vec3Df testRayOrigin;
Vec3Df testRayDestination;

//bounding box values
Vec3Df bmin;
Vec3Df bmax;

//Turn features on/off
const bool shadowFlag = true;
const bool reflectionFlag = true;
const bool refractionFlag = true;
const bool blinnPhongFlag = true;
const bool diffuseFlag = true;

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
	MyMesh.loadMesh("objects/coloredCube.obj", true);
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

//check if a ray intersects a triangle, 
//intersect returns the intersection point
float rayIntersectsTriangle(const float *p, const float *d,
	const float *v0, const float *v1, const float *v2, Vec3Df* intersect) {

	float e1[3], e2[3], h[3], s[3], q[3];
	float a, f, u, v;
	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct0(h, d, e2);
	a = innerProduct(e1, h);

	//if (a > -0.00001 && a < 0.00001)
	if (a == 0)
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

	if (t > 0) // ray intersection
	{

		*intersect = (1 - u - v)*Vec3Df(v0[0], v0[1], v0[2]) + u*Vec3Df(v1[0], v1[1], v1[2]) + v*Vec3Df(v2[0], v2[1], v2[2]); // intersect point of ray and triangle
		return(t);
	}

	else // this means that there is a line intersection
		// but not a ray intersection
		return (-1);

}


//check if a point is shaded i.e. no direct light
bool shadow(Vec3Df origin, Vec3Df dest)
{

	Vec3Df intersect;

	float *p = origin.p;
	float *d = dest.p;

	for (int i = 0; i < MyMesh.triangles.size(); i++){
		float *v0 = MyMesh.vertices[MyMesh.triangles[i].v[0]].p.p;
		float *v1 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p.p;
		float *v2 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p.p;

		float t = rayIntersectsTriangle(p, d, v0, v1, v2, &intersect);

		//intersect means above 0, rest is to filter out noise
		if (t > 0.0001)
		{
			return true;
		}

	}

	return false;
}

float getBlinnPhong(Vec3Df normal, Vec3Df H, float exp) {
	//blinn-phong = N dot H
	float blinnPhong = Vec3Df::dotProduct(normal, H);

	//clamp to zero
	if (blinnPhong < 0)
	{
		blinnPhong = 0;
	}
	else
	{
		//add exponent to equation
		blinnPhong = pow(blinnPhong, exp);
	}

	return blinnPhong;
}

float getDiffusion(Vec3Df normal, Vec3Df lightvector) {
	float diffusion = Vec3Df::dotProduct(normal, lightvector);

	if (diffusion < 0)
	{
		//this might need to change
		diffusion = 0.1;
	}

	return diffusion;
}

Vec3Df getRefraction(Vec3Df normal, Vec3Df ray, Vec3Df vertexPos) {
	Vec3Df result = Vec3Df(0, 0, 0);

	//variable to check if the light goes in or out the transparent object
	float inOut = Vec3Df::dotProduct(normal, ray);

	//when inOut < 0 the light goes into the sphere
	if (inOut < 0){

		//the fraction coefficient for glass
		float fractionCoGlass = 0.66;

		//calculate the the ray inside the object
		float dotPro = Vec3Df::dotProduct(normal, -vertexPos);
		float partUnderRoot = sqrtf(1 - fractionCoGlass * fractionCoGlass * (1 - dotPro*dotPro));
		if (partUnderRoot > 0){
			Vec3Df RayInside = (fractionCoGlass * dotPro - partUnderRoot) * normal - (fractionCoGlass * -ray);
		}
		// else the light goes out the sphere
	}
	else{

		//the fraction coefficient for glass
		float fractionCoGlass = 0.66;

		//calculate the the ray inside the object
		float dotPro = Vec3Df::dotProduct(-normal, -ray);
		float partUnderRoot = sqrtf(1 - fractionCoGlass * fractionCoGlass * (1 - dotPro*dotPro));
		if (partUnderRoot > 0){
			Vec3Df RayOutside = (fractionCoGlass * dotPro - partUnderRoot) * normal - (fractionCoGlass * -ray);
		}
	}
	// end refraction

	return result;
}

//function that gets the colour of reflection ray
//essentially a copy of getTriangleColour without reflections
Vec3Df getReflectionColour(int triangleind, Vec3Df ray, Vec3Df origin, int currentLightPos, Vec3Df brightness)
{
	Vec3Df result = Vec3Df(0, 0, 0);
	Vec3Df result2 = Vec3Df(0, 0, 0);

	Vec3Df vertexPos = ray;

	//make the normal
	Vec3Df edge01 = MyMesh.vertices[MyMesh.triangles[triangleind].v[1]].p - MyMesh.vertices[MyMesh.triangles[triangleind].v[0]].p;
	Vec3Df edge02 = MyMesh.vertices[MyMesh.triangles[triangleind].v[2]].p - MyMesh.vertices[MyMesh.triangles[triangleind].v[0]].p;
	Vec3Df normal = Vec3Df::crossProduct(edge01, edge02);
	normal.normalize();

	//light vector
	Vec3Df lightvector;

	//camera vector
	Vec3Df cameravector;

	//"halfway" vector, see wiki page for blinn-phong
	Vec3Df H;

	//get colour of material
	unsigned int triMat = MyMesh.triangleMaterials.at(triangleind);
	Vec3Df col = MyMesh.materials.at(triMat).Kd() * 0.5;

	float blinnPhongResult = 0;
	float diffuseResult = 0;

	lightvector = MyLightPositions[currentLightPos] - vertexPos;
	lightvector.normalize();

	cameravector = origin - vertexPos;
	cameravector.normalize();

	H = lightvector + cameravector;
	H.normalize();

	// Add Blinn-Phong effect
	if (blinnPhongFlag) {
		blinnPhongResult = getBlinnPhong(normal, H, 2);

		result2 += col * blinnPhongResult;
	}

	// Add diffuse effect
	if (diffuseFlag) {
		diffuseResult = getDiffusion(normal, lightvector);

		result2 += col * diffuseResult;
	}


	//check for shadow
	if (shadowFlag) {
		bool shaded = shadow(ray, MyLightPositions[currentLightPos]);
		if (shaded)
		{
			float shadowdepth = 0.2;

			result2 *= shadowdepth;
		}
	}

	result += result2;
	result *= brightness;
	return result;
}

Vec3Df getReflection(Vec3Df cameraPos, Vec3Df selectedPos, Vec3Df normal, int currentLightPos, Vec3Df currentResult) {
	Vec3Df reflectionColour = Vec3Df(0, 0, 0);
	//"camera" vector to origin
	Vec3Df V = cameraPos - selectedPos;
	V.normalize();

	float cosalpha = Vec3Df::dotProduct(V, normal);

	//FOR REFLECTION
	//if "camera" is on wrong side of surface (normal pointing other way, over 90 degrees) then cosalpha < 0
	if (cosalpha > 0) {
		Vec3Df r = ((2 * cosalpha) * normal) - V;
		r.normalize();

		//this will be the intersection point of reflection vector
		Vec3Df intersect;
		float mindistance = 10000;
		int triangleind = -1;

		Vec3Df reflray;

		const float *p = selectedPos.p;
		const float *d = r.p;

		//find closest triangle by looping through all
		for (int o = 0; o < MyMesh.triangles.size(); o++){

			float *v0 = MyMesh.vertices[MyMesh.triangles[o].v[0]].p.p;
			float *v1 = MyMesh.vertices[MyMesh.triangles[o].v[1]].p.p;
			float *v2 = MyMesh.vertices[MyMesh.triangles[o].v[2]].p.p;

			float t = rayIntersectsTriangle(p, d, v0, v1, v2, &intersect);

			//t < 0 means no intersect
			//0.0001 because slight noise filtering
			if (t > 0.0001)
			{
				//closest triangle
				if (mindistance > t)
				{
					reflray = intersect;
					mindistance = t;
					triangleind = o;
				}
			}


		}

		//if there was an intersection with a triangle
		if (triangleind >= 0)
		{
			reflectionColour = getReflectionColour(triangleind, reflray, selectedPos, currentLightPos, currentResult);

			////get colour of material
			unsigned int triMatr = MyMesh.triangleMaterials.at(triangleind);
			Vec3Df colr = MyMesh.materials.at(triMatr).Ks() * 0.5;
            reflectionColour *= colr;
		}

	}//end of the reflection part
	return reflectionColour;
}

//called when ray from origin intersects a triangle
//function that gets the colour
Vec3Df getTriangleColour(int i, Vec3Df ray, Vec3Df origin)
{
	Vec3Df result = Vec3Df(0, 0, 0);
	Vec3Df result2 = Vec3Df(0, 0, 0);

	Vec3Df vertexPos = ray;

	//make the normal
	Vec3Df edge01 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p - MyMesh.vertices[MyMesh.triangles[i].v[0]].p;
	Vec3Df edge02 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p - MyMesh.vertices[MyMesh.triangles[i].v[0]].p;
	Vec3Df normal = Vec3Df::crossProduct(edge01, edge02);
	normal.normalize();

	//light vector
	Vec3Df lightvector;

	//camera vector
	Vec3Df cameravector;

	//"halfway" vector, see wiki page for blinn-phong
	Vec3Df H;

	//get colour of material
	unsigned int triMat = MyMesh.triangleMaterials.at(i);
	Vec3Df col = MyMesh.materials.at(triMat).Kd() * 0.5;

	float blinnPhongResult = 0;
	float diffuseResult = 0;
	Vec3Df reflectionResult = Vec3Df(0, 0, 0);
	Vec3Df refractionResult = Vec3Df(0, 0, 0);



	//iterate through each light source
	for (int j = 0; j < MyLightPositions.size(); j++)
	{
		lightvector = MyLightPositions[j] - vertexPos;
		lightvector.normalize();

		cameravector = MyCameraPosition - vertexPos;
		cameravector.normalize();

		H = lightvector + cameravector;
		H.normalize();

		// Add Blinn-Phong effect
		if (blinnPhongFlag) {
			blinnPhongResult = getBlinnPhong(normal, H, 2);

			result2 += col * blinnPhongResult;
		}

		// Add diffuse effect
		if (diffuseFlag) {
			diffuseResult = getDiffusion(normal, lightvector);

			result2 += col * diffuseResult;
		}


		//check for shadow
		if (shadowFlag) {
			bool shaded = shadow(ray, MyLightPositions[j]);
			if (shaded)
			{
				float shadowdepth = 0.2;

				result2 *= shadowdepth;
			}
		}


		if (reflectionFlag) {
			reflectionResult = getReflection(origin, ray, normal, j, result2);

			result2 += reflectionResult;
		}

		//refraction for transparent objects
		if (refractionFlag) {
			refractionResult = getRefraction(normal, ray, vertexPos);
		}



		result += result2;
	}

	return result;
}



//type for bounding box
class Ray
{
public:
	Ray(Vec3Df orig, Vec3Df dir) : orig(orig),
		dir(dir),
		tmin(float(0)),
		tmax(std::numeric_limits<float>::max())
	{
		invdir[0] = float(1) / dir[0];
		invdir[1] = float(1) / dir[1];
		invdir[2] = float(1) / dir[2];
		sign[0] = (invdir[0] < 0);
		sign[1] = (invdir[1] < 0);
		sign[2] = (invdir[2] < 0);
	}
	Vec3D<float> orig, dir; /// ray orig and dir 
	mutable float tmin, tmax; /// ray min and max distances 
	Vec3D<float> invdir;
	int sign[3];
};

//bounding box
class Box
{
public:
	Box(Vec3D<float> vmin, Vec3D<float> vmax)
	{
		bounds[0] = vmin;
		bounds[1] = vmax;
	}
	Vec3D<float> bounds[2];

	bool intersect(const Ray &r) const {
		float tmin, tmax, tymin, tymax, tzmin, tzmax;
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


//the main function here
//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	//create ray and bounding box
	Ray ray1(origin, dest);
	Box box(bmin, bmax);

	//if the ray doesn't intersect the box, return colour (default black, change to something else for debugging bounding box)
	if (!(box.intersect(ray1)))
	{
		return Vec3Df(0, 0, 0);
	}

	//initialise to far, this float determines closest distance to origin
	float mindistance = 10000;

	//index of triangle in mesh
	int triangleind = -1;

	//used to get vertexPos, intersection point on triangle
	Vec3Df ray;

	Vec3Df intersect;

	const float *p = origin.p;
	const float *d = dest.p;

	//find closest triangle by looping through all
	for (int i = 0; i < MyMesh.triangles.size(); i++){

		float *v0 = MyMesh.vertices[MyMesh.triangles[i].v[0]].p.p;
		float *v1 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p.p;
		float *v2 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p.p;

		float t = rayIntersectsTriangle(p, d, v0, v1, v2, &intersect);

		//t < 0 means no intersect
		if (t >= 0)
		{
			//closest triangle
			if (mindistance > t)
			{
				ray = intersect;
				mindistance = t;
				triangleind = i;
			}
		}
	}


	//if there was an intersection with a triangle
	if (triangleind >= 0)
	{
		return getTriangleColour(triangleind, ray, origin);
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
	glColor3f(0, 1, 1);
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

	std::cout << t << " pressed! The mouse was in location " << x << "," << y << "!" << std::endl;
}