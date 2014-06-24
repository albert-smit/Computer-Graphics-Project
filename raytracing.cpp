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
#include <queue>
#include <functional>


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

Vec3Df getTriangleColour(int, Vec3Df, Vec3Df);

struct node
{
	float distance;
	int triangle;
	int operator<(const node& other)
	{
		return distance < other.distance;
	}
};

void makeBoundingBox()
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
	MyMesh.loadMesh("/Users/stephandumasy/Documents/cgprac/Computer-Graphics-Project/objects/town.obj", true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);

	//make a bounding box around the entire model
	makeBoundingBox();
}

//check if a ray intersects a triangle, 
//intersect returns the intersection point
//output is distance to point
//Moller & Trumbore method
float getTriangleIntersection(Vec3Df p, Vec3Df d, const float *v0, const float *v1, const float *v2, Vec3Df* intersect) {

	float det, invdet, u, v;
	Vec3Df edge1, edge2, pvec, qvec, tvec;

	//edge 1
	edge1[0] = v1[0] - v0[0];
	edge1[1] = v1[1] - v0[1];
	edge1[2] = v1[2] - v0[2];

	//edge 2
	edge2[0] = v2[0] - v0[0];
	edge2[1] = v2[1] - v0[1];
	edge2[2] = v2[2] - v0[2];

	//v0 to origin
	tvec[0] = p[0] - v0[0];
	tvec[1] = p[1] - v0[1];
	tvec[2] = p[2] - v0[2];

	pvec = Vec3Df::crossProduct(d, edge2);

	//determinant
	det = Vec3Df::dotProduct(edge1, pvec);

	//doesn't lie in plane of triangle
	if (det == 0)
	{
		return -1;
	}

	invdet = 1 / det;
	u = invdet * (Vec3Df::dotProduct(tvec, pvec));

	//test u bounds
	if (u < 0.0 || u > 1.0)
	{
		return -1;
	}

	qvec = Vec3Df::crossProduct(tvec, edge1);
	v = invdet * Vec3Df::dotProduct(d, qvec);

	//test v bounds
	if (v < 0.0 || u + v > 1.0)
	{
		return -1;
	}

	//find distance to intersection point
	float t = invdet * Vec3Df::dotProduct(edge2, qvec);

	if (t > 0)
	{
		//point of intersection
		*intersect = (1 - u - v)*Vec3Df(v0[0], v0[1], v0[2]) + u*Vec3Df(v1[0], v1[1], v1[2]) + v*Vec3Df(v2[0], v2[1], v2[2]);
		return t;
	}
	else
	{
		return -1;
	}

}


//check if a point is shaded i.e. no direct light
Vec3Df shadow(Vec3Df origin, Vec3Df dest, Vec3Df currentColour)
{
	float shadowResult = 1;
	Vec3Df intersect;
	Vec3Df matSpecColour = Vec3Df(0, 0, 0);
	Vec3Df matColour = Vec3Df(0, 0, 0);
	Vec3Df sResult = currentColour;
	float matTransparancy = 0;

	Vec3Df p = origin;
	Vec3Df d = dest;

	for (int i = 0; i < MyMesh.triangles.size(); i++){
		float *v0 = MyMesh.vertices[MyMesh.triangles[i].v[0]].p.p;
		float *v1 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p.p;
		float *v2 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p.p;

		float t = getTriangleIntersection(p, d, v0, v1, v2, &intersect);

		//intersect means above 0, rest is to filter out noise
		if (t > 0.0001)
		{
			unsigned int triMat = MyMesh.triangleMaterials.at(i);
			float matTransparancy = MyMesh.materials.at(triMat).Tr();
			float shadowdepth = 0.2;

			// 1 = non-transparent, 0 = fully transparent
			if (matTransparancy > 0 && matTransparancy <= 1) {
				Vec3Df matColour = MyMesh.materials.at(triMat).Kd();
				shadowResult = (1 / matTransparancy) * shadowdepth;
				// sResult = currentColour* shadowResult; Also working doesn't take transparent colour into account!!
				sResult = ((1 - matTransparancy) * currentColour + matTransparancy* matColour) * shadowResult;
			}
		}
	}

	return sResult;
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
	Vec3Df col = MyMesh.materials.at(triMat).Kd();

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
		result2 = shadow(ray, MyLightPositions[currentLightPos], result2);
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

		Vec3Df p = selectedPos;
		Vec3Df d = r;

		//find closest triangle by looping through all
		for (int o = 0; o < MyMesh.triangles.size(); o++){

			float *v0 = MyMesh.vertices[MyMesh.triangles[o].v[0]].p.p;
			float *v1 = MyMesh.vertices[MyMesh.triangles[o].v[1]].p.p;
			float *v2 = MyMesh.vertices[MyMesh.triangles[o].v[2]].p.p;

			float t = getTriangleIntersection(p, d, v0, v1, v2, &intersect);

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
			Vec3Df colr = MyMesh.materials.at(triMatr).Ks();
			reflectionColour *= colr;
		}

	}//end of the reflection part
	return reflectionColour;
}

Vec3Df getRefraction2(Vec3Df cameraPos, Vec3Df selectedPos, Vec3Df normal, int currentLightPos, Vec3Df currentResult) {
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

		Vec3Df p = selectedPos;
		Vec3Df d = r;

		//find closest triangle by looping through all
		for (int o = 0; o < MyMesh.triangles.size(); o++){

			float *v0 = MyMesh.vertices[MyMesh.triangles[o].v[0]].p.p;
			float *v1 = MyMesh.vertices[MyMesh.triangles[o].v[1]].p.p;
			float *v2 = MyMesh.vertices[MyMesh.triangles[o].v[2]].p.p;

			float t = getTriangleIntersection(p, d, v0, v1, v2, &intersect);

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
			Vec3Df colr = MyMesh.materials.at(triMatr).Ks();
			reflectionColour *= colr;
		}

	}//end of the reflection part
	return reflectionColour;
}


//called when ray from origin intersects a triangle
//function that gets the colour
Vec3Df getTriangleColour(int i, Vec3Df ray, Vec3Df origin, Vec3Df dest, std::priority_queue<node*> distances)
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
	Vec3Df col = MyMesh.materials.at(triMat).Kd();

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

		//refraction for transparent objects
		if (refractionFlag) {
			//refractionResult = getRefraction(dest, ray);
			unsigned int triMat = MyMesh.triangleMaterials.at(i);
			float matTransparancy = MyMesh.materials.at(triMat).Tr();
			if (distances.size() > 1 && matTransparancy < 1) {
				distances.pop();
				node* secondTriangle = distances.top();
				int triangleId = secondTriangle->triangle;

				unsigned int triMat = MyMesh.triangleMaterials.at(triangleId);
				Vec3Df refractionResult = matTransparancy * MyMesh.materials.at(triMat).Kd();
				Vec3Df colr = MyMesh.materials.at(triMat).Ks();
				result2 = (1 - matTransparancy) * result2 + refractionResult*colr;
			}
		}

		//check for shadow
		if (shadowFlag) {
			result2 = shadow(ray, MyLightPositions[j], result2);
		}


		if (reflectionFlag) {
			reflectionResult = getReflection(origin, ray, normal, j, result2);

			result2 += reflectionResult;
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

//bounding box intersection test
bool rayBoxIntersect(const Vec3Df & origin, const Vec3Df & dest)
{
	float min, max;
	min = bmin[0];
	max = bmax[0];

	Vec3Df invdest, bounds[2];
	int xsign, ysign, zsign;
	float xmin, xmax, ymin, ymax, zmin, zmax;

	bounds[0] = bmin;
	bounds[1] = bmax;

	invdest[0] = 1 / dest[0];
	invdest[1] = 1 / dest[1];
	invdest[2] = 1 / dest[2];

	xsign = (invdest[0] < 0);
	ysign = (invdest[1] < 0);
	zsign = (invdest[2] < 0);

	xmin = (bounds[xsign][0] - origin[0]) * invdest[0];
	xmax = (bounds[1 - xsign][0] - origin[0]) * invdest[0];
	ymin = (bounds[ysign][1] - origin[1]) * invdest[1];
	ymax = (bounds[1 - ysign][1] - origin[1]) * invdest[1];
	zmin = (bounds[zsign][2] - origin[2]) * invdest[2];
	zmax = (bounds[1 - zsign][2] - origin[2]) * invdest[2];

	if ((xmin > ymax) || (ymin > xmax))
	{
		return false;
	}

	if (ymin > xmin)
	{
		xmin = ymin;
	}

	if (ymax < xmax)
	{
		xmax = ymax;
	}

	if ((xmin > zmax) || (zmin > xmax))
	{
		return false;
	}

	if (zmin > xmin)
	{
		xmin = zmin;
	}

	if (zmax < xmax)
	{
		xmax = zmax;
	}

	if (xmin > min)
	{
		min = xmin;
	}

	if (xmax < max)
	{
		max = xmax;
	}

	return true;
}



Vec3Df performSubRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{

	//if the ray doesn't intersect the box, return colour (default black, change to something else for debugging bounding box)
	if (!rayBoxIntersect(origin, dest))
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

	Vec3Df p = origin;
	Vec3Df d = dest;

	//find closest triangle by looping through all
	std::priority_queue<node*> distances;
	//std::priority_queue<float, std::vector<float>, std::greater<float>> distances;

	for (int i = 0; i < MyMesh.triangles.size(); i++){

		float *v0 = MyMesh.vertices[MyMesh.triangles[i].v[0]].p.p;
		float *v1 = MyMesh.vertices[MyMesh.triangles[i].v[1]].p.p;
		float *v2 = MyMesh.vertices[MyMesh.triangles[i].v[2]].p.p;

		float distance = getTriangleIntersection(p, d, v0, v1, v2, &intersect);

		//t < 0 means no intersect
		if (distance >= 0)
		{
			node* n = new node;
			n->distance = distance;
			n->triangle = i;
			distances.push(n);

			//closest triangle
			if (mindistance > distance)
			{
				ray = intersect;
				mindistance = distance;
				triangleind = i;
			}
		}
	}


	//if there was an intersection with a triangle
	if (triangleind >= 0)
	{
		return getTriangleColour(triangleind, ray, origin, dest, distances);
	}

	return Vec3Df(0, 0, 0);

}

//the main function here
//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	//Distance between the pixels
	float sampleDistance = 0.0041;

	//Determine the direction of the ray
	Vec3Df rayvector = origin - dest;

	//Get the orthogonal of the ray and normalize it with the pixel distance
	Vec3Df vX = Vec3Df(-rayvector[1], rayvector[0], rayvector[2]);
	vX.normalize();
	vX *= sampleDistance;

	Vec3Df vY = Vec3Df(-rayvector[2], rayvector[1], rayvector[0]);
	vY.normalize();
	vY *= sampleDistance;

	//Get the color of every new ray and divide by the number of rays.
	Vec3Df resultTotal = performSubRayTracing(origin, dest) + performSubRayTracing(origin + vX, dest + vX) + performSubRayTracing(origin - vX, dest - vX) + performSubRayTracing(origin + vY, dest + vY) + performSubRayTracing(origin - vY, dest - vY);

	return resultTotal / 5;


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