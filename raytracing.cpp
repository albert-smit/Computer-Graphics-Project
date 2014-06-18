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
    MyMesh.loadMesh("cube.obj", true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);
}

//return the depth of the triangle
float triangleCollisionDistance(const Vec3Df & origin, const Vec3Df & dest, const Triangle & triangle)
{
    
    //printf("normal: %f", normal);
    
    Vec3Df v0 = MyMesh.vertices[triangle.v[0]].p;
    Vec3Df v1 = MyMesh.vertices[triangle.v[1]].p;
    Vec3Df v2 = MyMesh.vertices[triangle.v[2]].p;
    Vec3Df v01 = v1 - v0;
    Vec3Df v02 = v2 - v0;
    Vec3Df normal = Vec3Df::crossProduct(v01,v02);
    
    float D = Vec3Df::dotProduct(normal, v0);
    //std::cout" normal: "<<normal<<"!"<<std::endl;
    return -(Vec3Df::dotProduct(normal, origin) + D) / Vec3Df::dotProduct(normal, dest);
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
    float distance = (origin - dest).getLength();
    Triangle triangle;
    
    for(int i = 0; i < MyMesh.triangles.size(); i ++){
        float collision = triangleCollisionDistance(origin, dest, MyMesh.triangles[i]);
        if (collision < distance){
            distance = collision;
            triangle = MyMesh.triangles[i];
        }
    }
    
    return Vec3Df(triangle.v[0],triangle.v[1],triangle.v[2]);

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
