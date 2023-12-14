//Ioannis Kalaitzoglou AM 2982

// Include standard headers
#include <stdio.h>
#include <stdlib.h>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

#include <shader.hpp>
#include <iostream>
#include <string>
using namespace std;

#include <stb_image.h>
#include <ctime>
#include <vector>

//################ GLOBAL VARIABLES #########################
//General globals
const float PI = acos(-1);

//Small object global variables
int numberOfSmallObjects = 0;
//Vector that keeps all smallObjects buffers
std::vector < GLuint > smallObjects;
std::vector < unsigned int > indices_size;
//Vector that keeps small object direction 
std::vector < float > direction;
GLuint buffers_array[3];
int index_count = 0;
int smallcounter = 0;
//Vector that keeps all objects dimensions
std::vector < float > objectDim;
//Vector that keeps all objects velocities
std::vector < float > velocity;

//Big sphere global variables
int bigSphereIndicesCount;
std::vector < float >	BigSpherevertices;
std::vector < unsigned int > BigSphereindices;
std::vector < float > textureCoords;
std::vector < float > SpherePosition;
int BigSphereRadius = 15;
float BigSphereVelocity = 0.5f;
int textureOnOff = 0;

//Function to handle small object collision with wall
void smallObjectCollisionWithWall(int i, int index_count) {
	float dimension = objectDim[index_count];
	if (direction[i] + dimension > 100) {
		direction[i] = 100 - dimension;
		velocity[i] = -velocity[i];
	}
	else if (direction[i] - dimension < 0 - dimension ) {
		direction[i] = 0 + dimension;
		velocity[i] = -velocity[i];
	}
	else if (direction[i + 1] + dimension > 100) {
		direction[i + 1] = 100 - dimension;
		velocity[i + 1] = -velocity[i + 1];
	}
	else if (direction[i + 1] - dimension < 0 - dimension) {
		direction[i + 1] = 0 + dimension;
		velocity[i + 1] = -velocity[i + 1];
	}
	else if (direction[i + 2] + dimension > 100) {
		direction[i + 2] = 100 - dimension;
		velocity[i + 2] = -velocity[i + 2];
	}
	else if (direction[i + 2] - dimension < 0 - dimension) {
		direction[i + 2] = 0 + dimension;
		velocity[i + 2] = -velocity[i + 2];
	}
}
//Function to handle small object collision with another small object
void smallObjectCollisionWithSmallObject(int i, int index_count) {
	float dimension = objectDim[index_count];
	int SmallObjectCounter = 0;
	int ind_count = 0;
	while (SmallObjectCounter < numberOfSmallObjects) {
		if (SmallObjectCounter != i) {
			float smallObjectDistX = abs(direction[i] - direction[SmallObjectCounter]);
			float smallObjectDistY = abs(direction[i + 1] - direction[SmallObjectCounter + 1]);
			float smallObjectDistZ = abs(direction[i + 2] - direction[SmallObjectCounter + 2]);
			float distance = sqrtf(smallObjectDistX * smallObjectDistX + smallObjectDistY * smallObjectDistY + smallObjectDistZ * smallObjectDistZ);
			float sumOfRad = objectDim[ind_count] + dimension;
			if (distance <= sumOfRad) {
				velocity[i] = -velocity[i];
				velocity[i + 1] = -velocity[i + 1];
				velocity[i + 2] = -velocity[i + 2];
				velocity[SmallObjectCounter] = -velocity[SmallObjectCounter];
				velocity[SmallObjectCounter + 1] = -velocity[SmallObjectCounter + 1];
				velocity[SmallObjectCounter + 2] = -velocity[SmallObjectCounter + 2];
			}
		}
		SmallObjectCounter += 3;
		ind_count++;
	}
}
//Function to handle big sphere collision with wall
void BigSphereCollisionWithWall() {
	
	if (SpherePosition[0] + BigSphereRadius > 100) {
		SpherePosition[0] = 100 - BigSphereRadius;	
	}
	else if (SpherePosition[0] - BigSphereRadius < 0 ) {
		SpherePosition[0] = 0 + BigSphereRadius;
	}
	else if (SpherePosition[1] + BigSphereRadius > 100) {
		SpherePosition[1] = 100 - BigSphereRadius;
	}
	else if (SpherePosition[1] - BigSphereRadius < 0 ) {
		SpherePosition[1] = 0 + BigSphereRadius;
	}
	else if (SpherePosition[2] + BigSphereRadius > 100) {
		SpherePosition[2] = 100 - BigSphereRadius;
	}
	else if (SpherePosition[2] - BigSphereRadius < 0) {
		SpherePosition[2] = 0 + BigSphereRadius;
	}
}
//Function to handle small object collision with sphere
void smallObjectCollisionWithSphere(int i, int index_count) {
	float dimension = objectDim[index_count];
	float sphereObjectDistX = abs(SpherePosition[0] - direction[i]);
	float sphereObjectDistY = abs(SpherePosition[1] - direction[i + 1]);
	float sphereObjectDistZ = abs(SpherePosition[2] - direction[i + 2]);
	float distance = sqrtf(sphereObjectDistX * sphereObjectDistX + sphereObjectDistY * sphereObjectDistY + sphereObjectDistZ * sphereObjectDistZ);
	float sumOfRad =  BigSphereRadius + dimension;
	
	if (distance <= sumOfRad) {
		velocity[i] = -velocity[i];
		velocity[i + 1] = -velocity[i + 1];
		velocity[i + 2] = -velocity[i + 2];
	}
	
}
//Function to generate random starting small object direction between [0.1, 0.9]
glm::vec3 generateObjectDirection() {
	float vx, vy, vz;
	
	vx = (((float)rand()) / (float)RAND_MAX *0.8) - 0.1;
	vy = (((float)rand()) / (float)RAND_MAX * 0.8) - 0.1;
	vz = (((float)rand()) / (float)RAND_MAX * 0.8) - 0.1;

	return glm::vec3(vx, vy, vz);
}
//Function to generate random starting small object velocity
glm::vec3 generateObjectVelocity() {
	float vx, vy, vz;
	
	vx = (((float)rand()) / (float)RAND_MAX *0.9) - 0.1;
	vy = (((float)rand()) / (float)RAND_MAX * 0.9) - 0.1;
	vz = (((float)rand()) / (float)RAND_MAX * 0.9) - 0.1;

	return glm::vec3(vx, vy, vz);
}
//Function to create big sphere
void createBigSphere() {
	const float PI = acos(-1);
	int radius = BigSphereRadius;
	float x, y, z, xy;   // vertex position
	float sectorCount = 36.0f;
	float stackCount = 18.0f;
	float sectorStep = 2 * PI / sectorCount;
	float stackStep = PI / stackCount;
	float sectorAngle, stackAngle;
	float s, t;
	
	for (int i = 0; i <= stackCount; ++i)
	{
		stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
		xy = radius * cosf(stackAngle);             // r * cos(u)
		z = radius * sinf(stackAngle);              // r * sin(u)

		// add (sectorCount+1) vertices per stack
		// the first and last vertices have same position but different tex coords
		for (int j = 0; j <= sectorCount; ++j)
		{
			sectorAngle = j * sectorStep;           // starting from 0 to 2pi

			// vertex position (x, y, z)
			x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
			y = xy * sinf(sectorAngle);

			BigSpherevertices.push_back(x);
			BigSpherevertices.push_back(y);
			BigSpherevertices.push_back(z);

			s = (float)j / sectorCount;
			t = (float)i / stackCount;
			textureCoords.push_back(s);
			textureCoords.push_back(t);

		}
	}

	int k1, k2;
	for (int i = 0; i < stackCount; ++i)
	{
		k1 = i * (sectorCount + 1);     // beginning of current stack
		k2 = k1 + sectorCount + 1;      // beginning of next stack

		for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
		{
			// 2 triangles per sector excluding first and last stacks
			// k1 => k2 => k1+1
			if (i != 0)
			{
				BigSphereindices.push_back(k1);
				BigSphereindices.push_back(k2);
				BigSphereindices.push_back(k1 + 1);

			}

			// k1+1 => k2 => k2+1
			if (i != (stackCount - 1))
			{
				BigSphereindices.push_back(k1 + 1);
				BigSphereindices.push_back(k2);
				BigSphereindices.push_back(k2 + 1);

			}
		}
	}

	
	bigSphereIndicesCount = (unsigned int)BigSphereindices.size();
	SpherePosition.push_back(50.0f);
	SpherePosition.push_back(50.0f);
	SpherePosition.push_back(50.0f);
}
//Funstion to create a random small object (cube,sphere,cylinder)
void createRandomObject() {
	int dimension;
	int object_type;
	std::vector < float > colorbuffer;
	std::vector < float >	vertices;
	std::vector < unsigned int > indices;
	// random number generator. Generates numbers between 1 and 3 for defining the object type. (0=sphere, 1=cube, 2=cylinder)
	object_type = rand() % 3;
	// random number generator. Generated number between 1 and 10 for defining the dimensions of the object.
	dimension = rand() % 10 + 1;	
	//random R G B values for defining the color of the object
	float redValue = (float)rand() / (float)RAND_MAX;
	float greenValue = (float)rand() / (float)RAND_MAX;
	float blueValue = (float)rand() / (float)RAND_MAX;

	// Object is sphere
	if (object_type == 0) {
		float radius = float(dimension) / 2.0f;
		float x, y, z, xy;   // vertex position
		float nx, ny, nz;
		float lengthInv = 1.0f / radius;
		float sectorCount = 36.0f;
		float stackCount = 18.0f;
		float sectorStep = 2 * PI / sectorCount;
		float stackStep = PI / stackCount;
		float sectorAngle, stackAngle;

		for (int i = 0; i <= stackCount; ++i)
		{
			stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
			xy = radius * cosf(stackAngle);             // r * cos(u)
			z = radius * sinf(stackAngle);              // r * sin(u)

			// add (sectorCount+1) vertices per stack
			// the first and last vertices have same position but different tex coords
			for (int j = 0; j <= sectorCount; ++j)
			{
				sectorAngle = j * sectorStep;           // starting from 0 to 2pi

				// vertex position (x, y, z)
				x = xy * cosf(sectorAngle);             // r * cos(u) * cos(v)
				y = xy * sinf(sectorAngle);

				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);

			}
		}

		int k1, k2;
		for (int i = 0; i < stackCount; ++i)
		{
			k1 = i * (sectorCount + 1);     // beginning of current stack
			k2 = k1 + sectorCount + 1;      // beginning of next stack

			for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
			{
				// 2 triangles per sector excluding first and last stacks
				// k1 => k2 => k1+1
				if (i != 0)
				{
					indices.push_back(k1);
					indices.push_back(k2);
					indices.push_back(k1 + 1);

				}

				// k1+1 => k2 => k2+1
				if (i != (stackCount - 1))
				{
					indices.push_back(k1 + 1);
					indices.push_back(k2);
					indices.push_back(k2 + 1);

				}
			}
		}
		// Color buffer small sphere 
		for (int i = 0; i <= stackCount; ++i)
		{

			for (int j = 0; j <= sectorCount; ++j)
			{
				colorbuffer.push_back(redValue);
				colorbuffer.push_back(greenValue);
				colorbuffer.push_back(blueValue);
			}
		}
		//Small sphere buffers initialize
		GLuint vertexbufferSphere;
		glGenBuffers(1, &vertexbufferSphere);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferSphere);
		glBufferData(GL_ARRAY_BUFFER, (float)vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

		GLuint elementbuffer;
		glGenBuffers(1, &elementbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (unsigned int)indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

		GLuint spherecolorbuffer;
		glGenBuffers(1, &spherecolorbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spherecolorbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (float)colorbuffer.size() * sizeof(float), colorbuffer.data(), GL_STATIC_DRAW);

		buffers_array[0] = vertexbufferSphere;
		buffers_array[1] = elementbuffer;
		buffers_array[2] = spherecolorbuffer;
		smallObjects.push_back(buffers_array[0]);
		smallObjects.push_back(buffers_array[1]);
		smallObjects.push_back(buffers_array[2]);
		indices_size.push_back((unsigned int)indices.size());
		objectDim.push_back(radius);
	}
	else if (object_type == 1) {
		float cube_side = (float)dimension /2.0f;

		vertices.push_back(-cube_side); 
		vertices.push_back(-cube_side);
		vertices.push_back(cube_side);

		vertices.push_back(cube_side);
		vertices.push_back(-cube_side);
		vertices.push_back(cube_side);

		vertices.push_back(cube_side);
		vertices.push_back(cube_side);
		vertices.push_back(cube_side);

		vertices.push_back(-cube_side);
		vertices.push_back(cube_side);
		vertices.push_back(cube_side);
	
		vertices.push_back(-cube_side);
		vertices.push_back(-cube_side);
		vertices.push_back(-cube_side);

		vertices.push_back(cube_side);
		vertices.push_back(-cube_side);
		vertices.push_back(-cube_side);
		
		vertices.push_back(cube_side);
		vertices.push_back(cube_side);
		vertices.push_back(-cube_side);
		
		vertices.push_back(-cube_side);
		vertices.push_back(cube_side);
		vertices.push_back(-cube_side);

		for (int i = 0; i < 8; i++) {
			colorbuffer.push_back(redValue);
			colorbuffer.push_back(greenValue);
			colorbuffer.push_back(blueValue);
		}


		//front
		indices.push_back(0);
		indices.push_back(1);
		indices.push_back(2);

		indices.push_back(2);
		indices.push_back(3);
		indices.push_back(0);
		//right
		indices.push_back(1);
		indices.push_back(5);
		indices.push_back(6);

		indices.push_back(6);
		indices.push_back(2);
		indices.push_back(1);

		//back
		indices.push_back(7);
		indices.push_back(6);
		indices.push_back(5);

		indices.push_back(5);
		indices.push_back(4);
		indices.push_back(7);

		//left
		indices.push_back(4);
		indices.push_back(0);
		indices.push_back(3);
			
		indices.push_back(3);
		indices.push_back(7);
		indices.push_back(4);
		
		//bottom
		indices.push_back(4);
		indices.push_back(5);
		indices.push_back(1);

		indices.push_back(1);
		indices.push_back(0);
		indices.push_back(4);

		//top
		indices.push_back(3);
		indices.push_back(2);
		indices.push_back(6);

		indices.push_back(6);
		indices.push_back(7);
		indices.push_back(3);

		
		GLuint vertexbufferCube;
		glGenBuffers(1, &vertexbufferCube);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferCube);
		glBufferData(GL_ARRAY_BUFFER, (float)vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

		GLuint elementbuffer;
		glGenBuffers(1, &elementbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (unsigned int)indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

		GLuint cubecolorbuffer;
		glGenBuffers(1, &cubecolorbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubecolorbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (float)colorbuffer.size() * sizeof(float), colorbuffer.data(), GL_STATIC_DRAW);

		buffers_array[0] = vertexbufferCube;
		buffers_array[1] = elementbuffer;
		buffers_array[2] = cubecolorbuffer;
		smallObjects.push_back(buffers_array[0]);
		smallObjects.push_back(buffers_array[1]);
		smallObjects.push_back(buffers_array[2]);
		indices_size.push_back((unsigned int)indices.size());
		objectDim.push_back(cube_side);
		
	}
	//cyclinder 
	else if (object_type == 2) {
	float sectorCount = 36.0f;
	float stackCount = 18.0f;
	float sectorStep = 2 * PI / sectorCount;
	float sectorAngle;  
	float height = (float)dimension;
	float radius = (float)dimension;
	std::vector<float> circleVertices;
		for (int i = 0; i <= sectorCount; ++i)
		{
			sectorAngle = i * sectorStep;
			circleVertices.push_back(cos(sectorAngle)); // x
			circleVertices.push_back(sin(sectorAngle)); // y
			circleVertices.push_back(0);                // z
		}
	// put side vertices to arrays
		for (int i = 0; i < 2; ++i)
		{
			float h = -height / 2.0f + i * height;           // z value; -h/2 to h/2
			float t = 1.0f - i;                              // vertical tex coord; 1 to 0

			for (int j = 0, k = 0; j <= sectorCount; ++j, k += 3)
			{
				float ux = circleVertices[k];
				float uy = circleVertices[k + 1];
				float uz = circleVertices[k + 2];
				// position vector
				vertices.push_back(ux * radius);             // vx
				vertices.push_back(uy * radius);             // vy
				vertices.push_back(h);                       // vz
			}
		}

	// the starting index for the base/top surface
   // it is used for generating indices later
	int baseCenterIndex = (int)vertices.size() / 3;
	int topCenterIndex = baseCenterIndex + sectorCount + 1; // include center vertex

	// put base and top vertices to arrays
		for (int i = 0; i < 2; ++i)
		{
			float h = -height / 2.0f + i * height;           // z value; -h/2 to h/2
			float nz = -1 + i * 2;                           // z value of normal; -1 to 1

			// center point
			vertices.push_back(0); 
			vertices.push_back(0);
			vertices.push_back(h);		

			for (int j = 0, k = 0; j < sectorCount; ++j, k += 3)
			{
				float ux = circleVertices[k];
				float uy = circleVertices[k + 1];
				// position vector
				vertices.push_back(ux * radius);             // vx
				vertices.push_back(uy * radius);             // vy
				vertices.push_back(h);                       // vz
			}
		}

		int k1 = 0;                         // 1st vertex index at base
		int k2 = sectorCount + 1;           // 1st vertex index at top

		// indices for the side surface
		for (int i = 0; i < sectorCount; ++i, ++k1, ++k2)
		{
			// 2 triangles per sector
			// k1 => k1+1 => k2
			indices.push_back(k1);
			indices.push_back(k1 + 1);
			indices.push_back(k2);

			// k2 => k1+1 => k2+1
			indices.push_back(k2);
			indices.push_back(k1 + 1);
			indices.push_back(k2 + 1);
		}

		// indices for the base surface
		for (int i = 0, k = baseCenterIndex + 1; i < sectorCount; ++i, ++k)
		{
			if (i < sectorCount - 1)
			{
				indices.push_back(baseCenterIndex);
				indices.push_back(k + 1);
				indices.push_back(k);
			}
			else // last triangle
			{
				indices.push_back(baseCenterIndex);
				indices.push_back(baseCenterIndex + 1);
				indices.push_back(k);
			}
		}

		// indices for the top surface
		for (int i = 0, k = topCenterIndex + 1; i < sectorCount; ++i, ++k)
		{
			if (i < sectorCount - 1)
			{
				indices.push_back(topCenterIndex);
				indices.push_back(k);
				indices.push_back(k + 1);
			}
			else // last triangle
			{
				indices.push_back(topCenterIndex);
				indices.push_back(k);
				indices.push_back(topCenterIndex + 1);
			}
		}
		// Color buffer small sphere 
		for (int i = 0; i <= stackCount; ++i)
		{

			for (int j = 0; j <= sectorCount; ++j)
			{
				colorbuffer.push_back(redValue);
				colorbuffer.push_back(greenValue);
				colorbuffer.push_back(blueValue);
			}
		}
		GLuint vertexbufferCylinder;
		glGenBuffers(1, &vertexbufferCylinder);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferCylinder);
		glBufferData(GL_ARRAY_BUFFER, (float)vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

		GLuint elementbuffer;
		glGenBuffers(1, &elementbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (unsigned int)indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

		GLuint cylindercolorbuffer;
		glGenBuffers(1, &cylindercolorbuffer);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cylindercolorbuffer);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, (float)colorbuffer.size() * sizeof(float), colorbuffer.data(), GL_STATIC_DRAW);

		buffers_array[0] = vertexbufferCylinder;
		buffers_array[1] = elementbuffer;
		buffers_array[2] = cylindercolorbuffer;
		smallObjects.push_back(buffers_array[0]);
		smallObjects.push_back(buffers_array[1]);
		smallObjects.push_back(buffers_array[2]);
		indices_size.push_back((unsigned int)indices.size());
		objectDim.push_back(radius);
	}
	//Generate random direction and velocity for the small object and keep track of our small object number
	glm::vec3 directionsFromRng;
	directionsFromRng = generateObjectDirection();
	direction.push_back(directionsFromRng.x);
	direction.push_back(directionsFromRng.y);
	direction.push_back(directionsFromRng.z);
	glm::vec3 velocityRng;
	velocityRng = generateObjectVelocity();
	velocity.push_back(velocityRng.x);
	velocity.push_back(velocityRng.y);
	velocity.push_back(velocityRng.z);
	numberOfSmallObjects+=3;
}

//Draw small object enables small object buffers and draws object
void drawSmallObject(GLuint vertex_buffer,GLuint index_buffer,GLuint color_buffer, unsigned int indices_size) {
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
	glVertexAttribPointer(0,3,GL_FLOAT, GL_FALSE, 0, (void*)0	);

	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, color_buffer);
	glVertexAttribPointer(1,3,GL_FLOAT,	GL_FALSE,0,(void*)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
	
	glDrawElements(GL_TRIANGLES, indices_size, GL_UNSIGNED_INT, (void*)0);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
//Draw big sphere by enabling big sphere's buffer and  drawing our sphere
void drawBigSphereObject(GLuint vertex_buffer, GLuint index_buffer, GLuint texture_buffer, unsigned int indices_size) {
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);
	
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, texture_buffer);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

	glDrawElements(GL_TRIANGLES, indices_size, GL_UNSIGNED_INT, (void*)0);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
//Flag function used enable/disable big sphere texture
void EnableDisableTexture() {
	if (textureOnOff == 0) {
		textureOnOff = 1;
	}
	else {
		textureOnOff = 0;
	}
}
//Call back function that is called every time user presses either "SPACE" or "T"
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		createRandomObject();		
	}	
	if (key == GLFW_KEY_T && action == GLFW_PRESS) {
		EnableDisableTexture();
	}
}

//####################### MAIN ############################
int main(void)
{
	// Seed time to random number generator
	srand(time(NULL));
	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	}

	//Specifications for our window
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window with size 600x600 and create its OpenGL context
	window = glfwCreateWindow(600, 600, u8"Συγκρουόμενα", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window\n");
		glfwTerminate();
		return -1;
	}	
	glfwMakeContextCurrent(window);
	
	//Key callback  to handle user pressing "SPACE" for creating random objects and "T" to enable/disable texture 
	glfwSetKeyCallback(window, keyCallback);
	
	// Initialize GLEW
	glewExperimental = true; 
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		glfwTerminate();
		return -1;
	}
	// Black backround 
	glClearColor(0.0f,0.0f, 0.0f, 0.0f);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders 
	//This shader is used for all object instances except big sphere object
	GLuint programID = LoadShaders("TransformVertexShader.vertexshader", "ColorFragmentShader.fragmentshader");
	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");

	//Create GLSL program from the shader for our big sphere object and handle MVP, texture, textureFlag
	GLuint textureShaderID = LoadShaders("TexturedSphere.vertexshader", "TextureSpherefrag.fragmentshader");
	GLuint TextureID = glGetUniformLocation(textureShaderID, "texture1");
	GLuint Matrix2ID = glGetUniformLocation(textureShaderID, "MVP");
	GLuint OpacityID = glGetUniformLocation(textureShaderID, "textureFlag");
	
	
	// SC (Scene Cube) vertex array starting at (0,0,0) to (100,100,100)
	static const GLfloat scene_cube_vertex_buffer_data[] = {
		//Face 1
		0.0f , 0.0f , 0.0f,
		100.0f, 100.0f, 0.0f,
		100.0f, 0.0f, 0.0f,
		

		0.0f, 0.0f, 0.0f,
		0.0f, 100.0f, 0.0f,
		100.0f, 100.0f, 0.0f,
		//Face 2
		0.0f, 100.0f, 0.0f,		
		0.0f, 100.0f, 100.0f,
		100.0f, 100.0f, 0.0f,

		0.0f, 100.0f, 100.0f,
		100.0f, 100.0f, 100.0f,
		100.0f, 100.0f, 0.0f,

		//Face 3
		0.0f, 100.0f, 0.0f,
		0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 100.0f,

		0.0f, 100.0f, 0.0f,
		0.0f, 0.0f, 100.0f,
		0.0f, 100.0f, 100.0f,
		

		//Face 4
		0.0f, 100.0f, 100.0f,
		100.0f, 0.0f, 100.0f,
		100.0f, 100.0f, 100.0f,
		

		0.0f, 100.0f, 100.0f,
		0.0f, 0.0f, 100.0f,
		100.0f, 0.0f, 100.0f,
		

		//Face 5
		100.0f, 100.0f, 0.0f,
		100.0f, 0.0f, 100.0f,
		100.0f, 0.0f, 0.0f,
		

		100.0f, 100.0f, 0.0f,
		100.0f, 100.0f, 100.0f,
		100.0f, 0.0f, 100.0f,

		//Face 6
		0.0f, 0.0f, 0.0f,
		100.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 100.0f,

		100.0f, 0.0f, 0.0f,
		100.0f, 0.0f, 100.0f,
		0.0f, 0.0f, 100.0f,
	};

	//Scene cube color array initialized with random color and semi-transparent opacity
	static GLfloat scene_cube_color_buffer_data[12*3*4];	
	float red = (float)rand() / (float)RAND_MAX;
	float green = (float)rand() / (float)RAND_MAX;
	float blue = (float)rand() / (float)RAND_MAX;
	float opacity = 0.4f;
	for (int v = 0; v < 12 * 4  ; v++) {
		scene_cube_color_buffer_data[4 * v] = red;
		scene_cube_color_buffer_data[4 * v + 1] =  green;
		scene_cube_color_buffer_data[4 * v + 2] = blue;
		scene_cube_color_buffer_data[4 * v + 3] = opacity;
	}

	//Scene cube wireframe array. Lines created on cube edges
	static const GLfloat scene_cube_wireframe_vertex_buffer_data[] = {		
		0.0f , 0.0f , 0.0f,
		0.0f, 100.0f, 0.0f,

		0.0f , 0.0f , 0.0f,
		100.0f, 0.0f, 0.0f,

		0.0f , 0.0f , 0.0f,
		0.0f, 0.0f, 100.0f,

		0.0f, 100.0f, 0.0f,
		0.0f, 100.0f, 100.0f,

		0.0f, 100.0f, 0.0f,
		100.0f, 100.0f, 0.0f,
		
		0.0f, 100.0f, 100.0f,
		0.0f, 0.0f, 100.0f,

		0.0f, 0.0f, 100.0f,
		100.0f, 0.0f, 100.0f,

		100.0f, 0.0f, 100.0f,
		100.0f, 100.0f, 100.0f,

		100.0f, 100.0f, 100.0f,
		0.0f, 100.0f, 100.0f,

		100.0f, 100.0f, 100.0f,
		100.0f, 100.0f, 0.0f,
		
		100.0f, 100.0f, 0.0f,
		100.0f, 0.0f, 0.0f,

		100.0f, 0.0f, 0.0f,
		100.0f, 0.0f, 100.0f,
		
	};	

	//Scene cube wireframe color array. Initialized to have white color.
	static GLfloat scene_cube_wireframe_color_buffer_data[12 * 2 * 3];
	for (int v = 0; v <12 * 2; v++){
		
		scene_cube_wireframe_color_buffer_data[3 * v] = 1.0f;
		scene_cube_wireframe_color_buffer_data[3 * v + 1] = 1.0f;
		scene_cube_wireframe_color_buffer_data[3 * v + 2] = 1.0f;
	}

	//Z-buffer Depth test to display front and back fragment correctly
	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	//Createbigsphere() function call to create the big (red) sphere.
	createBigSphere();
	
	
	// Load and create a texture 
	// -------------------------
	unsigned int texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
	// set the texture wrapping parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	// set texture filtering parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	// load image, create texture and generate mipmaps
	int width, height, nrChannels;
	// The FileSystem::getPath(...) is part of the GitHub repository so we can find files on any IDE/platform; replace it with your own image path.
	unsigned char* data = stbi_load("texture-sphere.jpg", &width, &height, &nrChannels, 0);
	if (data)
	{
		//PixelStorei added to the code to handle (progressive JPEG) format error for our texture image. 
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
		glGenerateMipmap(GL_TEXTURE_2D);
	}
	else
	{

		std::cout << "Failed to load texture" << std::endl;
	}
	stbi_image_free(data);
	
	//############## Buffers ####################
	//Scene cube vertex buffer initialization
	GLuint sc_vertexbuffer;
	glGenBuffers(1, &sc_vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, sc_vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(scene_cube_vertex_buffer_data), scene_cube_vertex_buffer_data, GL_STATIC_DRAW);
	//Scene cube color buffer initialization
	GLuint sc_colorbuffer;
	glGenBuffers(1, &sc_colorbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, sc_colorbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(scene_cube_color_buffer_data), scene_cube_color_buffer_data, GL_STATIC_DRAW);

	//Scene cube wireframe vertex buffer initialization
	GLuint wireframe_vertexbuffer;
	glGenBuffers(1, &wireframe_vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, wireframe_vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(scene_cube_wireframe_vertex_buffer_data), scene_cube_wireframe_vertex_buffer_data, GL_STATIC_DRAW);
	//Scene cube wireframe color buffer initialization
	GLuint wireframe_colorbuffer;
	glGenBuffers(1, &wireframe_colorbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, wireframe_colorbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(scene_cube_wireframe_color_buffer_data), scene_cube_wireframe_color_buffer_data, GL_STATIC_DRAW);

	//Big sphere vertex buffer initialization
	GLuint vertexbufferBigSphere;
	glGenBuffers(1, &vertexbufferBigSphere);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbufferBigSphere);
	glBufferData(GL_ARRAY_BUFFER, (float)BigSpherevertices.size() * sizeof(float), BigSpherevertices.data(), GL_STATIC_DRAW);
	//Big sphere index buffer initialization
	GLuint elementbufferBigSphere;
	glGenBuffers(1, &elementbufferBigSphere);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbufferBigSphere);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, bigSphereIndicesCount * sizeof(unsigned int), BigSphereindices.data(), GL_STATIC_DRAW);
	//Big sphere texture buffer initialization
	GLuint texturebufferBigSphere;
	glGenBuffers(1, &texturebufferBigSphere);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, texturebufferBigSphere);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, (float)textureCoords.size() * sizeof(float), textureCoords.data(), GL_STATIC_DRAW);
	
	// Starting camera's position in World Space
	float cameraX = 180;
	float cameraY = 120;
	float cameraZ = 300;	
	float radius = 300.0f;	
	float theta = 0.45f;
	float phi = 0.45f;
	float a;
	float r;
	//Camera up-down orientation
	glm::vec3 cameraUpDown = glm::vec3(0, 1, 0);
	do {

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT |GL_DEPTH_BUFFER_BIT);

		// Use our programID shader
		glUseProgram(programID);

		// Projection matrix : 45° Field of View, 1:1 ratio, display range : 0.1 unit <-> 500 units
		glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 1.0f / 1.0f, 0.1f, 500.0f);
		
		// Camera matrix
		glm::mat4 View = glm::lookAt(glm::vec3(cameraX,cameraY,cameraZ), // Camera's position in world space
			glm::vec3(50,50,50), // where the camera is looking at (looks at the center of our scene,that is (50,50,50)
			cameraUpDown  // Head is up (0,1,0)
		);
		// Model matrix : an identity matrix (model will be at the origin)
		glm::mat4 Model = glm::mat4(1.0f);
		// Our ModelViewProjection : multiplication of our 3 matrices
		glm::mat4 MVP = Projection * View * Model; 
			
		
		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

		//camera right x
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {				
			theta += 0.01f;	
			cameraX = 50 + radius * sin(theta) * cos(phi);
			cameraZ = 50 + radius * cos(theta) * cos(phi);
		}
		//Camera left X
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
			theta -= 0.01f;
			cameraX = 50 + radius * sin(theta) * cos(phi);
			cameraZ = 50 + radius * cos(theta) * cos(phi);
		}

		//Camera Up 
		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {				
			phi += 0.01f;			
			cameraY = 50 + radius * sin(phi);
			cameraZ = 50 + radius * cos(theta) * cos(phi);
			cameraX = 50 + radius * sin(theta) * cos(phi);
			
		}
		//Camera Down
		if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
			phi -= 0.01f;
			cameraY = 50 + radius * sin(phi);
			cameraZ = 50 + radius * cos(theta) * cos(phi);
			cameraX = 50 + radius * sin(theta) * cos(phi);
			
		}

		//Zoom in
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
			radius -= 1.0f;
			cameraY = 50 + radius * sin(phi);
			cameraZ = 50 + radius * cos(theta) * cos(phi);
			cameraX = 50 + radius * sin(theta) * cos(phi);
		}
		//Zoom out
		if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {			
			radius += 1.0f;
			cameraY = 50 + radius * sin(phi);
			cameraZ = 50 + radius * cos(theta) * cos(phi);
			cameraX = 50 + radius * sin(theta) * cos(phi);
		}
		
		//Small Objects drawing , moving , collisions and refresh shader information
		//Loop goes through all small objects inside our global smallObjects vector
		for (int i = 0; i < numberOfSmallObjects; i += 3) {
			//In every frame for every object translating to the new position
			Model = glm::translate(Model, glm::vec3(direction[i], direction[i + 1], direction[i + 2]));
			//Informing shaders with Model's new position
			MVP = Projection * View * Model;
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
			//Call drawSmallObject function to draw small object
			drawSmallObject(smallObjects[i], smallObjects[i + 1], smallObjects[i + 2], indices_size[index_count]);
			//Handle small object velocity increase/decrease by pressing (SHIFT + '.') = " > " or (SHIFT + ',') = " < " (Smallest velocity = 0.05 and biggest = 10)
			if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS && glfwGetKey(window, GLFW_KEY_COMMA) == GLFW_PRESS) {
				if (velocity[i] > 0.05) {
					velocity[i] -= 0.01;
				}
				if (velocity[i + 1] > 0.05) {
					velocity[i + 1] -= 0.01;
				}
				if (velocity[i + 2] > 0.05) {
					velocity[i + 2] -= 0.01;
				}				
			}
			if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS && glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_PRESS) {
				if (velocity[i] < 10) {
					velocity[i] += 0.01;
				}
				if (velocity[i + 1] < 10) {
					velocity[i + 1] += 0.01;
				}
				if (velocity[i + 2] < 10) {
					velocity[i + 2] += 0.01;
				}
			}
			//Putting the new velocity in every frame in x,y and z coordinates.
			direction[i] += velocity[i];
			direction[i + 1] += velocity[i + 1];
			direction[i + 2] += velocity[i + 2];
			//The three following function handle small object collision with the wall, the big sphere and between all of them
			smallObjectCollisionWithWall(i, index_count);
			smallObjectCollisionWithSphere(i, index_count);
			smallObjectCollisionWithSmallObject(i, index_count);
			index_count++;	
			//At the end of the loop with translate model matrix back to 1.0 for the next object to translate;
			//If we used different model matrices we would not need to do this
			Model = glm::mat4(1.0f);
			MVP = Projection * View * Model;
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
		}
		index_count = 0;
		
		// Translate model matrix to take big sphere's position
		Model = glm::translate(Model, glm::vec3(SpherePosition[0], SpherePosition[1], SpherePosition[2]));
		MVP = Projection * View * Model;
		//Here we use the big sphere shader to handle sphere texture and mvp matrix
		glUseProgram(textureShaderID);
		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);
		//Send texture id to shader
		glUniform1i(TextureID, 0);
		//Send to shader flag to enable/disable texture
		glUniform1i(OpacityID,textureOnOff);
		//Send MVP matrix to 
		glUniformMatrix4fv(Matrix2ID, 1, GL_FALSE, &MVP[0][0]);
		//drawBigSphereObject draws our big sphere 
		drawBigSphereObject(vertexbufferBigSphere, elementbufferBigSphere,texturebufferBigSphere, bigSphereIndicesCount);
		
		//After we finish with the sphere texture we use again the other shader programID
		glUseProgram(programID);
		//Handle 3-axis sphere movement with keyboard
		// Move Up y
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
			SpherePosition[1] += BigSphereVelocity;
		}
		// Move Down y
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
			SpherePosition[1] -= BigSphereVelocity;
		}
		// Move right x
		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			SpherePosition[0] += BigSphereVelocity;
		}
		// Move left x
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
			SpherePosition[0] -= BigSphereVelocity;
		}// Move  forward z
		if (glfwGetKey(window, GLFW_KEY_KP_ADD) == GLFW_PRESS) {
			SpherePosition[2] += BigSphereVelocity;
		}
		// Move backward z
		if (glfwGetKey(window, GLFW_KEY_KP_SUBTRACT) == GLFW_PRESS) {
			SpherePosition[2] -= BigSphereVelocity;
		}
		//handles collision of big sphere with cube walls
		BigSphereCollisionWithWall();
		//Restore Model matrix so other objects can tranform it
		Model = glm::mat4(1.0f);
		MVP = Projection * View * Model;
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
		
		//Draw scene cube object
		//Gl blend called to handle cube transparency
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		//Enable vertex buffer and go to 0 location of our vertex shader
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, sc_vertexbuffer);
		glVertexAttribPointer(0,3, GL_FLOAT,GL_FALSE,0,(void*)0 );
		//Enable color buffer and go to 1 location of our vertex shader
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, sc_colorbuffer);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, (void*)0);
		//Draw cube triangles
		glDrawArrays(GL_TRIANGLES, 0, 12 * 3); // 12*3 indices starting at 0 -> 12 triangles
		glCullFace(GL_FRONT);
		glDisable(GL_CULL_FACE);
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		
		//Draw scene cube wireframe
		//Enable wireframe vertex buffer and go to 0 location of our vertex shader
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, wireframe_vertexbuffer);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
		//Enable wireframe color buffer and go to 1 location of our vertex shader
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, wireframe_colorbuffer);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);	
		//Draw line to create the edges of the cube
		glDrawArrays(GL_LINES, 0, 12*2); 	
		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO and shader
	glDeleteBuffers(1, &sc_vertexbuffer);
	glDeleteBuffers(1, &sc_colorbuffer);
	glDeleteBuffers(1, &wireframe_vertexbuffer);
	glDeleteBuffers(1, &wireframe_colorbuffer);
	glDeleteBuffers(1, &vertexbufferBigSphere);
	glDeleteBuffers(1, &elementbufferBigSphere);
	glDeleteBuffers(1, &texturebufferBigSphere);
	glDeleteProgram(programID);
	glDeleteProgram(textureShaderID);
	glDeleteVertexArrays(1, &VertexArrayID);
	for (int i = 0; i < numberOfSmallObjects; i += 3) {
		glDeleteBuffers(1 ,&smallObjects[i]);
		glDeleteBuffers(1 ,&smallObjects[i + 1]);
		glDeleteBuffers(1 ,&smallObjects[i + 2]);
	}
	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}

