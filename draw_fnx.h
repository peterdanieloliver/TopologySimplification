#pragma once

#include <vector>
#include "glError.h"
#include "polyline.h"


/*
draw a sphere
x, y, z are the coordiate of the dot
radius of the sphere
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawDot(double x, double y, double z, double radius = 0.1, double R = 1.0, double G = 0.0, double B = 0.0) {

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1., 1.);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	GLfloat mat_diffuse[4];

	{
		mat_diffuse[0] = R;
		mat_diffuse[1] = G;
		mat_diffuse[2] = B;
		mat_diffuse[3] = 1.0;
	}

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

	GLUquadric* quad = gluNewQuadric();

	glPushMatrix();
	glTranslatef(x, y, z);
	gluSphere(quad, radius, 50, 50);
	glPopMatrix();

	gluDeleteQuadric(quad);
}

/*
draw a line segment
width: the width of the line, should bigger than 0
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawLineSegment(LineSegment ls, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {

	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);
	glVertex3f(ls.start.x, ls.start.y, ls.start.z);
	glVertex3f(ls.end.x, ls.end.y, ls.end.z);
	glEnd();

	glDisable(GL_BLEND);
}

/*
draw a polyline
width: the width of the line, should bigger than 0
R: the red channel of the color, ranges [0, 1]
G: the green channel of the color, ranges [0, 1]
B: the blue channel of the color, ranges [0, 1]
*/
void drawPolyline(PolyLine pl, double width = 1.0, double R = 1.0, double G = 0.0, double B = 0.0) {

	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(width);

	glBegin(GL_LINES);
	glColor3f(R, G, B);

	for (int i = 0; i < pl.size(); i++) {
		glVertex3f(pl[i].start.x, pl[i].start.y, pl[i].start.z);
		glVertex3f(pl[i].end.x, pl[i].end.y, pl[i].end.z);
	}

	glEnd();

	glDisable(GL_BLEND);
}