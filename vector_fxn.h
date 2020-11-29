#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <vector>

#include "glError.h"
#include "gl/glew.h"
#include "gl/freeglut.h"
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "polyhedron.h"
#include "polyline.h"
#include "trackball.h"
#include "tmatrix.h"

const double TRAVEL_THRESHOLD = 0.2;
const double STEP_SIZE = 0.07;
const double SEP_OFFSET = 0.07;
const double WINDING_NUMBER_TOLERANCE = 6;
const double S_T_TOLERANCE = 0.06;

// maximum helper function
template <class T> const T& max3(const T& a, const T& b, const T& c)
{
	return std::max(std::max(a, b), c);
}

// maximum helper function
template <class T> const T& min3(const T& a, const T& b, const T& c)
{
	return std::min(std::min(a, b), c);
}

// recursive streamline stepping function (for forward parameter, true is forward and false is backward)
void streamline_step(double x0, double y0, Quad* quad, PolyLine* streamline, Polyhedron* poly, bool forward)
{
	// check for nullptr
	if (quad == nullptr) { return; }

	double x1, x2, y1, y2;
	double vx11, vx12, vx21, vx22, vy11, vy12, vy21, vy22;
	double vx, vy, z0;
	double x, y, z;
	Quad* new_quad;

	// extracting x and y variables
	x1 = min3(quad->verts[0]->x, quad->verts[1]->x, quad->verts[2]->x);
	x2 = max3(quad->verts[0]->x, quad->verts[1]->x, quad->verts[2]->x);
	y1 = min3(quad->verts[0]->y, quad->verts[1]->y, quad->verts[2]->y);
	y2 = max3(quad->verts[0]->y, quad->verts[1]->y, quad->verts[2]->y);
	z = quad->verts[0]->z;

	// extracting vector values
	for (Vertex* vert : quad->verts)
	{
		if ((vert->x == x1) && (vert->y == y1))
		{
			vx11 = vert->vx;
			vy11 = vert->vy;
		}
		else if (vert->x == x1)
		{
			vx12 = vert->vx;
			vy12 = vert->vy;
		}
		else if (vert->y == y1)
		{
			vx21 = vert->vx;
			vy21 = vert->vy;
		}
		else
		{
			vx22 = vert->vx;
			vy22 = vert->vy;
		}
	}

	// interpolate
	vx = ((x2 - x0) * (y2 - y0) * vx11) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y2 - y0) * vx21) / ((x2 - x1) * (y2 - y1))
		+ ((x2 - x0) * (y0 - y1) * vx12) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y0 - y1) * vx22) / ((x2 - x1) * (y2 - y1));
	vy = ((x2 - x0) * (y2 - y0) * vy11) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y2 - y0) * vy21) / ((x2 - x1) * (y2 - y1))
		+ ((x2 - x0) * (y0 - y1) * vy12) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y0 - y1) * vy22) / ((x2 - x1) * (y2 - y1));

	// check direction
	if (!forward)
	{
		vx = -1 * vx;
		vy = -1 * vy;
	}

	// normalize vectors
	double magnitude = sqrt((vx * vx) + (vy * vy));
	vx = vx / magnitude;
	vy = vy / magnitude;

	// calculate new position
	x = x0 + (STEP_SIZE * vx);
	y = y0 + (STEP_SIZE * vy);
	

	// check to see if the new point is outside the quad. If it is,
	// find the crossing point and make that the new point. Then 
	// determine which edge this crossing point lies on.
	if (x < x1 || x > x2 || y < y1 || y > y2)
	{
		// extracting edges from quad
		Edge* ex1;
		Edge* ex2;
		Edge* ey1;
		Edge* ey2;
		for (Edge* edge : quad->edges)
		{
			if ((edge->verts[0]->x == x1) && (edge->verts[1]->x == x1))
			{
				ex1 = edge;
			}
			else if ((edge->verts[0]->x == x2) && (edge->verts[1]->x == x2))
			{
				ex2 = edge;
			}
			else if (edge->verts[0]->y == y1)
			{
				ey1 = edge;
			}
			else
			{
				ey2 = edge;
			}
		}

		// find crossing point and corresponding edge
		double x_cross, y_cross;
		Edge* cross_edge = nullptr;
		if (vx > 0 && vy > 0)
		{
			x_cross = x0 + (((y2 - y0) / vy) * vx);
			y_cross = y0 + (((x2 - x0) / vx) * vy);
			if (x_cross < x2)
			{
				x = x_cross;
				y = y2;
				cross_edge = ey2;
			}
			else if (y_cross < y2)
			{
				x = x2;
				y = y_cross;
				cross_edge = ex2;
			}
			else
			{
				x = x2;
				y = y2;
				cross_edge = ex2;
			}
		}
		else if (vx > 0 && vy < 0)
		{
			x_cross = x0 + (((y1 - y0) / vy) * vx);
			y_cross = y0 + (((x2 - x0) / vx) * vy);
			if (x_cross < x2)
			{
				x = x_cross;
				y = y1;
				cross_edge = ey1;
			}
			else if (y_cross > y1)
			{
				x = x2;
				y = y_cross;
				cross_edge = ex2;
			}
			else
			{
				x = x2;
				y = y1;
				cross_edge = ex2;
			}
		}
		else if (vx < 0 && vy > 0)
		{
			x_cross = x0 + (((y2 - y0) / vy) * vx);
			y_cross = y0 + (((x1 - x0) / vx) * vy);
			if (x_cross > x1)
			{
				x = x_cross;
				y = y2;
				cross_edge = ey2;
			}
			else if (y_cross < y2)
			{
				x = x1;
				y = y_cross;
				cross_edge = ex1;
			}
			else
			{
				x = x1;
				y = y2;
				cross_edge = ex1;
			}
		}
		else if (vx < 0 && vy < 0)
		{
			x_cross = x0 + (((y1 - y0) / vy) * vx);
			y_cross = y0 + (((x1 - x0) / vx) * vy);
			if (x_cross > x1)
			{
				x = x_cross;
				y = y1;
				cross_edge = ey1;
			}
			else if (y_cross > y1)
			{
				x = x1;
				y = y_cross;
				cross_edge = ex1;
			}
			else
			{
				x = x1;
				y = y1;
				cross_edge = ex1;
			}
		}
		else if (vx == 0 && vy > 0)
		{
			x_cross = x0;
			y_cross = y2;
			cross_edge = ey2;
		}
		else if (vx == 0 && vy < 0)
		{
			x_cross = x0;
			y_cross = y1;
			cross_edge = ey1;
		}
		else if (vy == 0 && vx > 0)
		{
			x_cross = x2;
			y_cross = y0;
			cross_edge = ex2;
		}
		else if (vy == 0 && vx < 0)
		{
			x_cross = x1;
			y_cross = y0;
			cross_edge = ex1;
		}
		else
		{
			return;
		}
		
		// find new quad, if at edge of domain, set to nullptr
		if (cross_edge != nullptr)
		{
			if (cross_edge->nquads == 2)
			{
				if (cross_edge->quads[0] == quad)
				{
					new_quad = cross_edge->quads[1];
				}
				else
				{
					new_quad = cross_edge->quads[0];
				}
			}
			else
			{
				new_quad = nullptr;
			}
		}
	}
	else
	{
		new_quad = quad;
	}

	streamline->push_back(LineSegment(icVector3(x0,y0,z),icVector3(x, y, z)));

	if (magnitude*STEP_SIZE < TRAVEL_THRESHOLD)
	{
		return;
	}
	else
	{
		streamline_step(x, y, new_quad, streamline, poly, forward);
	}
}

// streamline building function
PolyLine* build_streamline(double x, double y, Polyhedron* poly)
{
	PolyLine* streamline = new PolyLine();
	Quad* q_start = nullptr;

	// loop through quads to find the starting point
	double x1, x2, y1, y2;
	Quad* q_temp = nullptr;
	for (int i = 0; i < poly->nquads; i++)
	{
		q_temp = poly->qlist[i];
		x1 = min3(q_temp->verts[0]->x, q_temp->verts[1]->x, q_temp->verts[2]->x);
		x2 = max3(q_temp->verts[0]->x, q_temp->verts[1]->x, q_temp->verts[2]->x);
		y1 = min3(q_temp->verts[0]->y, q_temp->verts[1]->y, q_temp->verts[2]->y);
		y2 = max3(q_temp->verts[0]->y, q_temp->verts[1]->y, q_temp->verts[2]->y);

		if (x >= x1 && x <= x2 && y >= y1 && y <= y2)
		{
			q_start = q_temp;
			break;
		}
	}

	if (q_start == nullptr)
	{
		return nullptr;
	}

	// build the forward part of the streamline
	streamline_step(x, y, q_start, streamline, poly, true);
	// build the backward part of the streamline
	streamline_step(x, y, q_start, streamline, poly, false);

	return streamline;
}

Quad* streamline_step2(icVector2* point, Quad* quad, PolyLine* streamline, Polyhedron* poly, bool forward)
{
	// check for nullptr
	if (quad == nullptr) { return nullptr; }

	double x0, y0;
	double x1, x2, y1, y2;
	double vx11, vx12, vx21, vx22, vy11, vy12, vy21, vy22;
	double vx, vy;
	double x, y, z;
	Quad* new_quad;

	x0 = point->x;
	y0 = point->y;

	// extracting x and y variables
	x1 = min3(quad->verts[0]->x, quad->verts[1]->x, quad->verts[2]->x);
	x2 = max3(quad->verts[0]->x, quad->verts[1]->x, quad->verts[2]->x);
	y1 = min3(quad->verts[0]->y, quad->verts[1]->y, quad->verts[2]->y);
	y2 = max3(quad->verts[0]->y, quad->verts[1]->y, quad->verts[2]->y);
	z = quad->verts[0]->z;

	// extracting vector values
	for (Vertex* vert : quad->verts)
	{
		if ((vert->x == x1) && (vert->y == y1))
		{
			vx11 = vert->vx;
			vy11 = vert->vy;
		}
		else if (vert->x == x1)
		{
			vx12 = vert->vx;
			vy12 = vert->vy;
		}
		else if (vert->y == y1)
		{
			vx21 = vert->vx;
			vy21 = vert->vy;
		}
		else
		{
			vx22 = vert->vx;
			vy22 = vert->vy;
		}
	}

	// interpolate
	vx = ((x2 - x0) * (y2 - y0) * vx11) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y2 - y0) * vx21) / ((x2 - x1) * (y2 - y1))
		+ ((x2 - x0) * (y0 - y1) * vx12) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y0 - y1) * vx22) / ((x2 - x1) * (y2 - y1));
	vy = ((x2 - x0) * (y2 - y0) * vy11) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y2 - y0) * vy21) / ((x2 - x1) * (y2 - y1))
		+ ((x2 - x0) * (y0 - y1) * vy12) / ((x2 - x1) * (y2 - y1))
		+ ((x0 - x1) * (y0 - y1) * vy22) / ((x2 - x1) * (y2 - y1));

	// check direction
	if (!forward)
	{
		vx = -1 * vx;
		vy = -1 * vy;
	}

	// normalize vectors
	double magnitude = sqrt((vx * vx) + (vy * vy));
	vx = vx / magnitude;
	vy = vy / magnitude;

	// calculate new position
	x = x0 + (STEP_SIZE * vx);
	y = y0 + (STEP_SIZE * vy);

	if (x < x1 || x > x2 || y < y1 || y > y2)
	{
		// extracting edges from quad
		Edge* ex1;
		Edge* ex2;
		Edge* ey1;
		Edge* ey2;
		for (Edge* edge : quad->edges)
		{
			if ((edge->verts[0]->x == x1) && (edge->verts[1]->x == x1))
			{
				ex1 = edge;
			}
			else if ((edge->verts[0]->x == x2) && (edge->verts[1]->x == x2))
			{
				ex2 = edge;
			}
			else if (edge->verts[0]->y == y1)
			{
				ey1 = edge;
			}
			else
			{
				ey2 = edge;
			}
		}

		// find crossing point and corresponding edge
		double x_cross, y_cross;
		Edge* cross_edge = nullptr;
		if (vx > 0 && vy > 0)
		{
			x_cross = x0 + (((y2 - y0) / vy) * vx);
			y_cross = y0 + (((x2 - x0) / vx) * vy);
			if (x_cross < x2)
			{
				x = x_cross;
				y = y2;
				cross_edge = ey2;
			}
			else if (y_cross < y2)
			{
				x = x2;
				y = y_cross;
				cross_edge = ex2;
			}
			else
			{
				x = x2;
				y = y2;
				cross_edge = ex2;
			}
		}
		else if (vx > 0 && vy < 0)
		{
			x_cross = x0 + (((y1 - y0) / vy) * vx);
			y_cross = y0 + (((x2 - x0) / vx) * vy);
			if (x_cross < x2)
			{
				x = x_cross;
				y = y1;
				cross_edge = ey1;
			}
			else if (y_cross > y1)
			{
				x = x2;
				y = y_cross;
				cross_edge = ex2;
			}
			else
			{
				x = x2;
				y = y1;
				cross_edge = ex2;
			}
		}
		else if (vx < 0 && vy > 0)
		{
			x_cross = x0 + (((y2 - y0) / vy) * vx);
			y_cross = y0 + (((x1 - x0) / vx) * vy);
			if (x_cross > x1)
			{
				x = x_cross;
				y = y2;
				cross_edge = ey2;
			}
			else if (y_cross < y2)
			{
				x = x1;
				y = y_cross;
				cross_edge = ex1;
			}
			else
			{
				x = x1;
				y = y2;
				cross_edge = ex1;
			}
		}
		else if (vx < 0 && vy < 0)
		{
			x_cross = x0 + (((y1 - y0) / vy) * vx);
			y_cross = y0 + (((x1 - x0) / vx) * vy);
			if (x_cross > x1)
			{
				x = x_cross;
				y = y1;
				cross_edge = ey1;
			}
			else if (y_cross > y1)
			{
				x = x1;
				y = y_cross;
				cross_edge = ex1;
			}
			else
			{
				x = x1;
				y = y1;
				cross_edge = ex1;
			}
		}
		else if (vx == 0 && vy > 0)
		{
			x_cross = x0;
			y_cross = y2;
			cross_edge = ey2;
		}
		else if (vx == 0 && vy < 0)
		{
			x_cross = x0;
			y_cross = y1;
			cross_edge = ey1;
		}
		else if (vy == 0 && vx > 0)
		{
			x_cross = x2;
			y_cross = y0;
			cross_edge = ex2;
		}
		else if (vy == 0 && vx < 0)
		{
			x_cross = x1;
			y_cross = y0;
			cross_edge = ex1;
		}
		else
		{
			return nullptr;
		}

		// find new quad, if at edge of domain, set to nullptr
		if (cross_edge != nullptr)
		{
			if (cross_edge->nquads == 2)
			{
				if (cross_edge->quads[0] == quad)
				{
					new_quad = cross_edge->quads[1];
				}
				else
				{
					new_quad = cross_edge->quads[0];
				}
			}
			else
			{
				new_quad = nullptr;
			}
		}
	}
	else
	{
		new_quad = quad;
	}

	streamline->push_back(LineSegment(icVector3(x0, y0, z), icVector3(x, y, z)));

	if (magnitude < TRAVEL_THRESHOLD)
	{
		return nullptr;
	}
	else
	{
		point->x = x;
		point->y = y;
		return new_quad;
	}
}

PolyLine* build_streamline2(double x, double y, Polyhedron* poly)
{
	PolyLine* streamline = new PolyLine();
	icVector2* p_start = new icVector2();
	p_start->x = x;
	p_start->y = y;
	Quad* q_start = nullptr;

	// loop through quads to find the starting point
	double x1, x2, y1, y2;
	Quad* q_temp = nullptr;
	for (int i = 0; i < poly->nquads; i++)
	{
		q_temp = poly->qlist[i];
		x1 = min3(q_temp->verts[0]->x, q_temp->verts[1]->x, q_temp->verts[2]->x);
		x2 = max3(q_temp->verts[0]->x, q_temp->verts[1]->x, q_temp->verts[2]->x);
		y1 = min3(q_temp->verts[0]->y, q_temp->verts[1]->y, q_temp->verts[2]->y);
		y2 = max3(q_temp->verts[0]->y, q_temp->verts[1]->y, q_temp->verts[2]->y);

		if (x >= x1 && x <= x2 && y >= y1 && y <= y2)
		{
			q_start = q_temp;
			break;
		}
	}

	if (q_start == nullptr)
	{
		return nullptr;
	}

	// forward streamline
	while (q_start != nullptr)
	{
		q_start = streamline_step2(p_start, q_start, streamline, poly, true);	
	}
	
	// backward streamline
	q_start = q_temp;
	p_start->x = x;
	p_start->y = y;
	while (q_start != nullptr)
	{
		q_start = streamline_step2(p_start, q_start, streamline, poly, false);
	}

	return streamline;
}

// calculates winding number over an edge
double winding_num(double vx1, double vy1, double vx2, double vy2)
{
	double alpha1, alpha2;
	double winding_num;

	// find direction of each vector
	alpha1 = atan2(vy1, vx1);
	alpha2 = atan2(vy2, vx2);

	// find the winding number
	winding_num = alpha2 - alpha1;
	if (winding_num > PI)
	{
		winding_num = winding_num - (2 * PI);
	}
	else if (winding_num < -PI)
	{
		winding_num = winding_num + (2 * PI);
	}
	
	return winding_num;
}

// 2x2 and 2x1 matrix multiplication
double * matrix_vect_mult(double m11, double m12, double m21, double m22, double v1, double v2)
{
	double* result = new double[2];

	result[0] = (m11 * v1) + (m12 * v2);
	result[1] = (m21 * v1) + (m22 * v2);

	return result;
}

// calculates major and minor eigenvectors of 2x2 matrix (major evect is first in array)
double * eigen_vects(double a, double b, double c, double d)
{
	double* eigen_vects = new double[4];
	double gamma_r, gamma_s, theta, phi;
	double* eigen_major;
	double* eigen_minor;

	gamma_r = (c - b) / 2;
	gamma_s = sqrt(((a - d) * (a - d)) + ((b + c) * (b + c))) / 2;

	theta = atan2((b + c), (a - d));
	phi = atan(gamma_r / gamma_s);
	
	double m11, m12, m21, m22, v1, v2;
	
	m11 = cos(theta / 2);
	m12 = -sin(theta / 2);
	m21 = sin(theta / 2);
	m22 = cos(theta / 2);

	v1 = sqrt(sin(phi + (PI / 4))) + sqrt(cos(phi + (PI / 4)));
	v2 = sqrt(sin(phi + (PI / 4))) - sqrt(cos(phi + (PI / 4)));

	eigen_major = matrix_vect_mult(m11, m12, m21, m22, v1, v2);
	eigen_minor = matrix_vect_mult(m11, m12, m21, m22, v2, v1);

	eigen_vects[0] = eigen_major[0];
	eigen_vects[1] = eigen_major[1];
	eigen_vects[2] = eigen_minor[0];
	eigen_vects[3] = eigen_minor[1];

	return eigen_vects;
}

// find singularity in a quad
void find_singularity(Quad* quad, std::vector<icVector3*>* source_sink_points, std::vector<icVector3*>* saddle_points, std::vector<icVector3*>* separatrix_points)
{
	double x1, x2, y1, y2;
	double vx11, vx12, vx21, vx22, vy11, vy12, vy21, vy22, f11, f12, f21, f22;
	double winding_number;
	bool source_sink = false;
	bool saddle = false;
	double a00, a10, a01, a11, b00, b10, b01, b11, c00, c10, c01;
	double s, t;
	double a, b, c;
	double discriminant;
	double x0, y0, z0;
	icVector3* singularity;

	// extracting vertex coordinates
	x1 = min3(quad->verts[0]->x, quad->verts[1]->x, quad->verts[2]->x);
	x2 = max3(quad->verts[0]->x, quad->verts[1]->x, quad->verts[2]->x);
	y1 = min3(quad->verts[0]->y, quad->verts[1]->y, quad->verts[2]->y);
	y2 = max3(quad->verts[0]->y, quad->verts[1]->y, quad->verts[2]->y);

	// extracting vertex vector values
	for (Vertex* vert : quad->verts)
	{
		if ((vert->x == x1) && (vert->y == y1))
		{
			vx11 = vert->vx;
			vy11 = vert->vy;
			f11 = vert->scalar;
		}
		else if (vert->x == x1)
		{
			vx12 = vert->vx;
			vy12 = vert->vy;
			f12 = vert->scalar;
		}
		else if (vert->y == y1)
		{
			vx21 = vert->vx;
			vy21 = vert->vy;
			f21 = vert->scalar;
		}
		else
		{
			vx22 = vert->vx;
			vy22 = vert->vy;
			f22 = vert->scalar;
		}
	}

	// compute winding number for quad 
	winding_number = winding_num(vx11, vy11, vx21, vy21)
				   + winding_num(vx21, vy21, vx22, vy22)
				   + winding_num(vx22, vy22, vx12, vy12)
				   + winding_num(vx12, vy12, vx11, vy11);
	
	if ((winding_number > ((2*PI) - WINDING_NUMBER_TOLERANCE)) && (winding_number < ((2*PI) + WINDING_NUMBER_TOLERANCE)))
	{
		// it is a source, sink, center or focus
		source_sink = true;
	}
	else if ((winding_number > ((-2 * PI) - WINDING_NUMBER_TOLERANCE)) && (winding_number < ((-2 * PI) + WINDING_NUMBER_TOLERANCE)))
	{
		// it is a saddle
		saddle = true;
	}
	else
	{
		return;
	}

	if (source_sink || saddle)
	{
		// computing coefficients
		a00 = vx11;
		a10 = vx21 - vx11;
		a01 = vx12 - vx11;
		a11 = vx11 - vx21 - vx12 + vx22;

		b00 = vy11;
		b10 = vy21 - vy11;
		b01 = vy12 - vy11;
		b11 = vy11 - vy21 - vy12 + vy22;

		c00 = (a11 * b00) - (a00 * b11);
		c10 = (a11 * b10) - (a10 * b11);
		c01 = (a11 * b01) - (a01 * b11);

		a = (-a11 * c10);
		b = (-a11 * c00) - (a01 * c10) + (a10 * c01);
		c = (a00 * c01) - (a01 * c00);

		// use quadratic formula to find s
		discriminant = (b * b) - (4 * a * c);
		if (discriminant == 0)
		{
			s = (-b / (2 * a));
			t = (-c00 / c01) - ((c10 / c01) * s);
			if (s <= (-S_T_TOLERANCE) || s >= (1 + S_T_TOLERANCE) || t <= (-S_T_TOLERANCE) || t >= (1 + S_T_TOLERANCE)) { return; }
		}
		else if (discriminant > 0)
		{
			s = (-b + sqrt(discriminant)) / (2 * a);
			if (s <= (-S_T_TOLERANCE) || s >= (1 + S_T_TOLERANCE))
			{
				s = (-b - sqrt(discriminant)) / (2 * a);
				if (s <= (-S_T_TOLERANCE) || s >= (1 + S_T_TOLERANCE)) { return; }
			}
			
			t = (-c00 / c01) - ((c10 / c01) * s);
			if (t <= (- S_T_TOLERANCE) || t >= (1 + S_T_TOLERANCE)) { return; }
		}
		else
		{
			return;
		}

		// use s and t to find singularity coordinates
		x0 = ((x2 - x1) * s) + x1;
		y0 = ((y2 - y1) * t) + y1;
		z0 = quad->verts[0]->z;

		singularity = new icVector3(x0, y0, z0);

		if (saddle)
		{
			double dfdx, dfdy, dgdx, dgdy;
			// calculate jacobian values
			dfdx = (-(y2 - y0) * vx11) / ((x2 - x1) * (y2 - y1))
				+ ((y2 - y0) * vx21) / ((x2 - x1) * (y2 - y1))
				+ (-(y0 - y1) * vx12) / ((x2 - x1) * (y2 - y1))
				+ ((y0 - y1) * vx22) / ((x2 - x1) * (y2 - y1));
			dfdy = (-(x2 - x0) * vx11) / ((x2 - x1) * (y2 - y1))
				+ (-(x0 - x1) * vx21) / ((x2 - x1) * (y2 - y1))
				+ ((x2 - x0) * vx12) / ((x2 - x1) * (y2 - y1))
				+ ((x0 - x1) * vx22) / ((x2 - x1) * (y2 - y1));
			dgdx = (-(y2 - y0) * vy11) / ((x2 - x1) * (y2 - y1))
				+ ((y2 - y0) * vy21) / ((x2 - x1) * (y2 - y1))
				+ (-(y0 - y1) * vy12) / ((x2 - x1) * (y2 - y1))
				+ ((y0 - y1) * vy22) / ((x2 - x1) * (y2 - y1));
			dgdy = (-(x2 - x0) * vy11) / ((x2 - x1) * (y2 - y1))
				+ (-(x0 - x1) * vy21) / ((x2 - x1) * (y2 - y1))
				+ ((x2 - x0) * vy12) / ((x2 - x1) * (y2 - y1))
				+ ((x0 - x1) * vy22) / ((x2 - x1) * (y2 - y1));

			double* evects = eigen_vects(dfdx, dfdy, dgdx, dgdy);

			icVector3* sep_out_1 = new icVector3((x0 + (evects[0] * SEP_OFFSET)), 
												 (y0 + (evects[1] * SEP_OFFSET)), 0.0);
			icVector3* sep_out_2 = new icVector3((x0 - (evects[0] * SEP_OFFSET)),
												 (y0 - (evects[1] * SEP_OFFSET)), 0.0);
			icVector3* sep_in_1 = new icVector3((x0 + (evects[2] * SEP_OFFSET)),
												(y0 + (evects[3] * SEP_OFFSET)), 0.0);
			icVector3* sep_in_2 = new icVector3((x0 - (evects[2] * SEP_OFFSET)),
												(y0 - (evects[3] * SEP_OFFSET)), 0.0);

			separatrix_points->push_back(sep_out_1);
			separatrix_points->push_back(sep_out_2);
			separatrix_points->push_back(sep_in_1);
			separatrix_points->push_back(sep_in_2);
			saddle_points->push_back(singularity);
		}
		else if (source_sink)
		{
			source_sink_points->push_back(singularity);
		}
		else
		{
			return;
		}
	}
}

void generate_topology(Polyhedron* poly, std::vector<icVector3*>* source_sink_points, std::vector<icVector3*>* saddle_points, std::vector<PolyLine*>* separatrices)
{
	std::vector<icVector3*>* separatrix_points = new std::vector<icVector3*>;

	Quad* quad;
	for (int i = 0; i < poly->nquads; i++)
	{
		quad = poly->qlist[i];

		if (i == 228)
		{
			int test = 1;
		}

		find_singularity(quad, source_sink_points, saddle_points, separatrix_points);
	}
	for (icVector3* sep_point : *separatrix_points)
	{
		PolyLine* separatrix = build_streamline2(sep_point->x, sep_point->y, poly);
		separatrices->push_back(separatrix);
	}
	return;
}
