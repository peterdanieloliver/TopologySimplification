#pragma once
#include "polyhedron.h"
#include "icVector.H"
#include <cmath>
#include <vector>
#include <math.h>

class Cluster
{
private:

	std::vector<icVector3*> singularities;		// set of singularities in the cluster
	std::vector<Quad*> cell;					// set of quads containing the singularity cluster

public:

	// constructor
	Cluster(std::vector<icVector3*> points, std::vector<Quad*> quads)
	{
		singularities = points;
		cell = quads;
	}
};

class SingClusterHandeler
{
private:
	
	Polyhedron* poly;
	std::vector<icVector3*> singularities;
	std::vector<Cluster*> clusters;
	double cluster_threshold;

public:

	// constructor
	SingClusterHandeler(Polyhedron* poly_in, std::vector<icVector3*> points, double threshold)
	{
		poly = poly_in;
		singularities = points;
		cluster_threshold = threshold;
	}

};