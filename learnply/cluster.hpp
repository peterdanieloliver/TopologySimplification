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

	std::vector<icVector3*> getSingularities() {
		return singularities;
	}

	std::vector<Quad*> getCell() {
		return cell;
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
		std::vector<Quad*> cell = std::vector<Quad*>();
		
		for (int i = 0; i < poly->nquads; ++i) {
			cell.push_back(poly->qlist[i]);
		}

		Cluster* mesh = new Cluster(singularities, cell);
		clusters = separateClusters(mesh);
	}

	std::vector<Cluster*> separateClusters(Cluster* c) {
		double s = DBL_MAX;
		Vertex* min_v = nullptr;

		//step 1 && step 2
		for (int i = 0; i < poly->nverts; ++i) {
			Vertex* v = poly->vlist[i];
			if (v->nquads == 4) {
				double sum = 0.0;
				for (icVector3* s : c->getSingularities()) {
					sum += sqrt(pow(v->x - s->x, 2) + pow(v->y - s->y, 2));
				}
				if (sum < s) {
					s = sum;
					min_v = v;
				}
			}
		}

		if (s > cluster_threshold) {
			//step 3: compute coordinate axis variance
			double variance_x = 0.0, variance_y = 0.0;
		
			for (icVector3* s : c->getSingularities()) {
				variance_x += pow(min_v->x - s->x, 2);
				variance_y += pow(min_v->y - s->y, 2);
			}
			
			std::vector<icVector3*> leftSings = std::vector<icVector3*>();
			std::vector<Quad*> leftCell = std::vector<Quad*>();
			std::vector<icVector3*> rightSings = std::vector<icVector3*>();
			std::vector<Quad*> rightCell = std::vector<Quad*>();

			if (variance_x > variance_y) {
				std::vector<icVector3*> clusterSings = c->getSingularities();
				for (icVector3* s : clusterSings) {
					if (s->x > min_v->x) {
						rightSings.push_back(s);
					}
					else {
						leftSings.push_back(s);
					}
				}

				std::vector<Quad*> clusterCell = c->getCell();
				for (Quad* q : clusterCell) {
					double avg = (q->verts[0]->x + q->verts[1]->x + q->verts[2]->x + q->verts[3]->x) / 4;
					if (avg > min_v->x) {
						rightCell.push_back(q);
					}
					else {
						leftCell.push_back(q);
					}
				}
			}
			else {
				std::vector<icVector3*> clusterSings = c->getSingularities();
				for (icVector3* s : clusterSings) {
					if (s->y > min_v->y) {
						rightSings.push_back(s);
					}
					else {
						leftSings.push_back(s);
					}
				}

				std::vector<Quad*> clusterCell = c->getCell();
				for (Quad* q : clusterCell) {
					double avg = (q->verts[0]->y + q->verts[1]->y + q->verts[2]->y + q->verts[3]->y) / 4;
					if (avg > min_v->y) {
						rightCell.push_back(q);
					}
					else {
						leftCell.push_back(q);
					}
				}
			}

			std::vector<Cluster*> left = separateClusters(new Cluster(leftSings, leftCell));
			std::vector<Cluster*> right = separateClusters(new Cluster(rightSings, rightCell));

			std::vector<Cluster*> join = std::vector<Cluster*>();
			join.reserve(left.size() + right.size());
			join.insert(join.end(), left.begin(), left.end());
			join.insert(join.end(), right.begin(), right.end());
			return join;
		}
		else {
			std::vector<Cluster*> ret = std::vector<Cluster*>();
			ret.push_back(c);
			return ret;
		}
	}

	std::vector<Cluster*> getClusters() {
		return clusters;
	}

};