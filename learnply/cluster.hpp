#pragma once
#include "polyhedron.h"
#include "icVector.H"
#include <cmath>
#include <vector>
#include <unordered_set>
#include <math.h>

class Cluster
{
private:

	std::vector<icVector3*> singularities;		// set of singularities in the cluster
	std::vector<Quad*> cell_quads;				// set of quads containing the singularity cluster
	std::unordered_set<Vertex*> perim_verts;	// set of vertices on the perimeter of the cluster cell
	std::unordered_set<Edge*> perim_edges;		// set of edges on the perimeter of the cluster cell
	double min_x, min_y, max_x, max_y;			// boundary coordinate values of the cluster cell
	double sing_max_x, sing_max_y, sing_min_x, sing_min_y; // singularity max values
	icVector3* mean_point;						// new singularity at singularity cluster mean

	// private helper functions

	void find_mean()
	{
		double xx = 0.0; double yy = 0.0;
		// find the mean point of the cluster
		for (icVector3* sing : singularities)
		{
			xx += sing->x;
			yy += sing->y;
			if (sing->x < sing_min_x) { sing_min_x = sing->x; }
			if (sing->x > sing_max_x) { sing_max_x = sing->x; }
			if (sing->y < sing_min_y) { sing_min_y = sing->y; }
			if (sing->y > sing_max_y) { sing_max_y = sing->y; }
		}
		xx = xx / singularities.size();
		yy = yy / singularities.size();

		mean_point = new icVector3(xx, yy, singularities[0]->z);
	}

	void trim()
	{
		double x1, x2, y1, y2;
		std::vector<Quad*> temp_quads;
		min_x = sing_min_x; max_x = sing_max_x; min_y = sing_min_y; max_y = sing_max_y;
		for (Quad* quad : cell_quads)
		{
			x1 = std::min(std::min(quad->verts[0]->x, quad->verts[1]->x), quad->verts[2]->x);
			x2 = std::max(std::max(quad->verts[0]->x, quad->verts[1]->x), quad->verts[2]->x);
			y1 = std::min(std::min(quad->verts[0]->y, quad->verts[1]->y), quad->verts[2]->y);
			y2 = std::max(std::max(quad->verts[0]->y, quad->verts[1]->y), quad->verts[2]->y);

			if (x2 > sing_min_x && x1 < sing_max_x && y2 > sing_min_y && y1 < sing_max_y)
			{
				temp_quads.push_back(quad);
				// update min and max values
				if (x1 < min_x) { min_x = x1; }
				if (x2 > max_x) { max_x = x2; }
				if (y1 < min_y) { min_y = y1; }
				if (y2 > max_y) { max_y = y2; }
			}
		}
		cell_quads.clear();
		cell_quads = temp_quads;
	}
	
	void perimeter()
	{
		Vertex* v_temp;
		Edge* e_temp;
		// find interior and exterior verts and perimeter edges
		for (Quad* quad : cell_quads)
		{
			for (int i = 0; i < 4; i++)
			{
				v_temp = quad->verts[i];
				if (v_temp->x == min_x ||
					v_temp->x == max_x ||
					v_temp->y == min_y ||
					v_temp->y == max_y)
				{
					perim_verts.insert(v_temp);
				}

				e_temp = quad->edges[i];
				if ((e_temp->verts[0]->x == min_x && e_temp->verts[1]->x == min_x) ||
					(e_temp->verts[0]->x == max_x && e_temp->verts[1]->x == max_x) ||
					(e_temp->verts[0]->y == min_y && e_temp->verts[1]->y == min_y) ||
					(e_temp->verts[0]->y == max_y && e_temp->verts[1]->y == max_y))
				{
					perim_edges.insert(e_temp);
				}
			}
			quad->in_cluster = true;
		}
	}

	void separatrices()
	{
		// todo: find separatrices using polar linear interpolation
	}

public:

	// constructor
	Cluster(std::vector<icVector3*> points, std::vector<Quad*> quads)
	{
		singularities = points;
		cell_quads = quads;
		min_x = min_y = sing_min_x = sing_min_y = DBL_MAX;
		max_x = max_y = sing_max_x = sing_max_y = -DBL_MAX;

		perim_verts = std::unordered_set<Vertex*>();
		perim_edges = std::unordered_set<Edge*>();
		mean_point = nullptr;

		// find min and max values
		Vertex* v_temp;
		for (Quad* quad : cell_quads)
		{
			for (int i = 0; i < 4; i++)
			{
				v_temp = quad->verts[i];
				if (v_temp->x < min_x) { min_x = v_temp->x; }
				if (v_temp->x > max_x) { max_x = v_temp->x; }
				if (v_temp->y < min_y) { min_y = v_temp->y; }
				if (v_temp->y > max_y) { max_y = v_temp->y; }
			}
		}
	}
	
	// prepare cluster for analysis
	void finalize()
	{
		find_mean();
		trim();
		perimeter();
		separatrices();
	}

	// getters

	std::vector<icVector3*> getSingularities() {
		return singularities;
	}

	std::vector<Quad*> getCell() {
		return cell_quads;
	}

	std::unordered_set<Edge*> getPerim()
	{
		return perim_edges;
	}

	double getMinx() { return min_x; }
	double getMiny() { return min_y; }
	double getMaxx() { return max_x; }
	double getMaxy() { return max_y; }

	double getSingMinx() { return sing_min_x; }
	double getSingMiny() { return sing_min_y; }
	double getSingMaxx() { return sing_max_x; }
	double getSingMaxy() { return sing_max_y; }
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
		Vertex* center = nullptr;
		double cen_x = (c->getMinx() + c->getMaxx()) / 2;
		double cen_y = (c->getMiny() + c->getMaxy()) / 2;

		// find centermost vertex
		double dist;
		double min_dist = DBL_MAX;
		Vertex* v_temp;
		for (Quad* quad : c->getCell())
		{
			for (int i = 0; i < 4; i++)
			{
				v_temp = quad->verts[i];
				dist = sqrt(pow(cen_x - v_temp->x, 2) + pow(cen_y - v_temp->y, 2));
				if (dist < min_dist)
				{
					min_dist = dist;
					center = v_temp;
				}
			}
		}

		// compute center vertex approximation error
		double sum = 0.0;
		double weights = 0.0;
		for (icVector3* sing : c->getSingularities())
		{
			sum += sqrt(pow(center->x - sing->x, 2) + pow(center->y - sing->y, 2));
			weights += 1;
		}
		s = sum / weights;

		// either recurse or add cluster to vector
		if (s > cluster_threshold && c->getCell().size() > 2) 
		{
			//step 3: compute coordinate axis variance
			double variance_x = 0.0, variance_y = 0.0;
		
			for (icVector3* s : c->getSingularities()) {
				variance_x += pow(center->x - s->x, 2);
				variance_y += pow(center->y - s->y, 2);
			}
			
			std::vector<icVector3*> leftSings = std::vector<icVector3*>();
			std::vector<Quad*> leftCell = std::vector<Quad*>();
			std::vector<icVector3*> rightSings = std::vector<icVector3*>();
			std::vector<Quad*> rightCell = std::vector<Quad*>();

			if (variance_x > variance_y) {
				std::vector<icVector3*> clusterSings = c->getSingularities();
				for (icVector3* s : clusterSings) {
					if (s->x > center->x) {
						rightSings.push_back(s);
					}
					else {
						leftSings.push_back(s);
					}
				}

				std::vector<Quad*> clusterCell = c->getCell();
				for (Quad* q : clusterCell) {
					double avg = (q->verts[0]->x + q->verts[1]->x + q->verts[2]->x + q->verts[3]->x) / 4;
					if (avg > center->x) {
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
					if (s->y > center->y) {
						rightSings.push_back(s);
					}
					else {
						leftSings.push_back(s);
					}
				}

				std::vector<Quad*> clusterCell = c->getCell();
				for (Quad* q : clusterCell) {
					double avg = (q->verts[0]->y + q->verts[1]->y + q->verts[2]->y + q->verts[3]->y) / 4;
					if (avg > center->y) {
						rightCell.push_back(q);
					}
					else {
						leftCell.push_back(q);
					}
				}
			}
			
			// recursive calls
			std::vector<Cluster*> join = std::vector<Cluster*>();
			if (leftSings.size() > 1)
			{
				std::vector<Cluster*> left = separateClusters(new Cluster(leftSings, leftCell));
				join.reserve(left.size());
				join.insert(join.end(), left.begin(), left.end());
			}
			if (rightSings.size() > 1)
			{
				std::vector<Cluster*> right = separateClusters(new Cluster(rightSings, rightCell));
				join.reserve(right.size());
				join.insert(join.end(), right.begin(), right.end());
			}
			return join;
		}
		else {
			std::vector<Cluster*> ret = std::vector<Cluster*>();
			c->finalize();
			ret.push_back(c);
			return ret;
		}
	}

	std::vector<Cluster*> getClusters() {
		return clusters;
	}

};