#include <vector>
#include <tuple>

#ifndef COLOR_H
#define COLOR_H

std::vector<std::tuple<double, double, double>> colors = {
	std::make_tuple<double, double, double>(0,0,0),
	std::make_tuple<double, double, double>(30,0,0),
	std::make_tuple<double, double, double>(60,0,0),
	std::make_tuple<double, double, double>(90,0,0),
	std::make_tuple<double, double, double>(120,0,0),
	std::make_tuple<double, double, double>(150,0,0),
	std::make_tuple<double, double, double>(180,0,0),
	std::make_tuple<double, double, double>(210,0,0),
	std::make_tuple<double, double, double>(240,0,0),
	std::make_tuple<double, double, double>(255,15,0),
	std::make_tuple<double, double, double>(240,30,0),
	std::make_tuple<double, double, double>(210,60,0),
	std::make_tuple<double, double, double>(180,90,0),
	std::make_tuple<double, double, double>(150,120,0),
	std::make_tuple<double, double, double>(120,150,0),
	std::make_tuple<double, double, double>(90,180,0),
	std::make_tuple<double, double, double>(60,210,0),
	std::make_tuple<double, double, double>(30,240,0),
	std::make_tuple<double, double, double>(0,255,0),
	std::make_tuple<double, double, double>(0,240,0),
	std::make_tuple<double, double, double>(0,210,0),
	std::make_tuple<double, double, double>(0,180,0),
	std::make_tuple<double, double, double>(0,150,0),
	std::make_tuple<double, double, double>(0,120,0),
	std::make_tuple<double, double, double>(0,90,0),
	std::make_tuple<double, double, double>(0,60,0),
	std::make_tuple<double, double, double>(0,30,0),
	std::make_tuple<double, double, double>(0,0,30),
	std::make_tuple<double, double, double>(0,0,60),
	std::make_tuple<double, double, double>(0,0,90),
	std::make_tuple<double, double, double>(0,0,120),
	std::make_tuple<double, double, double>(0,0,150),
	std::make_tuple<double, double, double>(0,0,180),
	std::make_tuple<double, double, double>(0,0,210),
	std::make_tuple<double, double, double>(0,0,240),
	std::make_tuple<double, double, double>(15,0,255),
	std::make_tuple<double, double, double>(30,0,240),
	std::make_tuple<double, double, double>(60,0,210),
	std::make_tuple<double, double, double>(90,0,180),
	std::make_tuple<double, double, double>(120,0,150),
	std::make_tuple<double, double, double>(150,0,120),
	std::make_tuple<double, double, double>(180,0,90),
	std::make_tuple<double, double, double>(210,0,60),
	std::make_tuple<double, double, double>(240,0,30),
	std::make_tuple<double, double, double>(30,30,30),
	std::make_tuple<double, double, double>(60,60,60),
	std::make_tuple<double, double, double>(90,90,90),
	std::make_tuple<double, double, double>(120,120,120),
	std::make_tuple<double, double, double>(150,150,150),
	std::make_tuple<double, double, double>(180,180,180),
	std::make_tuple<double, double, double>(210,210,210),
	std::make_tuple<double, double, double>(240,240,240)
};

double getColor(int index, int type) {
	if (type == 0) {
		return (std::get<0>(colors[index])) / 255;
	}
	else if (type == 1) {
		return (std::get<1>(colors[index])) / 255;
	}
	else if (type == 2) {
		return (std::get<2>(colors[index])) / 255;
	}
	else {
		return NULL;
	}
}

#endif

