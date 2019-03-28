#ifndef SMMAP_H
#define SMMAP_H
#include <gmapping/grid/map.h>
#include <gmapping/grid/harray2d.h>
#include <gmapping/utils/point.h>
#include <iostream>
#include <vector>
#define SIGHT_INC 1

namespace GMapping {

struct PointAccumulator{
	typedef point<float> FloatPoint;
	FloatPoint acc;
	int n, visits, num_glass_detected;
	inline double entropy() const;
	bool glass_detected;
	/* before
	PointAccumulator(int i=-1): acc(0,0), n(0), visits(0){assert(i==-1);}
	*/
	/*after begin*/
	PointAccumulator(): acc(0,0), n(0), visits(0), glass_detected(false), num_glass_detected(0) {}
	PointAccumulator(int i): acc(0,0), n(0), visits(0), glass_detected(false), num_glass_detected(0) {assert(i==-1);}
	/*after end*/
	inline operator double() const { return visits?(double)n*SIGHT_INC/(double)visits:-1; }

    inline void update(bool value, const Point& p=Point(0,0));
	inline Point mean() const {return 1./n*Point(acc.x, acc.y);}
	inline void add(const PointAccumulator& p) {acc=acc+p.acc; n+=p.n; visits+=p.visits; }
	static const PointAccumulator& Unknown();
	static PointAccumulator* unknown_ptr;
	bool isGlassDetected() { return glass_detected; }
	int getGlassNum() { return num_glass_detected; }
	bool resetGlass() { glass_detected = false; num_glass_detected = 0; }
	void updateGlass() { glass_detected=true; }
	void incGlass() { num_glass_detected++; }
};

void PointAccumulator::update(bool value, const Point& p){
	if (value) {
		acc.x+= static_cast<float>(p.x);
		acc.y+= static_cast<float>(p.y); 
		n++; 
		visits+=SIGHT_INC;
	} else
		visits++;
}

double PointAccumulator::entropy() const{
	if (!visits)
		return -log(.5);
	if (n==visits || n==0)
		return 0;
	double x=(double)n*SIGHT_INC/(double)visits;
	return -( x*log(x)+ (1-x)*log(1-x) );
}


typedef Map<PointAccumulator,HierarchicalArray2D<PointAccumulator> > ScanMatcherMap;

};

#endif 
