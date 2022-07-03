#include "model.cc"

//using namespace std;

void debug()
{
	model* m = new model();
	
	// set hodo ids
	//41, 42, 43, 44, 45, 46
	m->set_det(41);
	m->set_det(42);
	m->set_det(43);
	m->set_det(44);
  m->set_det(45);
  m->set_det(46);

	m->ana();
	m->done();
}
