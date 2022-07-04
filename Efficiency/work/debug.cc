#include "model.cc"
#include <vector>
#include <iostream>

using namespace std;

void debug()
{	
	// set hodo ids
	vector<int> hodoid = {41, 42, 43, 44, 45, 46};
	int nhodo = hodoid.size();

	//cout << nhodo << endl;

	// do ana
	model* m = new model();	
	for(int i = 0; i < nhodo; i++)
	{
		int id = hodoid.at(i);
		m->set_det(id);
	}
	m->ana();
	//m->done();
}
