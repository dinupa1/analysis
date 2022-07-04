#include "model.cc"
#include <vector>
#include <iostream>

using namespace std;

void debug()
{	
	// set hodo ids
	vector<int> hodo = {41, 42, 43, 44, 45, 46};
	int nhodo = hodo.size();

	//cout << nhodo << endl;

	// do ana
	model* m = new model();	
	for(int i = 0; i < nhodo; i++)
	{
		m->set_det(hodo.at(i));
		cout << hodo.at(i) << endl;
	}
	m->ana();
	m->done();
}
