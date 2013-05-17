#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <gmp.h>

using namespace std;

class EllipticCurve;
EllipticCurve init();

class Point
{
public:
	mpq_t x;
	mpq_t y;
	Point()
	{
		mpq_init(x);
		mpq_init(y);
	}
	~Point()
	{
		mpq_clear(x);
		mpq_clear(y);
	}
};

class EllipticCurve
{
public:
	int a2;
	int a4;
	int a6;
	int genx;
	int geny;
	EllipticCurve(int pa2, int pa4, int pa6, int gx, int gy)
	{
		a2=pa2;
		a4=pa4;
		a6=pa6;
		genx=gx;
		geny=gy;
	}
	EllipticCurve()
	{
		a2=a4=a6=0;
	}
};

int main()
{
	srand(time(NULL));
	EllipticCurve E=init();
	cout << E.a2 << " " << E.a4 << " " << E.a6 << " " << E.genx << " " << E.geny << endl;
}

EllipticCurve init()
{
	ifstream inf("processed_Ecurves.txt");
	vector<string> allCurves;
	while(inf)
	{
		allCurves.push_back("");
		getline(inf,allCurves[allCurves.size()-1]);
	}
	int rnd=rand()%allCurves.size();
	EllipticCurve E;

	int i=1;
	if(allCurves[rnd][i]=='-') //always -1,0,1
		E.a2=-1;
	else if(allCurves[rnd][i]=='1')
		E.a2=1;
	else if(allCurves[rnd][i]=='0')
		E.a2=0;

	while(allCurves[rnd][i]!=',')
		i++;
	i++;
	//now at second number
	int num=0;
	bool negative=false;
	while(allCurves[rnd][i]!=',') //while still in the number
	{
		if(allCurves[rnd][i]=='-')
		{
			negative=true;
			i++;
			continue;
		}
		num*=10;
		num+=(allCurves[rnd][i]-'0');
		i++;
	}
	if(negative==true)
		num*=-1;
	E.a4=num;

	i++; //get to third number
	num=0;
	negative=false;
	while(allCurves[rnd][i]!=']') //while still in the number
	{
		if(allCurves[rnd][i]=='-')
		{
			negative=true;
			i++;
			continue;
		}
		num*=10;
		num+=(allCurves[rnd][i]-'0');
		i++;
	}
	if(negative==true)
		num*=-1;
	E.a6=num;
	i+=3; //get to coordinates

	num=0;
	negative=false;
	while(allCurves[rnd][i]!=':') //while still in the number
	{
		if(allCurves[rnd][i]=='-')
		{
			negative=true;
			i++;
			continue;
		}
		num*=10;
		num+=(allCurves[rnd][i]-'0');
		i++;
	}
	if(negative==true)
		num*=-1;
	E.genx=num;
	i++; //get to final number

	num=0;
	negative=false;
	while(allCurves[rnd][i]!=']') //while still in the number
	{
		if(allCurves[rnd][i]=='-')
		{
			negative=true;
			i++;
			continue;
		}
		num*=10;
		num+=(allCurves[rnd][i]-'0');
		i++;
	}
	if(negative==true)
		num*=-1;
	E.geny=num;
	return E;
}

void addToSelf(EllipticCurve E, Point gen)
{
	mpq_t slope;
	mpq_init(slope);
	mpq_t p1,p2;
	mpq_init(p1);
	mpq_init(p2);

	mpz_mul_si(mpq_numref(p1),mpq_numref(gen.x),E.a2); //a*x1 == first coef*numerator of x
	mpz_mul(mpq_numref(p1),mpq_numref(p1),mpq_denref(gen.y)); //ans*y2 == ans * denominator of y

	mpz_mul(mpq_denref(p1),mpq_numref(gen.y),mpq_denref(gen.x));

	mpq_canonicalize(p1);

	mpz_mul_si(mpq_numref(p1),mpq_denref(gen.y),E.a4);
	mpz_mul_si(mpq_denref(p1),mpq_numref(gen.y),2);

	mpq_canonicalize(p2);

	mpq_add(slope,p1,p2);
	mpq_clear(p1);
	mpq_clear(p2);
	//slope is now correct
}
