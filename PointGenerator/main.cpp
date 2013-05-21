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

class EllipticCurve //y^2=x^3+a2x^2+a4x+a6
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

Point addToSelf(EllipticCurve E, Point gen) // returns 2p.
{
	mpq_t slope;
	mpq_init(slope);
	mpq_t p1,p2;
	mpq_init(p1);
	mpq_init(p2);

	mpz_mul_si(mpq_numref(p1),mpq_numref(gen.x),3); //3x_1
	mpz_mul(mpq_numref(p1),mpq_numref(p1),mpq_numref(gen.x)); //3x_1^2
	mpz_mul(mpq_numref(p1),mpq_numref(p1),mpq_denref(gen.y)); //3x_1^2y_2 == ans*y2 == ans * denominator of y

	mpz_mul(mpq_denref(p1),mpq_numref(gen.y),mpq_denref(gen.x)); //x_2y_1
	mpz_mul_si(mpq_denref(p1),mpq_denref(p1),2); //2*x_2y_1

	mpq_canonicalize(p1); //first fraction complete

	mpz_mul_si(mpq_numref(p1),mpq_denref(gen.y),E.a4); //Ay_2
	mpz_mul_si(mpq_denref(p1),mpq_numref(gen.y),2); //2y_1

	mpq_canonicalize(p2); //second fraction complete

	mpq_add(slope,p1,p2); //add
	//mpq_clear(p1);
	//mpq_clear(p2);
	//slope is now correct //Re-using variables because allocation is slow.
	mpq_set(p1,slope); //sets p1=slope
	mpq_set(p2,gen.x); //sets p2=gen.x
	mpq_mul(p1,slope,slope); //p1=slope^2
	mpz_mul_si(mpq_numref(p2),mpq_numref(gen.x),2); //changing the given point: multiplication by 2: alpha=2alpha
	mpq_canonicalize(p2);
	mpq_sub(p1,p1,p2); //p1=p1-p2: p1=slope^2-2alpha
	//slope is unchanged
	//p2 is useless, again. p1 contains the correct x value - to be returned.
	mpq_sub(p2,p1,gen.x); //p2=p1-x ... x-alpha
	mpq_mul(p2,p2,slope); //p2=slope*p2... p2=slope(p1-alpha)
	mpq_add(p2,p2,gen.y); //p2=p2+gen.y ... p2=slope(x-alpha)+beta == y coordinate.
	//curve is symmetric about the y axis. We want the other point.
	mpq_neg(p2,p2);

	Point ret;
	mpq_set(ret.x,p1);
	mpq_set(ret.y,p2);
	mpq_clear(slope);
	mpq_clear(p1);
	mpq_clear(p2);
	return ret;
}

Point add(EllipticCurve E, Point gen1, Point gen2) // returns gen1+gen2
{
	mpq_t slope;
	mpq_init(slope);
	mpq_t p1;
	mpq_init(p1);

	mpq_sub(slope,gen2.y,gen1.y); //numerator of slope: slope=y2-y1
	mpq_sub(p1,gen2.x,gen1.x); //denominator of slope: p1=x2-x1
	mpq_div(slope,slope,p1); //slope=slope/p1, proper slope
	//mpq_clear(p1);
	//slope is now correct //Re-using variables because allocation is slow.
	mpq_set(p1,slope); //sets p1=slope
	mpq_mul(p1,slope,slope); //p1=slope^2
	mpq_sub(p1,p1,gen1.x);
	mpq_sub(p1,p1,gen2.x); //p1=slope^2-alpha_1-alpha_2
	//net effect right now, from original slope, p1=slope^2-2alpha - correct x value
	//p2 is useless, again. p1 contains the correct x value - to be returned.
	Point ret;
	mpq_set(ret.x,p1); //p1 stored - x val stored
	mpq_sub(p1,p1,gen1.x); //p1=p1-x ... x-alpha
	mpq_mul(p1,p1,slope); //p1=slope*p1... p1=slope(x-alpha)
	mpq_add(p1,p1,gen1.y); //p2=p2+gen.y ... p1=slope(x-alpha)+beta == y coordinate.
	//curve is symmetric about the y axis. We want the other point.
	mpq_neg(p1,p1);

	mpq_set(ret.y,p1);
	mpq_clear(slope);
	mpq_clear(p1);
	return ret;
}
