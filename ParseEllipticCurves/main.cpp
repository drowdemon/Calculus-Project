#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

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
	ifstream inf("elliptic_curves.txt");
	ofstream outf("processed_Ecurves.txt");
	string input;
	vector<EllipticCurve> allCurves;
	if(!inf)
		cout << "No File found" << endl;
	while(inf)
	{
		getline(inf, input);
		int numspaces=0;
		EllipticCurve E;
		unsigned int i=0;
		for(i=0; i<input.size(); i++)
		{
			if(input[i]==' ')
				numspaces++;
			if(numspaces<3)
				continue;
			if(numspaces==3) //at third space
			{
				break;
			}
		}

		i+=2; //move up 2: at first number, always 1 or 0
		numspaces++;
		if(input[i]!='0') //We want a1=0
			continue; //next curve
		else
			i+=2;//next number
		if(input[i]=='-') //always -1,0,1
			E.a2=-1;
		else if(input[i]=='1')
			E.a2=1;
		else if(input[i]=='0')
			E.a2=0;
		while(input[i]!=',')
			i++;
		//now at comma

		i++; //get to the actual number
		if(input[i]!='0') //we want it to be 0
			continue; //next curve
		i+=2; //next number. number 4
		int num=0;
		bool negative=false;
		while(input[i]!=',') //while still in the number
		{
			if(input[i]=='-')
			{
				negative=true;
				i++;
				continue;
			}
			num*=10;
			num+=(input[i]-'0');
			i++;
		}
		if(negative==true)
			num*=-1;
		E.a4=num; //at comma for last number

		i++;
		num=0; //last number
		negative=false;
		while(input[i]!=']') //while still in the number
		{
			if(input[i]=='-')
			{
				negative=true;
				i++;
				continue;
			}
			num*=10;
			num+=(input[i]-'0');
			i++;
		}
		E.a6=num; //currently at bracket
		i+=2; //at next number: rank of elliptic curve
		if(input[i]=='0') //rank 0
			continue; //next curve
		int numbrackets=0;
		while(numbrackets<2)
		{
			if(input[i]=='[')
				numbrackets++;
			i++;
		}
		//at beginning of number after second bracket
		num=0; //coo of generator
		negative=false;
		while(input[i]!=':') //while still in the number
		{
			if(input[i]=='-')
			{
				negative=true;
				i++;
				continue;
			}
			num*=10;
			num+=(input[i]-'0');
			i++;
		}
		E.genx=num;
		i++; //at beginning of number, after colon

		num=0; //coo of second generator
		negative=false;
		while(input[i]!=':') //while still in the number
		{
			if(input[i]=='-')
			{
				negative=true;
				i++;
				continue;
			}
			num*=10;
			num+=(input[i]-'0');
			i++;
		}
		E.geny=num;
		i++; //skip colon
		if(input[i]!='1')
			continue; //next curve
		else
		{
			allCurves.push_back(E);
			outf << "[" << E.a2 << "," << E.a4 << "," << E.a6 << "] [" << E.genx << ":" << E.geny << "]\n";
			//cout << "y^2 = x^3 + " << E.a2 << "x^2 + " << E.a4 << "x + " << E.a6 << endl;
		}
	}
	cout << "Done" << endl;
	return 0;
}
