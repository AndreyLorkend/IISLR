#include <iostream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include <cstdlib>

const double Lferm = 5.71;
const double K = 1.0;
const double Fi = 4.5;
const int N = 10000;
double Random[N][2];

using namespace std;

class KoordPoint
{   public:
	double x;
	double y;
};
//KoordPoint koordpointxy[12] =
//{ {0,  0},{ 10, 0},{ 14, 4},{ 27, 4},
//{ 20, 7},{ 25, 7},{ 25, 0},{ 27, 0},
//{ 27, 4},{ 29, 4},{ 29, 8},{ 39, 8}};

KoordPoint koordpointxy[8] = { 
	{0,  0},{ 10, 0},{ 14, 4},{ 25, 4},
    { 25, 0},{ 32, 0},{ 38, 8},{ 48, 8}
};

double J(double x, double y, double z, double U)
{
	double S1 = 3.0 / (K * Fi);
	double Zxy = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); 
	double S2 = Zxy*(1.0 - 23.0 / (3 * Fi * K * Zxy + 10 - 2 * U * K * Zxy)) + S1;
	double FiOfZ = Fi - (U*(S1 + S2) / (2 * Zxy)) - (2.86 / (K*(S2 - S1))*log((S2*(Zxy - S1)) / (S1*(Zxy - S2))));
	double res = 1620.0*U*Lferm*exp(-1.025*Zxy*sqrt(FiOfZ));
	return res;
}

double MC_Int(double xL, double xH, double yL, double yH, double Z, double U)
{
	double b[2] = { xH, yH }; //Верхние пределы интегрирования
	double a[2] = { xL, yL }; //Нижние пределы интегрирования
	double result = 0.0; //Результат интегрирования
	double V = (b[0] - a[0])*(b[1] - a[1]); //Объем
	double X[2] = { 0.0, 0.0 }; //Буферный массив
	bool   Flag = true; //Флаг, определяющий попали ли мы в интервал

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 2; j++) X[j] = a[j] + (b[j] - a[j])*Random[i][j];
		Flag = true;
		for (int j = 0; j < 2; j++)
			if ((X[j] < a[j]) && (X[j]>b[j])) Flag = false;
		if (Flag) result += J(X[0], X[1], Z, U);
	}

	return result * V / N;
}

/* расчет нового значения туннельного тока */
double NewJ(double xp, double zp, double U)
{
	const int MAX = 100;

	double xH = 0.0, xL = 0.0, yH = 10.0, yL = -10.0;
	double Z = 0.0;
	double j = 0.0;
	bool Vid[7]; //Массив видимости

	for (int i = 0; i < 7; i++)
	{
		Vid[i] = false;
		if ((abs(koordpointxy[i].x - xp) <= MAX) && ((koordpointxy[i].y - zp) <= MAX))
		Vid[i] = true;
	}
	if (Vid[0] && (zp > 0))
	{
		xL = koordpointxy[0].x - xp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[1].x - xp;
		if (xH > MAX) xH = MAX;
		Z = zp;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[1] && (zp > (xp - koordpointxy[1].x)))
	{
		Z = (zp - (xp - koordpointxy[1].x))*cos(M_PI / 4);
		xL = (koordpointxy[1].x - (xp + Z*cos(M_PI / 4))) / cos(M_PI / 4);
		if (xL > MAX) xL = MAX;
		xH = (koordpointxy[2].x - (xp + Z*cos(M_PI / 4))) / cos(M_PI / 4);
		if (xH > MAX) xH = MAX;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[2] && (zp > koordpointxy[2].y))
	{
		xL = koordpointxy[2].x - xp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[3].x - xp;
		if (xH > MAX) xH = MAX;
		Z = zp - koordpointxy[2].y;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[3] && (xp > koordpointxy[3].x))
	{
		
		xL = koordpointxy[3].y - zp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[4].y - zp;
		if (xH > MAX) xH = MAX;
		Z = koordpointxy[3].x - xp;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[4] && (zp > koordpointxy[4].y))
	{
		xL = koordpointxy[4].x - xp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[5].x - xp;
		if (xH > MAX) xH = MAX;
		Z = zp - koordpointxy[4].y;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[5] && (zp > (xp - koordpointxy[5].x)))
	{
		Z = (zp - (xp - koordpointxy[5].x)) * cos(M_PI / 4);
		xL = (koordpointxy[5].x - (xp + Z * cos(M_PI / 4))) / cos(M_PI / 4);
		if (xL > MAX) xL = MAX;
		xH = (koordpointxy[6].x - (xp + Z * cos(M_PI / 4))) / cos(M_PI / 4);
		if (xH > MAX) xH = MAX;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[6] && (zp > koordpointxy[6].y))
	{
		xL = koordpointxy[6].x - xp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[7].x - xp;
		if (xH > MAX) xH = MAX;
		Z = zp - koordpointxy[7].y;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	//if (Vid[7] && (xp < koordpointxy[7].x))
	//{
	//	xL = koordpointxy[7].y - zp;
	//	if (xL > MAX) xL = MAX;
	//	xH = koordpointxy[8].y - zp;
	//	if (xH > MAX) xH = MAX;
	//	Z = koordpointxy[7].x - xp;
	//	j += MC_Int(xL, xH, yL, yH, Z, U);
	//}
	//if (Vid[8] && (zp > koordpointxy[8].y))
	//{
	//	xL = koordpointxy[8].x - xp;
	//	if (xL > MAX) xL = MAX;
	//	xH = koordpointxy[9].x - xp;
	//	if (xH > MAX) xH = MAX;
	//	Z = zp - koordpointxy[8].y;
	//	j += MC_Int(xL, xH, yL, yH, Z, U);
	//}
	/*if (Vid[9] && (xp < koordpointxy[9].x))
	{
		xL = koordpointxy[9].y - zp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[10].y - zp;
		if (xH > MAX) xH = MAX;
		Z = koordpointxy[9].x - xp;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}
	if (Vid[10] && (zp > koordpointxy[10].y))
	{
		xL = koordpointxy[10].x - xp;
		if (xL > MAX) xL = MAX;
		xH = koordpointxy[11].x - xp;
		if (xH > MAX) xH = MAX;
		Z = zp - koordpointxy[11].y;
		j += MC_Int(xL, xH, yL, yH, Z, U);
	}*/
	return j;
}

void BuildProf(double Z0, double U0, HDC hDC)
{
	const double e = 0.05;
	const double step = 0.1;

	double   j0 = 0.0;
	double   z[3] = { 0.0, 0.0, 0.0 };
	double   j[3] = { 0.0, 0.0, 0.0 };
	KoordPoint Diagramm[450];
	POINT    dot;

	Diagramm[0].x = 0.0;
	Diagramm[0].y = Z0;
	MoveToEx(hDC, 5, 280 - Diagramm[0].y * 11, &dot);
	for (int i = 0; i < 10; i++) j0 += NewJ(0, Z0, U0) / 10.0;
	for (int i = 1; i < 450; i++)
	{
		Diagramm[i].x = Diagramm[i - 1].x + step;
		z[2] = Diagramm[i - 1].y;
		j[1] = NewJ(Diagramm[i].x, z[2], U0);
		if ((j[1] < (j0 * (1 - e))) || (j[1] > (j0 * (1 + e))))
		{
			z[1] = z[2] + step * ((j[1] - j0) / abs(j[1] - j0));
			j[0] = NewJ(Diagramm[i].x, z[1], U0);
			if ((j[0] < (j0 * (1 - e))) || (j[0] > (j0 * (1 + e))))
			{
				int k = 0;
				do
				{
					z[0] = z[1] - (j[0] - j0) * (z[1] - z[2]) / (j[0] - j[1]);
					j[1] = j[0];
					j[0] = NewJ(Diagramm[i].x, z[0], U0);
					z[2] = z[1];
					z[1] = z[0];
					k++;
					if (k > 1000) break;
				} while (!(j[0] >= j0 * (1 - e) && (j[0] <= j0 * (1 + e))));
				Diagramm[i].y = z[0];
				if (k > 1000) Diagramm[i].y = Diagramm[i - 1].y;
			}
			else Diagramm[i].y = z[1];
		}
		else Diagramm[i].y = z[2];

		LineTo(hDC, Diagramm[i].x * 11 + 5, 280 - Diagramm[i].y * 11);
	}
}


int main()
{
	setlocale(LC_ALL, "Russian");
	int k = 0;
	string s = "";
	double res = J(2.0, 2.0, 2.0, 0.01);
	cout << res << endl;

	{ for (int i = 0; i < N; i++)
		for (int j = 0; j < 2; j++) Random[i][j] = (double)rand() / RAND_MAX;

	cout << MC_Int(-5, 5, -10, 10, 100.0, 0.01) << endl;

	cout << "Введите:" << endl;
	cout << "1 - Для вывода профилограмм (U = 0.01 В)" << endl;
	cout << "2 - Для вывода профилограмм (U = 0.1 В)" << endl;
	cout << "3 - Для выхода" << endl;
	cin >> k;
	

	/* Инициализация рисования */
	HWND  hWnd = GetConsoleWindow();
	HDC   hDC = GetDC(hWnd);
	POINT dot;
	RECT  consoleRect;
	SelectObject(hDC, GetStockObject(WHITE_PEN));
	GetClientRect(hWnd, &consoleRect);

	double Z0[4] = { 5, 7, 10, 15 };
	//for (int j = 0; j < 2; j++)
	system("cls");
	FillRect(hDC, &consoleRect, (HBRUSH)GetStockObject(BLACK_BRUSH));
	
	MoveToEx(hDC, 5, 280, &dot);
	for (int i = 1; i < 8; i++) LineTo(hDC, koordpointxy[i].x * 11 + 5, 280 - koordpointxy[i].y * 11);
	switch (k)
	{
	case 1:
		cout << "U0 = 0.01 В" << endl;
		for (int i = 0; i < 4; i++)
		{
			BuildProf(Z0[i], 0.01, hDC);
		}
		system("pause");
		return 0;
	case 2:
		cout << "U0 = 0.1 В" << endl;
		for (int i = 0; i < 4; i++)
		{
			BuildProf(Z0[i], 0.1, hDC);
		}
		system("pause");
		return 0;
	case 3:
		system("pause");
		return 0;
	}}
}