#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <cmath>

using namespace std;

//ESTRUCTURAS Y VARIABLES GLOBALES---------------------------------------------------

// Estructura para guardar los puntos primos del segundo problema
struct {
	float pPrimoX;
	float pPrimoY;
	float pX;
	float pY;
	float pZ;
}punto[10];
//Estructuras para superficies visibles
struct SPunto {
	float x, y, z;
};
struct Cara {
	vector<SPunto>  puntos;
};
// Definición de la estructura Punto en 3D
struct Punto {
	double x, y, z;
};
// Definición de la estructura para cuaterniones
struct Quaternion {
	double w, x, y, z;
};
// Definición de la función conjugate
Quaternion conjugate(const Quaternion& q) {
	Quaternion result;
	result.w = q.w;
	result.x = -q.x;
	result.y = -q.y;
	result.z = -q.z;
	return result;
}
//Constante PI
const double M_PI = 3.14159265358979323846;

//---------------------------------------------------------------------------------

//DECLARACION DE FUNCIONES-----------------------------------------------------------

void llenarMatriz(int** M, int fil, int col);
void GirarPuntoVector();
void GirarFiguraVector();
void sumarMatrices();
void restarMatrices();
void multiplicarMatrices();
void obtenerMatrizInversa();
void obtenerMatrizTranspuesta();
void obtenerMatrizAdjunta();
void rotarMatriz(float** M, int fil, int col, char eje, int angulo);
void trasladarMatriz();
void obtenerPuntoTranformaciones(int dx, int dy, int dz);
void llenarMatrizF(float** M, int fil, int col);
void otraFunc();
void imprimirMatriz(int** M, int fil, int col);
void imprimirMatrizInversa(float** M, int fil, int col);
void limpiarPantalla();
void obtenerPunto();
void mostrarMenu();
void trazoBresenham();
void reglaTrazo(int, int, int, int);
void curvasBresenham();
void reglaCurvas(int xc, int yc, int r);
void elipseBresenham();
void reglaElipse(int xc, int yc, int rx, int ry);
Quaternion conjugate(const Quaternion& q);
Quaternion normalizeQuaternion(const Quaternion& q);
Quaternion multiplyQuaternions(const Quaternion& q1, const Quaternion& q2);

//---------------------------------------------------------------------------------

//PRINCIPAL-------------------------------------------------------------------------

int main()
{
	int opcion;

	do {
		mostrarMenu();
		cin >> opcion;
		if (cin.fail() || cin.bad()) {
			opcion = 0;
			cin.clear();
			cin.ignore('\n', 80);
		}

		switch (opcion) {
		case 1:
			sumarMatrices();
			break;

		case 2:
			restarMatrices();
			break;

		case 3:
			multiplicarMatrices();
			break;

		case 4:
			obtenerMatrizInversa();
			break;

		case 5:
			obtenerMatrizTranspuesta();
			break;

		case 6:
			obtenerMatrizAdjunta();
			break;

		case 7:
			obtenerPunto();
			break;

		case 8:
			GirarPuntoVector();
			break;

		case 9:
			GirarFiguraVector();
			break;

		case 10:
			otraFunc();
			break;

		case 11:
			trasladarMatriz();
			break;

		case 12:
			elipseBresenham();
			break;

		case 13:
			trazoBresenham();
			break;

		case 14:
			curvasBresenham();
			break;

		case 15:
			cout << "\n\nPROGRAMA FINALIZADO\n\n";
			exit(0);
			break;

		default:
			cout << "\n\nOPCION NO VALIDA\n\n";
			break;
		}

	} while (opcion != 15);

	return 0;
}

//---------------------------------------------------------------------------------

//OTROS----------------------------------------------------------------------------

//limpiar pantalla
void limpiarPantalla() {
	cout << "\n\nPresione Enter para continuar...";
	getchar();
	getchar();
	system("cls");
}
//Llenar matriz con los datos ingresados (solo enteros)
void llenarMatriz(int** M, int fil, int col)
{
	cout << "\nRellenar la matriz:\n";
	for (int i = 0; i < fil; i++) {
		for (int j = 0; j < col; j++) {
			cout << "Elemento [" << i << "][" << j << "]: ";
			cin >> M[i][j];
		}
	}
}
//Llenar matriz con los datos ingresados (solo flotantes)
void llenarMatrizF(float** M, int fil, int col)
{
	cout << "\nRellenar la matriz:\n";
	for (int i = 0; i < fil; i++) {
		for (int j = 0; j < col; j++) {
			cout << "Elemento [" << i << "][" << j << "]: ";
			cin >> M[i][j];
		}
	}
}
//Imprimir matric (enteros)
void imprimirMatriz(int** M, int fil, int col)
{
	for (int i = 0; i < fil; i++) {
		cout << "| ";
		for (int j = 0; j < col; j++) {
			cout << setw(4) << M[i][j] << " ";  // Ajusta el ancho del campo a 4 caracteres
		}
		cout << "|\n";
	}
	cout << endl;
}
//Imprimir matric (flotantes)
void imprimirMatrizInversa(float** M, int fil, int col)
{
	//cout << fixed << setprecision(4);
	for (int i = 0; i < fil; i++) {
		cout << "| ";
		for (int j = 0; j < col; j++) {
			cout << setw(9) << M[i][j] << " ";  // Ajusta el ancho del campo a 4 caracteres
		}
		cout << "|\n";
	}
	cout << endl;
}

//---------------------------------------------------------------------------------

//MENU-----------------------------------------------------------------------------

void mostrarMenu()
{
	cout << "\n>>> CALCULADORA DE MATRICES <<<"
		<< "\n-------------------------------"
		<< "\n[1] Sumar Matrices"
		<< "\n[2] Restar Matrices"
		<< "\n[3] Multiplicar Matrices"
		<< "\n[4] Obtener Matriz Inversa"
		<< "\n[5] Obtener Transpuesta"
		<< "\n[6] Obtener Adjunta"
		<< "\n[7] Obtener punto primo despues de realizar operaciones"
		<< "\n[8] Girar un punto por un vector"
		<< "\n[9] Girar una figura por un vector"
		<< "\n[10] Rotar matriz"
		<< "\n[11] Trasladar matriz"
		<< "\n[12] Trazo de elipses"
		<< "\n[13] Trazo de lineas"
		<< "\n[14] Trazo de curvas"
		<< "\n[15] Salir"
		<< "\n-------------------------------"
		<< "\nIngrese el numero de la opcion deseada: ";
}

//---------------------------------------------------------------------------------

//OPCIONES DE MENU-------------------------------------------------------------------

//1 * BIEN *
void sumarMatrices()
{
	int fil, col;

	cout << "\nLAS MATRICES DEBEN TENER LA MISMA DIMENSION\n"
		<< "\nDimensiones de la matriz:";
	cout << "\nFilas de la matriz: "; cin >> fil;
	cout << "Columnas de la matriz: "; cin >> col;

	int** A = new int* [fil];
	for (int i = 0; i < fil; i++)
		A[i] = new int[col];

	int** B = new int* [fil];
	for (int i = 0; i < fil; i++)
		B[i] = new int[col];

	int** C = new int* [fil];
	for (int i = 0; i < fil; i++)
		C[i] = new int[col];

	cout << "\nLlenar la matriz A: ";
	llenarMatriz(A, fil, col);

	cout << "\nLlenar la matriz B: ";
	llenarMatriz(B, fil, col);

	for (int i = 0; i < fil; i++)
		for (int j = 0; j < col; j++)
			C[i][j] = A[i][j] + B[i][j];

	cout << "\nMatriz A:\n";
	imprimirMatriz(A, fil, col);
	cout << "\nMatriz B:\n";
	imprimirMatriz(B, fil, col);
	cout << "\nSuma de las matrices (A + B):\n";
	imprimirMatriz(C, fil, col);
	limpiarPantalla();
}
//2 * BIEN *
void restarMatrices()
{
	int fil, col;

	cout << "\nLAS MATRICES DEBEN TENER LA MISMA DIMENSION\n"
		<< "\nDimensiones de la matriz:";
	cout << "\nFilas de la matriz: "; cin >> fil;
	cout << "Columnas de la matriz: "; cin >> col;

	int** A = new int* [fil];
	for (int i = 0; i < fil; i++)
		A[i] = new int[col];

	int** B = new int* [fil];
	for (int i = 0; i < fil; i++)
		B[i] = new int[col];

	int** C = new int* [fil];
	for (int i = 0; i < fil; i++)
		C[i] = new int[col];

	cout << "\nLlenar la matriz A: ";
	llenarMatriz(A, fil, col);

	cout << "\nLlenar la matriz B: ";
	llenarMatriz(B, fil, col);

	for (int i = 0; i < fil; i++)
		for (int j = 0; j < col; j++)
			C[i][j] = A[i][j] - B[i][j];

	cout << "\nMatriz A:\n";
	imprimirMatriz(A, fil, col);
	cout << "\nMatriz B:\n";
	imprimirMatriz(B, fil, col);
	cout << "\nResta de las matrices (A - B):\n";
	imprimirMatriz(C, fil, col);
	limpiarPantalla();
}
//3 * BIEN *
void multiplicarMatrices()
{
	int filA, colA, filB, colB;
	int rotar, angulo;

	cout << "\nEL NUMERO DE COLUMNAS DE LA MATRIZ A"
		<< "\nDEBE COINCIDIR CON EL DE FILAS DE LA MATRIZ B\n";

	cout << "\nDimensiones de la matriz A:";
	cout << "\nFilas de la matriz: "; cin >> filA;
	cout << "Columnas de la matriz: "; cin >> colA;

	int** A = new int* [filA];
	for (int i = 0; i < filA; i++)
		A[i] = new int[colA];

	llenarMatriz(A, filA, colA);

	cout << "\nDimensiones de la matriz B:";
	cout << "\nFilas de la matriz: "; cin >> filB;
	cout << "Columnas de la matriz: "; cin >> colB;

	int** B = new int* [filB];
	for (int i = 0; i < filB; i++)
		B[i] = new int[colB];

	llenarMatriz(B, filB, colB);

	int** C = new int* [filA];
	for (int i = 0; i < filA; i++)
		C[i] = new int[colB];

	if (colA == filB) {
		for (int i = 0; i < filA; ++i) {
			for (int j = 0; j < colB; ++j) {
				C[i][j] = 0;
				for (int z = 0; z < colA; ++z)
					C[i][j] += A[i][z] * B[z][j];
			}
		}

		cout << "\nMatriz A:\n";
		imprimirMatriz(A, filA, colA);
		cout << "\nMatriz B:\n";
		imprimirMatriz(B, filB, colB);
		cout << "\nMultiplicacion de las matrices (A * B):\n";
		imprimirMatriz(C, filA, colB);
		limpiarPantalla();
	}
	else
		cout << "\nNO SE PUEDEN MULTIPLICAR"
		<< "\nEL NUMERO DE COLUMNAS DE LA MATRIZ A"
		<< "\nDEBE COINCIDIR CON EL DE FILAS DE LA MATRIZ B";
}
//4 * BIEN *
void obtenerMatrizInversa()
{
	int fil, col;

	cout << "\nIngrese el numero de filas de la matriz: ";
	cin >> fil;
	cout << "\nIngrese el numero de columnas de la matriz: ";
	cin >> col;
	if (fil != col) {
		cout << "\nLa matriz debe ser cuadrada para tener una inversa.\n";
		return;
	}

	// Crear una matriz extendida [A | I]
	float** A = new float* [fil];
	for (int i = 0; i < fil; i++) {
		A[i] = new float[2 * col];
	}

	// Llenar la matriz A y la matriz identidad I
	cout << "\nIngrese los elementos de la matriz:\n";
	for (int i = 0; i < fil; i++) {
		for (int j = 0; j < col; j++) {
			cout << "Elemento [" << i << "][" << j << "]: ";
			cin >> A[i][j];
			if (i == j) {
				A[i][j + col] = 1.0; // Inicializar la matriz identidad
			}
			else {
				A[i][j + col] = 0.0;
			}
		}
	}

	// Aplicar el algoritmo de Gauss-Jordan
	for (int i = 0; i < fil; i++) {
		float pivot = A[i][i];
		for (int j = 0; j < 2 * col; j++) {
			A[i][j] /= pivot;
		}
		for (int k = 0; k < fil; k++) {
			if (k != i) {
				float factor = A[k][i];
				for (int j = 0; j < 2 * col; j++) {
					A[k][j] -= factor * A[i][j];
				}
			}
		}
	}

	// Extraer la matriz inversa
	float** inversa = new float* [fil];
	for (int i = 0; i < fil; i++) {
		inversa[i] = new float[col];
		for (int j = 0; j < col; j++) {
			inversa[i][j] = A[i][j + col];
		}
	}

	// Imprimir la matriz inversa
	cout << "\nMatriz Inversa:\n";
	imprimirMatrizInversa(inversa, fil, col);

	// Liberar memoria
	for (int i = 0; i < fil; i++) {
		delete[] A[i];
		delete[] inversa[i];
	}
	delete[] A;
	delete[] inversa;

	limpiarPantalla();
}
//5 * BIEN *
void obtenerMatrizTranspuesta()
{
	int fil, col;

	cout << "\nIngrese el numero de filas de la matriz: ";
	cin >> fil;
	cout << "\nIngrese el numero de columnas de la matriz: ";
	cin >> col;

	int** M = new int* [fil];
	for (int i = 0; i < fil; i++) {
		M[i] = new int[col];
	}

	cout << "\nIngrese los elementos de la matriz:\n";
	llenarMatriz(M, fil, col);

	int** transpuesta = new int* [col];
	for (int i = 0; i < col; i++) {
		transpuesta[i] = new int[fil];
	}

	for (int i = 0; i < fil; i++) {
		for (int j = 0; j < col; j++) {
			transpuesta[j][i] = M[i][j];
		}
	}

	cout << "\nMatriz Original:\n";
	imprimirMatriz(M, fil, col);
	cout << "\nMatriz Transpuesta:\n";
	imprimirMatriz(transpuesta, col, fil);

	// Liberar memoria
	for (int i = 0; i < fil; i++) {
		delete[] M[i];
	}
	delete[] M;

	for (int i = 0; i < col; i++) {
		delete[] transpuesta[i];
	}
	delete[] transpuesta;
	limpiarPantalla();
}
//6 * BIEN *
int obtenerDeterminante(int** matriz, int fil, int col)
{
	if (fil == 2 && col == 2) {
		return matriz[0][0] * matriz[1][1] - matriz[1][0] * matriz[0][1];
	}

	int determinante = 0;

	for (int i = 0; i < col; i++) {
		int** submatriz = new int* [fil - 1];
		for (int j = 0; j < fil - 1; j++) {
			submatriz[j] = new int[col - 1];
		}

		for (int j = 1; j < fil; j++) {
			for (int k = 0; k < col; k++) {
				if (k < i) {
					submatriz[j - 1][k] = matriz[j][k];
				}
				else if (k > i) {
					submatriz[j - 1][k - 1] = matriz[j][k];
				}
			}
		}

		int signo = (i % 2 == 0) ? 1 : -1;

		determinante += signo * matriz[0][i] * obtenerDeterminante(submatriz, fil - 1, col - 1);

		for (int j = 0; j < fil - 1; j++) {
			delete[] submatriz[j];
		}
		delete[] submatriz;
	}

	return determinante;
}
void obtenerMatrizAdjunta()
{
	int fil, col;

	cout << "\nIngrese el numero de filas de la matriz: ";
	cin >> fil;
	cout << "\nIngrese el numero de columnas de la matriz: ";
	cin >> col;

	int** A = new int* [fil];
	for (int i = 0; i < fil; i++) {
		A[i] = new int[col];
	}

	cout << "\nIngrese los elementos de la matriz:\n";
	llenarMatriz(A, fil, col);

	int** adjunta = new int* [fil];
	for (int i = 0; i < fil; i++) {
		adjunta[i] = new int[col];
	}

	for (int i = 0; i < fil; i++) {
		for (int j = 0; j < col; j++) {
			// Calcular el cofactor de A[i][j]
			int** cofactor = new int* [fil - 1];
			for (int k = 0; k < fil - 1; k++) {
				cofactor[k] = new int[col - 1];
			}

			int signo = ((i + j) % 2 == 0) ? 1 : -1;
			int subfil = 0, subcol = 0;

			for (int k = 0; k < fil; k++) {
				for (int l = 0; l < col; l++) {
					if (k != i && l != j) {
						cofactor[subfil][subcol++] = A[k][l];

						if (subcol == col - 1) {
							subfil++;
							subcol = 0;
						}
					}
				}
			}

			// Calcular el determinante del cofactor
			int determinanteCofactor = 0;
			if (fil == 2 && col == 2) {
				determinanteCofactor = cofactor[0][0] * cofactor[1][1] - cofactor[1][0] * cofactor[0][1];
			}
			else {
				determinanteCofactor = obtenerDeterminante(cofactor, fil - 1, col - 1);
			}

			// Asignar el cofactor multiplicado por el signo al elemento de la adjunta
			adjunta[i][j] = signo * determinanteCofactor;

			// Liberar memoria del cofactor
			for (int k = 0; k < fil - 1; k++) {
				delete[] cofactor[k];
			}
			delete[] cofactor;
		}
	}

	// Transponer la matriz adjunta
	int** transpuesta = new int* [col];
	for (int i = 0; i < col; i++) {
		transpuesta[i] = new int[fil];
		for (int j = 0; j < fil; j++) {
			transpuesta[i][j] = adjunta[j][i];
		}
	}

	// Imprimir la matriz adjunta
	cout << "\nMatriz Adjunta:\n";
	imprimirMatriz(transpuesta, col, fil);

	// Liberar memoria
	for (int i = 0; i < fil; i++) {
		delete[] A[i];
		delete[] adjunta[i];
	}
	delete[] A;
	delete[] adjunta;

	for (int i = 0; i < col; i++) {
		delete[] transpuesta[i];
	}
	delete[] transpuesta;
}
//7 * BIEN *
void obtenerPunto() {
	float x, y, z, zp, nP;
	float xp = 0, yp = 0;

	cout << "\nIngrese el numero de puntos a evaluar: ";
	cin >> nP;
	cout << "\nIngrese perspectiva z : ";
	cin >> zp;

	for (int i = 0; i < nP; i++) {
		cout << "\nPunto " << i + 1 << ":" << endl;
		cout << "\nIngrese valor de x: ";
		cin >> x;
		cout << "\nIngrese valor de y: ";
		cin >> y;
		cout << "\nIngrese valor de z: ";
		cin >> z;

		xp = x * (zp / (zp - z));
		yp = y * (zp / (zp - z));

		punto[i].pPrimoX = xp;
		punto[i].pPrimoY = yp;
		cout << "\n";

	}

	for (int i = 0; i < nP; i++) {
		cout << "Punto " << i + 1 << ":" << endl;
		cout << fixed << setprecision(4);
		cout << "\nxp = " << punto[i].pPrimoX;
		cout << "\nyp = " << punto[i].pPrimoY;
		cout << "\n\n";
	}
	limpiarPantalla();
}
//8 * BIEN *
void rotatePoint(const Quaternion& point, const Quaternion& vector, double angle, Quaternion& result) {
	// Convertir el ángulo a radianes
	angle = angle * acos(-1.0) / 180.0;

	// Normalizar el vector de rotación
	Quaternion axis = normalizeQuaternion(vector);

	// Calcular el cuaternión de rotación
	Quaternion rotation;
	rotation.w = cos(angle / 2.0);
	rotation.x = axis.x * sin(angle / 2.0);
	rotation.y = axis.y * sin(angle / 2.0);
	rotation.z = axis.z * sin(angle / 2.0);

	// Calcular el cuaternión resultante
	result = multiplyQuaternions(multiplyQuaternions(rotation, point), conjugate(rotation));
}
Quaternion normalizeQuaternion(const Quaternion& q) {
	double norm = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
	Quaternion result;
	result.w = q.w / norm;
	result.x = q.x / norm;
	result.y = q.y / norm;
	result.z = q.z / norm;
	return result;
}
Quaternion multiplyQuaternions(const Quaternion& q1, const Quaternion& q2) {
	Quaternion result;
	result.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
	result.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
	result.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
	result.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
	return result;
}
void GirarPuntoVector()
{
	// Pedir al usuario las coordenadas del punto P
	cout << "Ingrese las coordenadas del punto P (x y z): ";
	double px, py, pz;
	cin >> px >> py >> pz;

	// Pedir al usuario las coordenadas del vector V
	cout << "Ingrese las coordenadas del vector V (x y z): ";
	double vx, vy, vz;
	cin >> vx >> vy >> vz;

	// Pedir al usuario el ángulo de rotación
	cout << "Ingrese el angulo de rotacion en grados: ";
	double angle;
	cin >> angle;

	// Crear cuaterniones para el punto y el vector
	Quaternion point = { 0, px, py, pz };
	Quaternion vector = { 0, vx, vy, vz };

	// Rotar el punto alrededor del vector
	Quaternion result;
	rotatePoint(point, vector, angle, result);

	// Mostrar el resultado
	cout << "El punto rotado es: (" << result.x << ", " << result.y << ", " << result.z << ")" << endl;
	limpiarPantalla();
}
//9 *MAL - HAY QUE HALLAR MATRIZ COMPUESTA *
Punto rotarPunto(const Punto& punto, const Punto& eje1, const Punto& eje2, double grados) {
	// Conversión de grados a radianes
	double radianes = grados * acos(-1.0) / 180.0;

	// Cálculo del vector de rotación
	Punto vRotacion;
	vRotacion.x = eje2.x - eje1.x;
	vRotacion.y = eje2.y - eje1.y;
	vRotacion.z = eje2.z - eje1.z;

	// Normalización del vector de rotación
	double longitud = sqrt(vRotacion.x * vRotacion.x + vRotacion.y * vRotacion.y + vRotacion.z * vRotacion.z);
	vRotacion.x /= longitud;
	vRotacion.y /= longitud;
	vRotacion.z /= longitud;

	// Cálculo de la matriz de rotación
	double cosTheta = cos(radianes);
	double sinTheta = sin(radianes);
	double unoMenosCosTheta = 1 - cosTheta;

	double matrizRotacion[3][3] = {
		{cosTheta + vRotacion.x * vRotacion.x * unoMenosCosTheta, vRotacion.x * vRotacion.y * unoMenosCosTheta - vRotacion.z * sinTheta, vRotacion.x * vRotacion.z * unoMenosCosTheta + vRotacion.y * sinTheta},
		{vRotacion.y * vRotacion.x * unoMenosCosTheta + vRotacion.z * sinTheta, cosTheta + vRotacion.y * vRotacion.y * unoMenosCosTheta, vRotacion.y * vRotacion.z * unoMenosCosTheta - vRotacion.x * sinTheta},
		{vRotacion.z * vRotacion.x * unoMenosCosTheta - vRotacion.y * sinTheta, vRotacion.z * vRotacion.y * unoMenosCosTheta + vRotacion.x * sinTheta, cosTheta + vRotacion.z * vRotacion.z * unoMenosCosTheta}
	};

	// Aplicación de la matriz de rotación al punto
	Punto puntoRotado;
	puntoRotado.x = punto.x * matrizRotacion[0][0] + punto.y * matrizRotacion[0][1] + punto.z * matrizRotacion[0][2];
	puntoRotado.y = punto.x * matrizRotacion[1][0] + punto.y * matrizRotacion[1][1] + punto.z * matrizRotacion[1][2];
	puntoRotado.z = punto.x * matrizRotacion[2][0] + punto.y * matrizRotacion[2][1] + punto.z * matrizRotacion[2][2];

	return puntoRotado;
}
void GirarFiguraVector()
{
	int n;
	cout << "Ingrese la cantidad de puntos: ";
	cin >> n;

	// Vector que almacenará los puntos del prisma
	vector<Punto> puntos;

	cout << "Ingrese los puntos del prisma (x y z separados por espacio):" << endl;
	for (int i = 0; i < n; ++i) {
		Punto punto;
		cin >> punto.x >> punto.y >> punto.z;
		puntos.push_back(punto);
	}

	Punto vp1, vp2;
	cout << "Ingrese los puntos vp1 y vp2 para definir el eje de rotacion (x y z separados por espacio para cada punto):" << endl;
	cin >> vp1.x >> vp1.y >> vp1.z;
	cin >> vp2.x >> vp2.y >> vp2.z;

	double grados;
	cout << "Ingrese el angulo de rotacion en grados: ";
	cin >> grados;

	// Rotar cada punto del prisma alrededor del eje definido por vp1 y vp2
	for (int i = 0; i < n; ++i) {
		Punto puntoRotado = rotarPunto(puntos[i], vp1, vp2, grados);
		cout << "Punto original: (" << puntos[i].x << ", " << puntos[i].y << ", " << puntos[i].z << ") ";
		cout << "Punto rotado: (" << puntoRotado.x << ", " << puntoRotado.y << ", " << puntoRotado.z << ")" << endl;
	}
}
//10 *  *
//************************
void rotarMatriz(float** M, int fil, int col, char eje, int angulo)
{
	// Convertir el ángulo a radianes
	double radianes = angulo * M_PI / 180.0;

	// Matrices de rotación tridimensional
	float matrizRotacionX[3][3] = {
		{1, 0, 0},
		{0, cos(radianes), -sin(radianes)},
		{0, sin(radianes), cos(radianes)}
	};

	float matrizRotacionY[3][3] = {
		{cos(radianes), 0, sin(radianes)},
		{0, 1, 0},
		{-sin(radianes), 0, cos(radianes)}
	};

	float matrizRotacionZ[3][3] = {
		{cos(radianes), -sin(radianes), 0},
		{sin(radianes), cos(radianes), 0},
		{0, 0, 1}
	};

	// Realizar la rotación según el eje especificado
	switch (eje)
	{
	case 'x':
		for (int i = 0; i < fil; ++i) {
			float yOriginal = M[i][1];
			float zOriginal = M[i][2];

			M[i][1] = matrizRotacionX[1][1] * yOriginal + matrizRotacionX[1][2] * zOriginal;
			M[i][2] = matrizRotacionX[2][1] * yOriginal + matrizRotacionX[2][2] * zOriginal;
		}
		break;
	case 'y':
		for (int i = 0; i < fil; ++i) {
			float xOriginal = M[i][0];
			float zOriginal = M[i][2];

			M[i][0] = matrizRotacionY[0][0] * xOriginal + matrizRotacionY[0][2] * zOriginal;
			M[i][2] = matrizRotacionY[2][0] * xOriginal + matrizRotacionY[2][2] * zOriginal;
		}
		break;
	case 'z':
		for (int i = 0; i < fil; ++i) {
			float xOriginal = M[i][0];
			float yOriginal = M[i][1];

			M[i][0] = matrizRotacionZ[0][0] * xOriginal + matrizRotacionZ[0][1] * yOriginal;
			M[i][1] = matrizRotacionZ[1][0] * xOriginal + matrizRotacionZ[1][1] * yOriginal;
		}
		break;
	default:
		cout << "Eje de rotacion no valido. Use 'x', 'y' o 'z'.";
		return;
	}

	// Imprimir la matriz original y la matriz rotada
	cout << "\nMatriz Original:\n";
	imprimirMatrizInversa(M, fil, col);

	limpiarPantalla();
}
void otraFunc() {
	int fil, col, angulo;
	float** M;

	cout << "\nIngrese el numero de filas de la matriz: ";
	cin >> fil;
	cout << "Ingrese el numero de columnas de la matriz: ";
	cin >> col;

	M = new float* [fil];
	for (int i = 0; i < fil; i++)
		M[i] = new float[col];

	cout << "\nIngrese los elementos de la matriz:\n";
	llenarMatrizF(M, fil, col);

	cout << "\nIngrese el angulo de rotacion en grados: ";
	cin >> angulo;

	char eje;
	cout << "Ingrese el eje de rotacion ('x', 'y', o 'z'): ";
	cin >> eje;

	// Realizar la rotación
	rotarMatriz(M, fil, col, eje, angulo);

	// Liberar memoria
	for (int i = 0; i < fil; i++)
		delete[] M[i];
	delete[] M;
}
//************************
//11 * NI BIEN NI MAL - NO DEBERIA SER OPCION *
void obtenerPuntoTranformaciones(int dx, int dy, int dz) {
	int fil = 4, col = 4;
	int nP;
	float x, y, z, zp;
	float xp, yp;

	cout << "\nIngrese el numero de puntos a evaluar: ";
	cin >> nP;
	cout << "Ingrese z perspectiva: ";
	cin >> zp;

	for (int i = 0; i < nP; i++) {
		cout << "\nPunto " << i + 1 << ":" << endl;
		cout << "\nIngrese valor de x: ";
		cin >> x;
		cout << "\nIngrese valor de y: ";
		cin >> y;
		cout << "\nIngrese valor de z: ";
		cin >> z;

		punto[i].pX = x + dx;
		punto[i].pY = y + dy;
		punto[i].pZ = z + dz;

		xp = x * (zp / (zp - z));
		yp = y * (zp / (zp - z));

		punto[i].pPrimoX = xp;
		punto[i].pPrimoY = yp;
	}

	cout << "\nPuntos Obtenidos:\n";

	for (int i = 0; i < nP; i++) {
		cout << "Punto " << i + 1 << ": (" << punto[i].pX << ", " << punto[i].pY << ", " << punto[i].pZ << " )\n\n";
	}

	cout << "\nPuntos transformados: \n";

	for (int i = 0; i < nP; i++) {
		cout << "Punto " << i + 1 << ":" << endl;
		cout << fixed << setprecision(4);
		cout << "\nxp = " << punto[i].pPrimoX;
		cout << "\nyp = " << punto[i].pPrimoY;
		cout << "\n\n";
	}
}
void trasladarMatriz() {
	float dx, dy, dz;
	int fil = 4, col = 4, nP = 0;

	float** M = new float* [fil];
	for (int i = 0; i < fil; i++)
		M[i] = new float[col];

	cout << "Desplazamiento en X: ";
	cin >> dx;
	cout << "Desplazamiento en Y: ";
	cin >> dy;
	cout << "Desplazamiento en Z: ";
	cin >> dz;

	// Crear una matriz de traslación
	float** matrizTraslacion = new float* [4];
	for (int i = 0; i < 4; i++) {
		matrizTraslacion[i] = new float[4];
		for (int j = 0; j < 4; j++) {
			matrizTraslacion[i][j] = 0;
		}
	}

	matrizTraslacion[0][0] = 1;
	matrizTraslacion[1][1] = 1;
	matrizTraslacion[2][2] = 1;
	matrizTraslacion[3][3] = 1;
	matrizTraslacion[0][3] = dx;
	matrizTraslacion[1][3] = dy;
	matrizTraslacion[2][3] = dz;

	obtenerPuntoTranformaciones(dx, dy, dz);

	// Liberar la memoria de las matrices auxiliares
	for (int i = 0; i < 4; i++) {
		delete[] matrizTraslacion[i];
	}
	delete[] matrizTraslacion;

	limpiarPantalla();
}
//12 * MAL - CORREGIR *
//**********************
void elipseBresenham() {
	int xc, yc, rx, ry;
	cout << "Ingrese las coordenadas del centro separados por un espacio (xc & yc): ";
	cin >> xc >> yc;
	cout << "Ingrese los radios en x & en y: ";
	cin >> rx >> ry;
	reglaElipse(xc, yc, rx, ry);
}
void reglaElipse(int xc, int yc, int rx, int ry) {
	int x = 0, y = ry;
	int rxSq = rx * rx;
	int rySq = ry * ry;
	int fx = 2 * rySq - 2 * rxSq * ry + rxSq;
	int fy = 2 * rxSq - 2 * rySq * rx + rySq;

	// Crear una matriz para representar la gráfica
	vector<vector<char>> grafica(yc + ry + 3, vector<char>(xc + rx + 3, ' '));

	cout << "\n\nPuntos Generados:\n\n";
	grafica[yc + 1][xc + 1] = 'o';  // Marcar el punto inicial

	while (x <= y) {
		if (fx < 0) {
			fx += 2 * rySq * (2 * x + 3);
			x++;
		}
		else if (fy > 0) {
			fy -= 2 * rxSq * (2 * y - 3);
			y--;
		}
		else {
			fx += 2 * rySq * (2 * x + 3);
			fy -= 2 * rxSq * (2 * y - 3);
			x++;
			y--;
		}
		cout << "(" << xc + x << ", " << yc + y << ")" << endl;
		grafica[yc + y + 1][xc + x + 1] = 'o';  // Marcar cada punto generado
	}

	// Imprimir la gráfica
	for (int i = grafica.size() - 1; i >= 0; i--) {
		for (int j = 0; j < grafica[i].size(); j++) {
			cout << grafica[i][j];
		}
		cout << "\n";
	}
	limpiarPantalla();
}
//**********************
//13 * BIEN - CORREGIDO *
//**********************
void trazoBresenham() {
	int x1, y1, x2, y2;
	cout << "Ingrese las coordenadas del primer punto separados por un espacio (x1 y  y1): ";
	cin >> x1 >> y1;

	cout << "Ingrese las coordenadas del segundo punto separados por un espacio (x2 y y2): ";
	cin >> x2 >> y2;

	reglaTrazo(x1, x2, y1, y2);
}
void reglaTrazo(int x1, int x2, int y1, int y2) {
	int dx, dy, pO, pK, repetir;

	dx = x2 - x1;
	dy = y2 - y1;
	repetir = dx - 1;

	pO = 2 * (dy)-dx;

	// Crear una matriz para representar la gráfica
	vector<vector<char>> grafica(y2 + 3, vector<char>(x2 + 3, ' '));

	cout << "\n\nPuntos Generados:\n\n(" << x1 << ", " << y1 << ")" << endl;
	grafica[y1 + 1][x1 + 1] = 'o';  // Marcar el punto inicial

	for (int i = 0; i < repetir; i++) {
		if (pO < 0) {
			pK = pO + 2 * dy;
		}
		else {
			pK = pO + 2 * dy - 2 * dx;
			y1++;
		}
		x1++;
		pO = pK;
		cout << "(" << x1 << ", " << y1 << ")" << endl;
		grafica[y1 + 1][x1 + 1] = 'o';  // Marcar cada punto generado
	}

	// Imprimir la gráfica
	for (int i = grafica.size() - 1; i >= 0; i--) {
		for (int j = 0; j < grafica[i].size(); j++) {
			cout << grafica[i][j];
		}
		cout << "\n";
	}

	limpiarPantalla();
}
//************************
//14  * BIEN - CORREGIDO *
//**********************
void curvasBresenham() {
	int xc, yc, r;
	cout << "Ingrese las coordenadas del centro separados por un espacio (xc y yc): ";
	cin >> xc >> yc;
	cout << "Ingrese el radio: ";
	cin >> r;
	reglaCurvas(xc, yc, r);
}
void reglaCurvas(int xc, int yc, int r) {
	int x = 0, y = r;
	int d = 3 - 2 * r;

	// Crear una matriz para representar la gráfica
	vector<vector<char>> grafica(yc + r + 3, vector<char>(xc + r + 3, ' '));

	cout << "\n\nPuntos Generados:\n\n(" << xc << ", " << yc << ")" << endl;
	grafica[yc + 1][xc + 1] = 'o';  // Marcar el punto inicial

	while (y >= x) {
		x++;
		if (d > 0) {
			y--;
			d = d + 4 * (x - y) + 10;
		}
		else {
			d = d + 4 * x + 6;
		}
		cout << "(" << xc + x << ", " << yc + y << ")" << endl;
		grafica[yc + y + 1][xc + x + 1] = 'o';  // Marcar cada punto generado
	}

	// Imprimir la gráfica
	for (int i = grafica.size() - 1; i >= 0; i--) {
		for (int j = 0; j < grafica[i].size(); j++) {
			cout << grafica[i][j];
		}
		cout << "\n";
	}
	limpiarPantalla();
}
//************************

//---------------------------------------------------------------------------------