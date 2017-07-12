#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

using namespace std;

#include <GL/freeglut.h>

int N = 10;                      // # de moleculas
double L = 10;                    // Lado del cuadrado
double kT = 1;                   // Energia cinetica inicial
int skip = 0;                    // Salto de tiempo para la simulación

void ingredar_datos() {
    cout << "Dinamica Molecular de Lennard-Jones\n";
    cout << "Numero de moleculas N = ";
    cin >> N;
    cout << "Lado del cuadrado L = ";
    cin >> L;
    cout << "Energia iniciar cinetica por molecular: ";
    cin >> kT;
    cout << "Saltos de tiempo :";
    cin >> skip;
}

double *x, *y;    //Posiciones de las moleculas
double *vx, *vy;  //Velocidades de las moleculas
double *ax, *ay;  //Aceleracion de las moleculas

void calcularAceleraciones()
{
    for (int i = 1; i <= N; i++)
        ax[i] = ay[i] = 0;

    for (int i = 1; i < N; i++)
    for (int j = i + 1; j <= N; j++) {
        double dx = x[i] - x[j];
        double dy = y[i] - y[j];
        //Uso de la imagen periodica más cercana

        if (abs(dx) > 0.5 * L)
            dx *= 1 - L / abs(dx);
        if (abs(dy) > 0.5 * L)
            dy *= 1 - L / abs(dy);
        double dr = sqrt(dx * dx + dy * dy);
        double f = 48 * pow(dr, -13.0) - 24 * pow(dr, -7.0);
        ax[i] += f * dx / dr;
        ay[i] += f * dy / dr;
        ax[j] -= f * dx / dr;
        ay[j] -= f * dy / dr;

    }
}

double t = 0;                   // time
double dt = 0.01;               // integration time step
int step = 0;                   // step number
double T;                       // temperature
double Tsum;                    // to compute average T
int step0;                      // starting step for computing average

int nBins = 50;                 // Numero de velocidad de todas las moleculas
double *vBins;                  // para pa distribución de Maxwell-Boltzmann
double vMax = 4;                // Velocidad maxima total
double dv = vMax / nBins;       // tamaño de bin


void inicializacion() {

    // Inicializamos espacio para las posiciones
    x = new double [N+1];
    y = new double [N+1];
    // Inicializamos espacio para  las velocidades
    vx = new double [N+1];
    vy = new double [N+1];
    // Inicializamos espacio para  las aceleraciones
    ax = new double [N+1];
    ay = new double [N+1];
    vBins = new double[nBins+1];
    //inicializamos las posiciones de las moleculas

    int m = int(ceil(sqrt(double(N))));  // molecules in each direction
    double d = L / m;                    // espacio entre moleculas
    int n = 0;                           // molecules placed so far
    for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++) {
        if (n < N) {
            ++n;
            x[n] = (i + 0.5) * d;
            y[n] = (j + 0.5) * d;
        }
    }

    // inicializando velocidades con direcciones aleatorias
    double pi = 4 * atan(1.0);
    double v = sqrt(2 * kT); 
    for (int n = 1; n <= N; n++) {
        double theta = 2 * pi * rand() / double(RAND_MAX);
        vx[n] = v * cos(theta);
        vy[n] = v * sin(theta);
    }

    calcularAceleraciones();

    T = kT;
}

void timeStep () {

    t += dt;
    ++step;
    double K = 0;
    for (int i = 1; i <= N; i++) {

        // Integrando velocidades con el algoritmo de Verlet
        x[i] += vx[i] * dt + 0.5 * ax[i] * dt * dt;
        y[i] += vy[i] * dt + 0.5 * ay[i] * dt * dt;

        // periodic boundary conditions
        if (x[i] < 0) x[i] += L;
        if (x[i] > L) x[i] -= L;
        if (y[i] < 0) y[i] += L;
        if (y[i] > L) y[i] -= L;

        vx[i] += 0.5 * ax[i] * dt;
        vy[i] += 0.5 * ay[i] * dt;

        calcularAceleraciones();

        vx[i] += 0.5 * ax[i] * dt;
        vy[i] += 0.5 * ay[i] * dt;

        double v = sqrt(vx[i] * vx[i] + vy[i] * vy[i]);
        K += 0.5 * v * v;
        int bin = (int) (nBins * v / vMax);
        if (bin > 0 && bin <= nBins)
            ++vBins[bin];
    }
    Tsum += K / N;
    T = Tsum / (step - step0);
}



int mainWindow, ventana_moleculas;
int margin = 10;

void takeStep() {
    for (int i = 0; i < skip; i++)
        timeStep();
    timeStep();
    glutSetWindow(ventana_moleculas);
    glutPostRedisplay();
    if (step % 10 == 0) {

        glutPostRedisplay();
    }
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
}

void drawText(const string &str, double x, double y) {
    glRasterPos2d(x, y);
    int len = str.find('\0');
    for (int i = 0; i < len; i++)
       glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]);
}

void displayMols() {
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3ub(0, 0, 255);
    double pi = 4 * atan(1.0);
    // Dibuja la molecula a partir de triangulos
    for (int n = 1; n <= N; n++) {
        glBegin(GL_TRIANGLE_FAN);
            glVertex2d(x[n], y[n]);
            double phi = 2 * pi / 24;
            for (int j = 0; j < 25; j++)
                glVertex2d(x[n] + 0.5 * cos(phi*j), y[n] + 0.5 * sin(phi*j));
        glEnd();
    }

    glColor3ub(255, 255, 255); //Color
    ostringstream os;
    os << "Step No: " << step << "   Time t = " << t << ends;
    drawText(os.str(), L / 50, L / 50);
    os.seekp(0);
    glutSwapBuffers();
}

void reshape(int w, int h) {

    glutSetWindow(ventana_moleculas);
    int width = w - 2*margin;
    int height= h - 2*margin;
    glutReshapeWindow(width, height);
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, L, 0, L);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

bool running = false;

void molsMouse(int button, int state, int x, int y) {
    switch (button) {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            if (running)
            {
                glutIdleFunc(NULL); // Pausa cuando se hace click
                running = false;
            }
            else
            {
                glutIdleFunc(takeStep); // Se reanuda el programa
                cout<<"play"<<endl;
                running = true;
            }
        }
    }
}

void generarVentana() {
    // main window
    glutInitWindowSize(400, 400);
    glutInitWindowPosition(100, 100);
    mainWindow = glutCreateWindow("Dinamica Molecular - Lennard Jones");

    glShadeModel(GL_FLAT);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    // Ventana para las moleculas
    ventana_moleculas = glutCreateSubWindow(mainWindow, margin, margin,400 - 2*margin, 400 - 2*margin);
    glClearColor(1.0, 0.0, 0.0, 0.0);
    glutDisplayFunc(displayMols);
    glutMouseFunc(molsMouse);
}





int main(int argc, char *argv[])
{
    glutInit(&argc, argv);

    if (argc == 1)
        ingredar_datos(); // Ingresar el numero de moleculas, lado del cuadrado, energia cinetica y salto de tiempo
    inicializacion();  // Inicializamos velocidades, aceleraciones, posiciones y direcciones
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    generarVentana();
    glutMainLoop();
}
