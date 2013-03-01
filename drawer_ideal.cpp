#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "ideal_sliding.h"
#include "matrix_surface.h"
#include "free_deformation.h"
#include "membrane.h"

using namespace std;

#define SCENE_BOUNDING_SPHERE_RADIUS (3.1)
#define CAMERA_GAP  (0.9)
#define CAMERA_STEP (0.5)
#define UP (1)
#define DOWN (0)

double h0 = 0.02;
double q = 2650.0/88.3/1000000;
double n = 3.4;
double a = 1.0;
double SigmaB = 1.0;//88.3*1000000;

Membrane m(h0, q, n, SigmaB, a, 0.1, 999, 99);
FreeDeformation f(h0, q, n);
IdealSliding is(m);


double po(double alpha){
  return a/sin(alpha);
}

double center(double alpha){
  return -1*a/tan(alpha);
}

double h(double alpha){
  return sin(alpha)/alpha;
}

GLfloat spin=-9.0;
GLUquadricObj *quadObj;

typedef struct rgba {
  GLfloat r, g, b, a;
} rgba_t;

typedef struct {
  struct {
    GLint width, height;
  } screen;

  struct {
    GLdouble rotation_x, 
             rotation_y, 
             rotation_z;
    GLdouble distance;
  } camera;

  struct {
    GLint x, y;
    GLint active;
    GLdouble sensetivity;
  } mouse;

  rgba_t background;
  
  int animate;
  int direction;
} global_t;

global_t global;

void renderAxes(){
  glColor3ub(255,0,0);
  glBegin(GL_LINES); // red x
    glVertex3f(500,0,0);
    glVertex3f(0,0,0);
    
    glColor3ub(0,0,255); //blue y
    glVertex3f(0,500,0);
    glVertex3f(0,0,0);
  
    glColor3ub(0,250,0); //green z
    glVertex3f(0,0,500);
    glVertex3f(0,0,0);
  glEnd();
}

void set_camera() {

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(50.0, 1.0, 0.5, 200);

 
}

void draw_matrix(){
  glColor4f(0.2, 0.3, 0.32, 0.3);
  double x= -1;
  for(x = -1.0; x < 1; x+=0.01){
    glBegin(GL_POLYGON);
      glVertex3d(x, m.m_surface_(x), 0);
      glVertex3d(x, m.m_surface_(x), a);
      glVertex3d(x+0.01, m.m_surface_(x+0.01), a);
      glVertex3d(x+0.01, m.m_surface_(x+0.01), 0);
    glEnd(); 
  }
}



void draw_free(alpha){
  for(x = -1.0; x < 1; x+=0.01){
    glBegin(GL_POLYGON);
      glVertex3d(x, center(alpha), 0);
      glVertex3d(x, center(alpha), a);
      glVertex3d(x+0.01, center(x+0.01), a);
      glVertex3d(x+0.01, m.m_surface_(x+0.01), 0);
    glEnd(); 
  }
}

void display(void) {
  double t, x_0, alpha;
  
  renderAxes();
  
  draw_matrix();
  while(cin >> t >> alpha && t > 0){
    draw_free(alpha);
  }
  glutSwapBuffers();

}

void init(void) {
  GLfloat mat_specular[]={0.8,0.8,0.8,0.0};
	GLfloat mat_shininess[]={200.0};
	//GLfloat light_position[]={1.0,1.0,1.0,0.0};
	GLfloat white_light[]={0.8,0.8,0.8,1.0};
  //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glClearColor(0.0,0.0,0.0,0.0);
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_LIGHTING);
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_FLAT);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);

  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_NORMALIZE);
  glEnable(GL_AUTO_NORMAL);

  glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
  glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, white_light);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective( /* field of view in degree */ 50.0,
        /* aspect ratio */ 1.0,
        /* Z near */ 1.0,
        /* Z far */ 100.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(8, 1, 10, 0, 0, 0, 0, 1,0);
}

void reshape(GLint width, GLint height) {
  global.screen.width  = width;
  global.screen.height = height;

  glViewport(0, 0, width, height);
}

void spinDisplay(void){
  glutPostRedisplay();
}

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH|GLUT_STENCIL);
  glutInitWindowSize(512,512);
  glutInitWindowPosition(100,100);
  glutCreateWindow("Membrane");
  init();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(spinDisplay);
  glutMainLoop();
  return 0;  
}