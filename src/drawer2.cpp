#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

using namespace std;

#define SCENE_BOUNDING_SPHERE_RADIUS (3.1)
#define CAMERA_GAP  (0.9)
#define CAMERA_STEP (0.5)
#define UP (1)
#define DOWN (0)

double a = 1.0;
double l = 3.0;
double h0 = 0.02;
double dx = 0.001;
const double kB = 4.5;// heigth of matrix

ifstream free_data("data/free_new.dat");
ifstream constrained_data_x("data/constrained_x.dat");
ifstream constrained_data_y("data/constrained_y.dat");

int stadia = 1;

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

void draw_matrix(){
  //right side
  glPushMatrix();
    glTranslatef(1.05, kB/2, 0);
    glScalef(0.1, kB, 3);
    glutSolidCube(1);
  glPopMatrix();
  //left side
  glPushMatrix();
    glTranslatef(-1.05, kB/2, 0);
    glScalef(0.1, kB, 3);
    glutSolidCube(1);
  glPopMatrix();
  //bottom
  glPushMatrix();
    glTranslatef(0, kB, 0);
    glScalef(2.2, 0.1, 3);
    glutSolidCube(1);
  glPopMatrix();
}

void set_camera() {

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(50.0, 1.0, 0.5, 200);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-1.0, -1.0, -global.camera.distance);

  glRotatef(global.camera.rotation_x, 1.0, 0.0, 0.0);
  glRotatef(global.camera.rotation_y, 0.0, 1.0, 0.0);
  glRotatef(global.camera.rotation_z, 0.0, 0.0, 1.0);
}

void half_cylinder(double center, double rho, double h, GLUquadricObj *quadObj1, double cliping){
  double clip_plane2[] = {0, 1.0, 0 , -cliping};
  glEnable(GL_CLIP_PLANE1);

  glPushMatrix();
    glTranslatef(0, center, -l/2.0);
    gluCylinder(quadObj1, rho, rho, l, 100, 100);
  glPopMatrix();

  glClipPlane(GL_CLIP_PLANE1,clip_plane2);

  glPushMatrix();
    glTranslatef(0, center, -l/2.0);
    gluCylinder(quadObj1, rho-h, rho-h, l, 100, 100);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(0, center, l/2.0);
    gluDisk(quadObj1, rho - h, rho, 100, 100);
  glPopMatrix();


  glPushMatrix();
    glTranslatef(0, center, -l/2.0);
    gluDisk(quadObj1, rho - h, rho, 100, 100);
  glPopMatrix();

  glDisable(GL_CLIP_PLANE1);
}

void quarter_cylinder(double center_x, double center_y, double rho, double h, GLUquadricObj *quadObj1, double *clip_plane1, double cliping_y){
  double clip_plane2[] = {0.0, 1.0, 0.0, -cliping_y};
  // double clip_plane1[] = {1.0, 0.0, 0.0, -cliping_x};

  glEnable(GL_CLIP_PLANE1);
  glEnable(GL_CLIP_PLANE2);


  glPushMatrix();
    glTranslatef(center_x, center_y, -l/2.0);
    gluCylinder(quadObj1, rho, rho, l, 100, 100);
  glPopMatrix();

  glClipPlane(GL_CLIP_PLANE2,clip_plane2);
  glClipPlane(GL_CLIP_PLANE1,clip_plane1);

  glPushMatrix();
    glTranslatef(center_x, center_y, -l/2.0);
    gluCylinder(quadObj1, rho-h, rho-h, l, 100, 100);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(center_x, center_y, l/2.0);
    gluDisk(quadObj1, rho - h, rho, 100, 100);
  glPopMatrix();


  glPushMatrix();
    glTranslatef(center_x, center_y, -l/2.0);
    gluDisk(quadObj1, rho - h, rho, 100, 100);
  glPopMatrix();

  glDisable(GL_CLIP_PLANE1);
  glDisable(GL_CLIP_PLANE2);
}

void borders(double y, double h){
  glPushMatrix();
    glTranslatef(a-h, y/2.0, 0);
    glScalef(h, y, l);
    glutSolidCube(1);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(-a+h, y/2.0, 0);
    glScalef(h, y, l);
    glutSolidCube(1);
  glPopMatrix();
}

void display(void) {
  double t, var, h, rho, center;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glColor3ub(255,255,0);

  set_camera();
  renderAxes();

  glColor4f(0.3, 0.3, 0.3,0.5);

  GLUquadricObj *quadObj1;

  quadObj1 = gluNewQuadric();

  // ********** fixed parts ***********//
  glColor4f(0.5,0.5,0.5, 1.0);
  glPushMatrix();
    glTranslatef(-a*1.6, 0, 0);
    glScalef(a, h0, l);
    glutSolidCube(1.0); 
    glTranslatef(3.1*a, 0, 0);
    glutSolidCube(1.0); 
  glPopMatrix();


 
  glColor4f(0.7,0.3,0.9, 1.0);

  if (free_data.fail())
    stadia = 2;

  if (constrained_data_y.fail())
    stadia = 3;

  if (constrained_data_x.fail()){
    sleep(10000);
    glutSwapBuffers();
    return; 
  }
    

  if(stadia == 1)
    free_data >> t >> var >> h >> rho >> center;
  else if(stadia == 2){
    constrained_data_y >> t >> var >> h >> rho >> center;
  } else {
    constrained_data_x >> t >> var >> h >> rho >> center;
  }

  if(stadia == 1){ 
    half_cylinder(center, rho, h, quadObj1, h0/2.0);

  } else if (stadia==2) {
    borders(var, h);
    half_cylinder(center, rho, h, quadObj1, var);
    // usleep(10000);
  } else {
    double cpx[] = {1.0, 0.0, 0.0,  -var};
    double cpx2[] = {-1.0, 0.0, 0.0, -var};
    borders(kB-(a-var), h);

    // bottom part
    glPushMatrix();
      glTranslatef(0, kB-h, 0);
      glScalef(var, h, l);
      glutSolidCube(1);
    glPopMatrix();

    glPopMatrix();

    quarter_cylinder(var,  kB-(a-var), a-var, h, quadObj1, cpx, kB-(a-var));
    quarter_cylinder(-var, kB-(a-var), a-var, h, quadObj1, cpx2, kB-(a-var));
    usleep(100000);
  }

  glColor4f(0.3, 0.3, 0.3,0.5);
  draw_matrix();


  gluDeleteQuadric(quadObj1);
  glutSwapBuffers();

}

void init_global(){
  global.camera.distance = 20.0;
  global.camera.rotation_x = 0.0;
  global.camera.rotation_y = 0.0;
  global.camera.rotation_z = 0.0;
  global.mouse.sensetivity = 0.25;
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
  
  init_global();

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_NORMALIZE);
  glEnable(GL_AUTO_NORMAL);
  
  quadObj = gluNewQuadric();
  glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
  glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, white_light);
}

void reshape(GLint width, GLint height) {
  global.screen.width  = width;
  global.screen.height = height;

  glViewport(0, 0, width, height);
}

void motion(int x, int y) {
  global.camera.rotation_y += global.mouse.sensetivity * (x - global.mouse.x); 
  global.camera.rotation_x += global.mouse.sensetivity * (y - global.mouse.y);

  global.mouse.x = x;
  global.mouse.y = y;

  glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
  const static unsigned char KEY_ESCAPE = '\033';

  if (key == '-') {
    global.camera.distance += CAMERA_STEP;
  }

  if (key == '+') {
    global.camera.distance -= CAMERA_STEP;
    if (global.camera.distance < SCENE_BOUNDING_SPHERE_RADIUS + CAMERA_GAP)
      global.camera.distance = SCENE_BOUNDING_SPHERE_RADIUS + CAMERA_GAP;
  }

  if (key == KEY_ESCAPE) {
    exit(EXIT_SUCCESS);
  }

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  global.mouse.x = x;
  global.mouse.y = y;

  glutPostRedisplay();
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
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(spinDisplay);
  glutMainLoop();
  return 0;  
}