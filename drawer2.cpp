#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

#define SCENE_BOUNDING_SPHERE_RADIUS (3.1)
#define CAMERA_GAP  (0.9)
#define CAMERA_STEP (0.5)
#define UP (1)
#define DOWN (0)

double a = 1.0;
double l = 3.0;
double h0 = 0.1;

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

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-1.0, -1.0, -global.camera.distance);

  glRotatef(global.camera.rotation_x, 1.0, 0.0, 0.0);
  glRotatef(global.camera.rotation_y, 0.0, 1.0, 0.0);
  glRotatef(global.camera.rotation_z, 0.0, 0.0, 1.0);
}

void display(void) {
  double t, alpha;
  double clip_plane1[]={0.0,1.0,0.0,-h0/2.0};
  double clip_plane2[]={0.0,1.0,0.0,h0/2.0};

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glColor3ub(255,255,0);

  set_camera();
  renderAxes();
  
  GLUquadricObj *quadObj1;

  quadObj1 = gluNewQuadric();

  glColor4f(0.5,0.5,0.5, 1.0);
  glPushMatrix();
    glTranslatef(-a*1.5, 0, 0);
    glScalef(a, h0, l);
    glutSolidCube(1.0); 
    glTranslatef(3*a, 0, 0);
    glutSolidCube(1.0); 
  glPopMatrix();

  glClipPlane(GL_CLIP_PLANE1,clip_plane1);
  glEnable(GL_CLIP_PLANE1);

  glColor4f(0.5,0.5,0.5, 1.0);

  if (cin.fail()) return;

  cin >> t >> alpha;
  // cerr << po(alpha) << " " << center(alpha) << endl;
  

  glPushMatrix();
    glTranslatef(0, center(alpha), -l/2.0);
    gluCylinder(quadObj1, po(alpha), po(alpha), l, 100, 100);
  glPopMatrix();

  glClipPlane(GL_CLIP_PLANE1,clip_plane2);

  glPushMatrix();
    glTranslatef(0, center(alpha), -l/2.0);
    gluCylinder(quadObj1, po(alpha)-h(alpha), po(alpha)-h(alpha), l, 100, 100);
  glPopMatrix();

  glPushMatrix();
    glTranslatef(0, center(alpha), l/2.0);
    gluDisk(quadObj1, po(alpha) - h(alpha), po(alpha), 100, 100);
  glPopMatrix();


  glPushMatrix();
    glTranslatef(0, center(alpha), -l/2.0);
    gluDisk(quadObj1, po(alpha) - h(alpha), po(alpha), 100, 100);
  glPopMatrix();

  // usleep(1000);
  glDisable(GL_CLIP_PLANE1);
  
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