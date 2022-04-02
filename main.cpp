#include<GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include<math.h>
#include<bits/stdc++.h>
#include <vector>
#include<time.h>
#include "include/BmpLoader.h"
#include "include/light.h"

using namespace std;

unsigned int ID[100];

double Exval, Eyval, Ezval, Cxval, Cyval, Czval;
double rangle =0, rotval = 0, rotval2 = 0, handrot = 130.0, doll_y=0.0,doll_x = 0.0,ufo_z=0.0;
double Exval2, Eyval2, Ezval2, Cxval2, Cyval2, Czval2;
double windowHeight=700, windowWidth=1200;
GLboolean Light1 = true, Light2 = false, Light3 = false, amb = true, dif =true, spec =true, eyeFix = false;
int disp = 1, near_button = 0, button_press = 0, slide_h17= 0;
int game_won = 0,  wheel_rot = 0, black_rot = 1, blue_key = 0;
int hint = 0, hint_checked = 2, take_time= 1, hand =0,leg = 0;
double slide_x = 13, button_x =0;
double Play_time, passed_time, t_time, duration;
double count_down = 11;
int Best_score, Game_level, ending_animation=0;
int key_array[3] = {0,0,0};
double legrot=20;

const double PI = 3.14159265389;


/* GLUT callback Handlers */


int anglex= 0, angley = 0, anglez = 0;          //rotation angles
int window;
int wired=0;
int shcpt=0;
const int L1=18;
const int L2=19;
const int L3=10;
const int dgre=3;
int ncpt=L1+1;
const int nt = 100;				//number of slices along x-direction
const int nt1 = 10;				//number of slices along x-direction
const int ntheta = 100;
const int ntheta1 = 5;


static GLfloat v_cube[8][3] =
{
    {0.0,0.0,0.0},
    {0.0,1.0,0.0},
    {0.0,1.0,1.0},
    {0.0,0.0,1.0},
    {1.0,0.0,0.0},
    {1.0,0.0,1.0},
    {1.0,1.0,1.0},
    {1.0,1.0,0.0},
};

static GLfloat v_pyramid[5][3] =
{
    {0.0, 0.0, 0.0},
    {0.0, 0.0, 2.0},
    {2.0, 0.0, 2.0},
    {2.0, 0.0, 0.0},
    {1.0, 4.0, 1.0}
};

static GLfloat v_trapizium[8][3] =
{
    {0.0,0.0,0.0},
    {0.0,1.0,0.0},
    {-1.0,2.0,-1.0},
    {-1.0,-1.0,-1.0},
    {1.0,0.0,0.0},
    {2.0,-1.0,-1.0},
    {2.0,2.0,-1.0},
    {1.0,1.0,0.0},
};

static GLubyte p_Indices[4][3] =
{
    {4, 1, 2},
    {4, 2, 3},
    {4, 3, 0},
    {4, 0, 1}

};

static GLubyte q_Indices[6][4] =
{
    {0,3,2,1},
    {4,5,3,0},
    {0,1,7,4},
    {3,5,6,2},
    {4,7,6,5},
    {1,2,6,7}
};

GLfloat ctrlpoints1[L1+1][3] =
{
    {1.5,1.4,0},
    {2.0,2.8,0},
    {2.5,3.5,0},
    {3.0,4.2,0},
    {3.75,4.6,0},
    {3.5,4.6,0},
    {4.0,4.2,0},
    {4.5,2.8,0},
    {5,1.4,0},
    {5.5,0,0},
    {5,-1.4,0},
    {4.5,-2.8,0},
    {4.0,-4.2,0},
    {3.5,-4.6,0},
    {3.75,-4.6,0},
    {3.0,-4.2,0},
    {2.5,-3.5,0},
    {2.0,-2.8,0},
    {1.5,-1.4,0}
};

GLfloat ctrlpoints2[L2+1][3] =
{
    {1.5,1.4,0},
    {1.0,2.4,0},
    {0.5,2.4,0},
    {0.0,2.0,0},
    {-0.5,2.0,0},
    {-0.7,2.0,0},
    {-0.9,2.0,0},
    {-1.1,2.0,0},
    {-1.3,2.0,0},
    {-1.5,1.0,0},
    {-1.5,-1.0,0},
    {-1.3,-2.0,0},
    {-1.1,-2.0,0},
    {-0.9,-2.0,0},
    {-0.7,-2.0,0},
    {-0.5,-2.0,0},
    {0.0,-2.0,0},
    {0.5,-2.4,0},
    {1.0,-2.4,0},
    {1.5,-1.4,0}
};

GLfloat ctrlpoints3[L3+1][3] =
{
    {3.0,4,0},
    {3.5,4.6,0},
    {4.0,3.2,0},
    {4.5,2.8,0},
    {5,1.4,0},
    {5.5,0,0},
    {5,-1.4,0},
    {4.5,-2.8,0},
    {4.0,-3.2,0},
    {3.5,-4.6,0},
    {3.0,-4,0}
};

static void getNormal3p(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3, int neg = 0)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    if(neg == 1)
        glNormal3f(-Nx,-Ny,-Nz);
    else
        glNormal3f(Nx,Ny,Nz);
}

void material_prop(double dr, double dg, double db, int amb = 0, bool emi = false, int side = 0)
{
    GLfloat mat_ambient[] = { dr/2.0, dg/2.0, db/2.0, 1.0 };
    if(amb == 1)
        GLfloat mat_ambient[] = { dr, dg, db, 1.0 };
    GLfloat mat_diffuse[] = { dr, dg, db, 1.0 };
    GLfloat mat_specular[] = { .5, .5, .5, 1.0 };
    GLfloat mat_shininess[] = {60};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    if(side == 2)
    {
        glMaterialfv( GL_BACK, GL_AMBIENT, mat_ambient);
        glMaterialfv( GL_BACK, GL_DIFFUSE, mat_diffuse);
        glMaterialfv( GL_BACK, GL_SPECULAR, mat_specular);
        glMaterialfv( GL_BACK, GL_SHININESS, mat_shininess);
    }

    if(emi ==true)
    {
        GLfloat mat_emission[] = { dr, dg, db, 1.0 };
        glMaterialfv( GL_FRONT, GL_EMISSION, mat_emission);
    }
    else{
        GLfloat mat_emission[] = { 0, 0, 0, 0.0 };
        glMaterialfv( GL_FRONT, GL_EMISSION, mat_emission);
    }

}

void drawpyramid(GLfloat dr, GLfloat dg, GLfloat db)
{
    material_prop(dr,dg,db);
    //glColor3f(r,g,b);
    glBegin(GL_TRIANGLES);
    for (GLint i = 0; i <4; i++)
    {

        getNormal3p(v_pyramid[p_Indices[i][0]][0], v_pyramid[p_Indices[i][0]][1], v_pyramid[p_Indices[i][0]][2],
                    v_pyramid[p_Indices[i][1]][0], v_pyramid[p_Indices[i][1]][1], v_pyramid[p_Indices[i][1]][2],
                    v_pyramid[p_Indices[i][2]][0], v_pyramid[p_Indices[i][2]][1], v_pyramid[p_Indices[i][2]][2]);

        glVertex3fv(&v_pyramid[p_Indices[i][0]][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][1]][0]);
        glVertex3fv(&v_pyramid[p_Indices[i][2]][0]);
    }
    glEnd();

    glBegin(GL_QUADS);
    getNormal3p(v_pyramid[0][0], v_pyramid[0][1], v_pyramid[0][2],
                v_pyramid[3][0], v_pyramid[3][1], v_pyramid[3][2],
                v_pyramid[2][0], v_pyramid[2][1], v_pyramid[2][2]);

    glVertex3fv(&v_pyramid[0][0]);
    glVertex3fv(&v_pyramid[3][0]);
    glVertex3fv(&v_pyramid[2][0]);
    glVertex3fv(&v_pyramid[1][0]);
    glEnd();
}

void drawcube(double dr, double dg, double db)
{


    material_prop(dr,dg,db);

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {

        getNormal3p(v_cube[q_Indices[i][0]][0], v_cube[q_Indices[i][0]][1], v_cube[q_Indices[i][0]][2],
                    v_cube[q_Indices[i][1]][0], v_cube[q_Indices[i][1]][1], v_cube[q_Indices[i][1]][2],
                    v_cube[q_Indices[i][2]][0], v_cube[q_Indices[i][2]][1], v_cube[q_Indices[i][2]][2]);

        glVertex3fv(&v_cube[q_Indices[i][0]][0]);glTexCoord2f(1,1);
        glVertex3fv(&v_cube[q_Indices[i][1]][0]);glTexCoord2f(1,0);
        glVertex3fv(&v_cube[q_Indices[i][2]][0]);glTexCoord2f(0,0);
        glVertex3fv(&v_cube[q_Indices[i][3]][0]);glTexCoord2f(0,1);
    }
    glEnd();
}

void drawtrap(GLfloat dr, GLfloat dg, GLfloat db)
{
    material_prop(dr,dg,db);

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_trapizium[q_Indices[i][3]][0], v_trapizium[q_Indices[i][3]][1], v_trapizium[q_Indices[i][3]][2],
                    v_trapizium[q_Indices[i][2]][0], v_trapizium[q_Indices[i][2]][1], v_trapizium[q_Indices[i][2]][2],
                    v_trapizium[q_Indices[i][1]][0], v_trapizium[q_Indices[i][1]][1], v_trapizium[q_Indices[i][1]][2]);

        glVertex3fv(&v_trapizium[q_Indices[i][3]][0]);
        glVertex3fv(&v_trapizium[q_Indices[i][2]][0]);
        glVertex3fv(&v_trapizium[q_Indices[i][1]][0]);
        glVertex3fv(&v_trapizium[q_Indices[i][0]][0]);
    }
    glEnd();
}


void LoadTexture(const char*filename, int id)
{
    glGenTextures(100, &ID[id]);
    glBindTexture(GL_TEXTURE_2D, ID[id]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID[id]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}


void drawText(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, void *font, const char *textstring, int length, double dr, double dg, double db)
{
    glPushMatrix();
    material_prop(dr, dg, db, 1);

    glRasterPos3f(trans_x,trans_y,trans_z);
    for(int i=0; i< length; i++)
    {
        glutBitmapCharacter(font, (int)textstring[i]);
    }
    glPopMatrix();
}


void cube_obj(GLfloat scale_x, GLfloat scale_y, GLfloat scale_z, GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, int tex_id)
{
    glPushMatrix();
    if(tex_id >=1 && tex_id<=100)
    {
        glBindTexture(GL_TEXTURE_2D, ID[tex_id]);
        glEnable(GL_TEXTURE_2D);
    }

    glTranslatef(trans_x, trans_y, trans_z);
    glScalef(scale_x, scale_y, scale_z);
    if(tex_id >= 1 && tex_id <=100)
    {
        drawcube(1.0,1.0,1.0);
        glDisable(GL_TEXTURE_2D);
    }
    else if(tex_id == 0) // black
    {
        drawcube(0.0,0.0,0.0);
    }
    else if(tex_id == -1) // green
    {
        drawcube(0,154/255.0,23/255.0);
    }
    else if(tex_id == 101) // white
    {
        drawcube(1.0,1.0,1.0);
    }
    glPopMatrix();

}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2], int L, int cp)
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        if(cp == 1)
        {
            x+=coef*ctrlpoints1[i][0];
            y+=coef*ctrlpoints1[i][1];
        }
        else if(cp == 2)
        {
            x+=coef*ctrlpoints2[i][0];
            y+=coef*ctrlpoints2[i][1];
        }
        else if(cp == 3)
        {
            x+=coef*ctrlpoints3[i][0];
            y+=coef*ctrlpoints3[i][1];
        }

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

void Bezier(int L,int cp, GLfloat ctrlpoints[][3])
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t, xy, L,cp);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy, L,cp);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                getNormal3p(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z,1);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

void Bezier_nt1(int L,int cp, GLfloat ctrlpoints[][3])
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta1;		//angular step size

    float t=0;
    float dt=1.0/nt1;
    float xy[2];
    BezierCurve( t, xy, L,cp);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt1; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy, L,cp);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta1; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                getNormal3p(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z,1);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

void _button(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, int tex_id)
{
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z);
    glRotatef(-90,0,1,0);
    glScalef(1,1,1);

    GLUquadricObj *cyli;
    cyli = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[tex_id]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(cyli, true);
    gluCylinder(cyli, .15,.3,.3,20,20);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans_x-.3, trans_y, trans_z);
    glRotatef(-90,0,1,0);
    glScalef(1,1,1);

    GLUquadricObj *disk;
    disk = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[tex_id]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,.3,20,20);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

void _spotlight(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, double dr, double dg, double db)
{
    //holder
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z);
    glRotatef(-180,0,1,0);
    glScalef(1,1,1);

    material_prop(dr, dg, db);

    GLUquadricObj *cyli;
    cyli = gluNewQuadric();
    gluCylinder(cyli, 0,.3,.5,20,20);

    glPopMatrix();

    //bulb
    glPushMatrix();
    glTranslatef(trans_x , trans_y , trans_z - .3);
    glScalef(1,1,1);

    material_prop(255/255.0, 131/255.0, 0,0,Light3);

    GLUquadricObj *bulb;
    bulb = gluNewQuadric();
    gluSphere(bulb, 0.1,100,20);

    glPopMatrix();


}

void _color_wheel(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, int tex_id, int nonrot=0)
{
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z);
    if(nonrot == 0)
    {
        glRotatef(rotval,1,0,0);
    }
    glRotatef(-90,0,1,0);
    glScalef(1,1,1);

    GLUquadricObj *disk;
    disk = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[tex_id]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,.5,20,20);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

void blackhole(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, int tex_id, int side = 0)
{
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z);
    if(side == 0)
    {
        glRotatef(rotval2,0,1,0);
        glRotatef(-90,1,0,0);
    }
    else if(side == 1)
    {
        glRotatef(rotval2,1,0,0);
        glRotatef(-90,0,1,0);
    }
    glScalef(2,2,2.5);

    GLUquadricObj *disk;
    disk = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[tex_id]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,1.5,100,100);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

void doll(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z)
{

    //Eye left
    glPushMatrix();
    glTranslatef(trans_x-.2, trans_y-.54, trans_z+1.15);
    glRotatef(-120,1,0,0);
    glScalef(1,1,1);

    GLUquadricObj *disk;
    disk = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[18]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,.1,100,100);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    //eye right
    glPushMatrix();
    glTranslatef(trans_x+.2, trans_y-.54, trans_z+1.15);
    glRotatef(-120,1,0,0);
    glScalef(1,1,1);

    glBindTexture(GL_TEXTURE_2D, ID[18]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,.1,100,100);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();

    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

    glTranslated(trans_x,trans_y,trans_z);
    glRotatef( anglex, 1.0, 0.0, 0.0);
    glRotatef( angley, 0.0, 1.0, 0.0);         	//rotate about y-axis
    glRotatef( anglez, 0.0, 0.0, 1.0);

    // head
    glPushMatrix();
    material_prop(0/225.0,119/255.0,190/255.0,0,false,2);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(.3,.15,.15);
    Bezier(L1,1,ctrlpoints1);
    glPopMatrix();

    //body
    glPushMatrix();
    material_prop(0/225.0,119/255.0,190/255.0,0,false,2);
    glTranslatef(0,0,0);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(.4,.2,.2);
    Bezier(L2,2,ctrlpoints2);
    glPopMatrix();


    glPopMatrix();

    //mouth
    glPushMatrix();
    glTranslatef(trans_x , trans_y-.55 , trans_z+.9);
    glScalef(1,1,1);

    material_prop(153/255.0, 0/255.0, 0);

    GLUquadricObj *mouth;
    mouth = gluNewQuadric();
    gluSphere(mouth, 0.08,100,20);
    glPopMatrix();

    // hand left
    glPushMatrix();
    glTranslatef(trans_x-.25, trans_y-.1, trans_z+.5);
    glRotatef(-handrot,0,1,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    GLUquadricObj *cyli;
    cyli = gluNewQuadric();
    gluCylinder(cyli, .05,.05,.5,20,20);

    glPopMatrix();

    //hand right
    glPushMatrix();
    glTranslatef(trans_x+.25, trans_y-.1, trans_z+.5);
    glRotatef(handrot,0,1,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .05,.05,.5,20,20);

    glPopMatrix();

    //leg left
    glPushMatrix();
    glTranslatef(trans_x-.25, trans_y-.1, trans_z-.75);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .07,.07,.5,20,20);

    glPopMatrix();

    //leg right
    glPushMatrix();
    glTranslatef(trans_x+.25, trans_y-.1, trans_z-.75);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .07,.07,.5,20,20);

    glPopMatrix();

}

void doll2(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z)
{

    //Eye left
    glPushMatrix();
    glTranslatef(trans_x-.54, trans_y-.2, trans_z+1.15);
    glRotatef(120,0,1,0);
    glScalef(1,1,1);

    GLUquadricObj *disk;
    disk = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[18]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,.1,100,100);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    //eye right
    glPushMatrix();
    glTranslatef(trans_x-.54, trans_y+.2, trans_z+1.15);
    glRotatef(120,0,1,0);
    glScalef(1,1,1);

    glBindTexture(GL_TEXTURE_2D, ID[18]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(disk, true);
    gluDisk(disk,0,.1,100,100);
    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();

    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

    glTranslated(trans_x,trans_y,trans_z);
    glRotatef( anglex, 1.0, 0.0, 0.0);
    glRotatef( angley, 0.0, 1.0, 0.0);         	//rotate about y-axis
    glRotatef( anglez, 0.0, 0.0, 1.0);

    // head
    glPushMatrix();
    material_prop(0/225.0,119/255.0,190/255.0,0,false,2);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(.3,.15,.15);
    Bezier(L1,1,ctrlpoints1);
    glPopMatrix();

    //body
    glPushMatrix();
    material_prop(0/225.0,119/255.0,190/255.0,0,false,2);
    glTranslatef(0,0,0);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(.4,.2,.2);
    Bezier(L2,2,ctrlpoints2);
    glPopMatrix();


    glPopMatrix();

    //mouth
    glPushMatrix();
    glTranslatef(trans_x-.5 , trans_y , trans_z+.9);
    glScalef(1,1,1);

    material_prop(153/255.0, 0/255.0, 0);

    GLUquadricObj *mouth;
    mouth = gluNewQuadric();
    gluSphere(mouth, 0.08,100,20);
    glPopMatrix();

    // hand left
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y+.25, trans_z+.5);
    glRotatef(-130,1,0,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    GLUquadricObj *cyli;
    cyli = gluNewQuadric();
    gluCylinder(cyli, .05,.05,.5,20,20);

    glPopMatrix();

    //hand right
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y-.25, trans_z+.5);
    glRotatef(130,1,0,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .05,.05,.5,20,20);

    glPopMatrix();

    //leg left
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y+.25, trans_z-.85);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .07,.07,.8,20,20);

    glPopMatrix();

    //leg right
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y-.25, trans_z-.85);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .07,.07,.8,20,20);

    glPopMatrix();

}

void animated_doll(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z)
{
    glPushMatrix();

    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

    glTranslated(trans_x,trans_y,trans_z);

    // head
    glPushMatrix();
    material_prop(0/225.0,119/255.0,190/255.0,0,false,2);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(.3,.15,.15);
    Bezier(L1,1,ctrlpoints1);
    glPopMatrix();

    //body
    glPushMatrix();
    material_prop(0/225.0,119/255.0,190/255.0,0,false,2);
    glTranslatef(0,0,0);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(.4,.2,.2);
    Bezier(L2,2,ctrlpoints2);
    glPopMatrix();


    glPopMatrix();

    // hand left
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y+.25, trans_z+.5);
    glRotatef(-130,1,0,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    GLUquadricObj *cyli;
    cyli = gluNewQuadric();
    gluCylinder(cyli, .05,.05,.5,20,20);

    glPopMatrix();

    //hand right
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y-.25, trans_z+.5);
    glRotatef(130,1,0,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .05,.05,.5,20,20);

    glPopMatrix();

    //leg left
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y+.25, trans_z-.85);
    glRotatef(-legrot,0,1,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .07,.07,.8,20,20);

    glPopMatrix();

    //leg right
    glPushMatrix();
    glTranslatef(trans_x-.1, trans_y-.25, trans_z-.85);
    glRotatef(legrot,0,1,0);
    glScalef(1,1,1);

    material_prop(0/225.0,119/255.0,190/255.0);

    gluCylinder(cyli, .07,.07,.8,20,20);

    glPopMatrix();

}


void key(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, int color)
{
    glPushMatrix();
    // circular part
    glPushMatrix();
    if(color == 0)
    {
        material_prop(0.0,64/255.0,1.0);
        glTranslatef(trans_x,trans_y,trans_z);
        glRotatef(90,0,1,0);
        glScalef(.5,.5,.5);
        glutSolidTorus(0.2,0.8,50,50);
    }
    else if(color == 1)
    {
        glTranslatef(trans_x,trans_y,trans_z);
        glRotatef(90,0,1,0);
        glScalef(.3,.3,.5);
        drawtrap(1.0,36/255.0,0.0);
    }
    glPopMatrix();

    //key hole shape
    glPushMatrix();
    if(color == 0)
    {
        glTranslatef(trans_x-.25,trans_y+1.5,trans_z+.15);
        glRotatef(-90,1,0,0);
        glScalef(.25,.25,.25);
        drawpyramid(0.0,64/255.0,1.0);
    }
    else if(color == 1)
    {
        glTranslatef(trans_x+1.15,trans_y,trans_z-.05);
        glRotatef(-90,1,0,0);
        glScalef(.25,.25,.25);
        drawpyramid(1.0,36/255.0,0.0);
    }
    glPopMatrix();

    // key stick
    glPushMatrix();
    if(color == 0)
    {
        material_prop(0.0,64/255.0,1.0);
        glTranslatef(trans_x,trans_y+1.9,trans_z);
        glRotatef(90,1,0,0);
        GLUquadricObj *cyli;
        cyli = gluNewQuadric();
        gluCylinder(cyli, .15,.15,1.7,20,20);
    }
    else if(color == 1)
    {
        material_prop(1.0,36/255.0,0.0);
        glTranslatef(trans_x-.1,trans_y+.15,trans_z-.15);
        glRotatef(90,0,1,0);
        GLUquadricObj *cyli;
        cyli = gluNewQuadric();
        gluCylinder(cyli, .15,.15,1.7,20,20);
    }
    glPopMatrix();


    glPopMatrix();
}

void ufo(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z, double dr, double dg, double db)
{
    //lowerpart
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z);
    glScalef(1,1,1);

    material_prop(dr, dg, db);

    GLUquadricObj *cyli;
    cyli = gluNewQuadric();
    gluCylinder(cyli, 0,10,7,20,20);

    glPopMatrix();

    //middle part
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z+7);
    glScalef(1,1,1);


    material_prop(dr, dg, db);

    GLUquadricObj *disk;
    disk = gluNewQuadric();
    gluDisk(disk,0,10,100,100);
    glPopMatrix();

    //upper part
    glPushMatrix();
    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }
    glPushMatrix();
    material_prop(178/255.0, 151/255.0, 0/255.0,0,false,2);
    glTranslatef(trans_x, trans_y, trans_z-2.2);
    glRotatef( -90, 0.0, 1.0, 0.0);
    glScalef(3,1.5,1.5);
    Bezier_nt1(L3,3,ctrlpoints3);
    glPopMatrix();

    glPopMatrix();

}

void _maze_horizontal_walls(int level)
{
    //hori 1
    cube_obj(15,.5,9,8,8,0,4);

    //hori 2
    cube_obj(11,.5,9,30,8,0,4);

    //hori 3
    cube_obj(33,.5,9,8,-8,0,4);

    //hori 5
    cube_obj(14,.5,9,0,16,0,4);

    // hori 6
    cube_obj(13,.5,9,20,20,0,4);

    //hori 7
    cube_obj(13,1,9,28,28,0,4);

    //hori 8
    cube_obj(13,.5,9,20,36,0,4);

    //hori 9
    cube_obj(6,.5,9,27,48,0,4);

    //hori 10
    cube_obj(10,.5,9,23,56,0,4);

    //hori 11
    cube_obj(4,.5,9,19,46,0,4);

    //hori 12
    cube_obj(8,.5,9,19,40,0,4);

    //hori 13
    cube_obj(8,.5,9,33,78,0,4);

    //hori 14
    cube_obj(12,.5,9,21,60,0,4);

    //hori 15
    cube_obj(12,.5,9,21,68,0,4);

    //hori 18
    cube_obj(13,.5,9,13,88,0,4);

    //hori 19
    cube_obj(10,.5,9,20,96,0,4);

    //hori 21
    cube_obj(12,.5,9,21,78,0,4);

    //hori 24
    cube_obj(12,.5,9,0,30,0,4);

    //hori 25
    cube_obj(6,.5,9,0,38,0,4);

    //hori 26
    cube_obj(6,.5,9,0,52,0,4);

    //hori 27
    cube_obj(14.5,.5,9,-8,60,0,4);

    //hori 28
    cube_obj(12.5,.5,9,0,68,0,4);

    //hori 29
    cube_obj(13,.5,9,0,88,0,4);

    //hori 30
    cube_obj(20,.5,9,0,96,0,4);

    //hori 32
    cube_obj(33.5,.5,9,0,100,0,4);

    if(level == 1)
    {
        //hori 4
        cube_obj(7.5,.5,9,26,16,0,5);

        //hori 16
        cube_obj(8,.5,9,13,54,0,4);

        //hori 17
        cube_obj(8,.5,9,slide_x,74,0,6);

        //hori 20
        cube_obj(3.5,.5,9,30,88,0,5);

        //hori 22
        cube_obj(20,.5,9,21,108,0,5);

        //hori 23
        cube_obj(21,.5,9,-8,108,0,5);

        //hori 31
        cube_obj(8,.5,9,-8,98,0,6);
    }
    else if(level ==2)
    {
        //hori 4
        cube_obj(7.5,.5,9,26,16,0,4);

        //hori 20
        cube_obj(3.5,.5,9,30,88,0,4);

        //hori 33
        cube_obj(49,.5,9,-8,108,0,4);

        //hori 34
        cube_obj(16,.5,9,-8,-8,0,4);
    }
}


void _maze_vertical_walls(int level)
{
    //vert1
    cube_obj(.5,56,9,40.5,-8,0,4);

    //vert 2
    cube_obj(.5,100,9,-8,8,0,4);

    // vert 3
    cube_obj(.5,4,9,33,16.5,0,4);

    //vert 4
    cube_obj(.5,16,9,20,20,0,4);

    //vert 5
    cube_obj(.5,12.5,9,33,36,0,4);

    //vert 6
    cube_obj(.5,8,9,27,40,0,4);

    //vert 7
    cube_obj(.5,10,9,23,46,0,4);

    //vert 11
    cube_obj(.5,10,9,40.5,68,0,4);

    //vert 12
    cube_obj(.5,10,9,33,68,0,4);

    //vert 13
    cube_obj(.5,6,9,21,54,0,4);

    //vert 14
    cube_obj(.5,6,9,13,54,0,4);

    //vert 16
    cube_obj(.5,6,9,21,68,0,4);

    //vert 17
    cube_obj(.5,4,9,21,74,0,4);

    //vert 18
    cube_obj(.5,14,9,13,74,0,4);

    //vert 19
    cube_obj(.5,8,9,30,88,0,4);

    //vert 22
    cube_obj(.5,20,9,40.5,88,0,4);

    //vert 23
    cube_obj(.5,12,9,33,88.5,0,4);

    //vert 24
    cube_obj(.5,14,9,0,16,0,4);

    //vert 25
    cube_obj(.5,14,9,0,38,0,4);

    //vert 26
    cube_obj(.5,22,9,6,38,0,4);

    //vert 27
    cube_obj(.5,38,9,12,30,0,4);

    //vert 28
    cube_obj(.5,20,9,0,68,0,4);

    //vert 29
    cube_obj(.5,4,9,0,96,0,4);

    if(level == 1)
    {
        //vert 9
        cube_obj(.5,4.5,9,33,56,0,5);

        //vert 10
        cube_obj(.5,20,9,40.5,48,0,5);

        //vert 15
        cube_obj(.5,14,9,13,60,0,5);

        //vert 20
        cube_obj(.5,8,9,20,88,0,4);

        //vert 21
        cube_obj(.5,10,9,40.5,78,0,5);
    }
    else if(level == 2)
    {
        //vert 9
        cube_obj(.5,4.5,9,33,56,0,4);

        //vert 10
        cube_obj(.5,20,9,40.5,48,0,4);

        //vert 15
        cube_obj(.5,14,9,13,60,0,4);

        //vert 21
        cube_obj(.5,10,9,40.5,78,0,4);

        //vert 30
        cube_obj(.5,16,9,-8,-8,0,4);
    }
}

void maze1()
{
    glPushMatrix();

    glPushMatrix();
    //room wall back
    cube_obj(.5,16,9,-8,-8,0,1);

    //room wall front
    cube_obj(.5,7.5,9,8,-8,0,1); //right
    cube_obj(.5,3,2.5,8,-.5,6.5,14); //mid
    cube_obj(.5,6,9,8,2.5,0,1); //left

    //room wall right
    cube_obj(16,.5,9,-8,-8,0,1);

    //room wall left
    cube_obj(16,.5,9,-8,8,0,1);

    //ceiling
    cube_obj(18,18, .5, -9,-9,9,2);

    //floor
    cube_obj(50,117.5,.5,-9,-9,0,3);

    glPopMatrix();

    glPushMatrix();
    //horizontal walls
    _maze_horizontal_walls(1);

    //vertical walls
    _maze_vertical_walls(1);
    glPopMatrix();

    glPushMatrix();
    //black hole
    blackhole(16.5, 32, 4,15);
    glPopMatrix();

    glPushMatrix();
    //button
    _button(21.1+ button_x,73,5,8);
    glPopMatrix();

    glPushMatrix();
    //spotlight
    _spotlight(7,-7,9,33/255.0,33/255.0,105/255.0);
    glPopMatrix();

    glPushMatrix();
    // rotating colorwheel
    _color_wheel(40.4,12,5,13);
    _color_wheel(20.6,24,5,13);
    _color_wheel(40.4,32,5,13);
    _color_wheel(40.4,98,5,13);
    glPopMatrix();

    glPushMatrix();
    // non rotating colorwheel
    _color_wheel(-7.4,12,5,13,1);
    _color_wheel(-7.4,64,5,13,1);
    _color_wheel(11.9,34,5,13,1);
    glPopMatrix();

    //dialoue box for button
    if(near_button == 1)
    {
        glPushMatrix();
        cube_obj(.1,1,.8,20.9,71.2,4.65,0);
        glPopMatrix();

        glPushMatrix();
        std::string text;
        text = "Press X";

        drawText(20.8,71.8,5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
        glPopMatrix();
    }

    glPopMatrix();
}

void maze2()
{
    glPushMatrix();
    //animation part
    if(ending_animation == 1 && game_won == 1 && doll_x<=31.5)
    {
        glPushMatrix();
        animated_doll(-21.5-doll_x,50,1.5);
        glPushMatrix();
    }
    glPopMatrix();

    if(ufo_z<=60)
    {
        glPushMatrix();
        ufo(-55,50,.5+ufo_z,106/255.0,108/255.0,109/255.0);
        glPushMatrix();
    }

    glPopMatrix();

    glPushMatrix();
    //floor
    cube_obj(50,117.5,.5,-9,-9,0,3);
    glPopMatrix();

    glPushMatrix();
    //horizontal walls
    _maze_horizontal_walls(2);

    //vertical walls
    _maze_vertical_walls(2);
    glPopMatrix();

    glPushMatrix();
    //black hole1
    blackhole(16.5, 32, 4,15);
    glPopMatrix();

    glPushMatrix();
    //black hole2
    blackhole(-10, 50, 6,15,1);
    glPopMatrix();

    glPushMatrix();
    // rotating colorwheel
    _color_wheel(-7.4,12,5,13);
    _color_wheel(-7.4,64,5,13);
    _color_wheel(11.9,34,5,13);
    _color_wheel(40.4,12,5,13);
    _color_wheel(20.6,24,5,13);
    _color_wheel(40.4,32,5,13);
    _color_wheel(40.4,98,5,13);
    glPopMatrix();

    //red key appear, nullify blue
    glPushMatrix();
    if(Exval>=-8 && Exval<=0 && Eyval>=98.5 && Eyval<=101.5 && disp == 10 && key_array[1] == 0)
    {
        key_array[1] = 1;
        blue_key = 0;
    }
    else if(key_array[1] == 0 && disp == 10)
        key(-4,100,4,1);
    glPopMatrix();

    //blue key appear-reappear
    glPushMatrix();
    if(Exval>=18.5 && Exval<=21.5 && Eyval>=88 && Eyval<=96 && disp == 10 && blue_key==0)
    {
        blue_key = 1;
    }
    else if(blue_key==0 && disp ==10)
    {
        key(20,91.5,4,0);
    }
    glPopMatrix();

    //red key appear,nullifies the blue
    glPushMatrix();
    if(Exval>=13 && Exval<=21 && Eyval>=52.5 && Eyval<=55.5 && disp == 10 && key_array[2] == 0)
    {
        key_array[2] = 1;
        blue_key = 0;
    }
    else if(key_array[2] == 0 && disp == 10)
    {
        key(17,54,4,1);
    }
    glPopMatrix();

    glPopMatrix();
}


void maze3()
{
    //front_wall
    glPushMatrix();
    cube_obj(30,1,9,2,15,0,19);

    //right wall
    cube_obj(1,24,9,2,15,0,19);

    //left wall
    cube_obj(1,24,9,32,15,0,19);

    //back wall
    cube_obj(30,1,9,2,39,0,19);

    //floor
    cube_obj(31,25,.5,2,15,0,19);

    //ceiling
    cube_obj(31,25,.5,2,15,9,19);
    glPopMatrix();

    //blackholes
    glPushMatrix();
    blackhole(6.5,16.1, 4.0,15);
    blackhole(13.5,16.1, 4.0,15);
    blackhole(20.5,16.1, 4.0,15);
    blackhole(27.5,16.1, 4.0,15);
    glPopMatrix();

    // blackhole number
    glPushMatrix();
    std::string text;
    text = "4";
    drawText(6.5,16.2,4, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    glPopMatrix();

    glPushMatrix();
    text = "1";
    drawText(13.5,16.2,4, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    glPopMatrix();

    glPushMatrix();
    text = "2";
    drawText(20.5,16.2,4, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    glPopMatrix();

    glPushMatrix();
    text = "3";
    drawText(27.5,16.2,4, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    glPopMatrix();
}

void sky(int tex_id)
{
    glPushMatrix();
    cube_obj(1000,1000,800, -483.5,-450,-2,tex_id);
    glPopMatrix();
}

void tree(GLfloat trans_x, GLfloat trans_y, GLfloat trans_z)
{
    // tree trunk
    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z);
    glScalef(1,1,1);

    GLUquadricObj *tree;
    tree = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[11]);
    glEnable(GL_TEXTURE_2D);

    gluQuadricTexture(tree, true);
    gluCylinder(tree, .8,.2,8,20,20);

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    //leaves
    glPushMatrix();
    glTranslatef(trans_x-1, trans_y-1, trans_z+7.5);
    glScalef(1,1,1);

    GLUquadricObj *leaf;
    leaf = gluNewQuadric();
    glBindTexture(GL_TEXTURE_2D, ID[12]);
    glEnable(GL_TEXTURE_2D);

    gluQuadricTexture(leaf, true);
    gluSphere(leaf,2,100,20);

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans_x+1, trans_y-1, trans_z+7.5);
    glScalef(1,1,1);

    glBindTexture(GL_TEXTURE_2D, ID[12]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(leaf, true);
    gluSphere(leaf,2,100,20);

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans_x-1, trans_y+1, trans_z+7.5);
    glScalef(1,1,1);

    glBindTexture(GL_TEXTURE_2D, ID[12]);
    glEnable(GL_TEXTURE_2D);
    gluQuadricTexture(leaf, true);
    gluSphere(leaf,2,100,20);

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans_x+1, trans_y+1, trans_z+7.5);
    glScalef(1,1,1);

    gluQuadricTexture(leaf, true);
    glBindTexture(GL_TEXTURE_2D, ID[12]);
    glEnable(GL_TEXTURE_2D);

    gluSphere(leaf,2,100,20);

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(trans_x, trans_y, trans_z+9.5);
    glScalef(1,1,1);

    glBindTexture(GL_TEXTURE_2D, ID[12]);
    glEnable(GL_TEXTURE_2D);

    gluQuadricTexture(leaf, true);
    gluSphere(leaf,2,100,20);

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();


}

void surrounding()
{
    glPushMatrix();
    for(int i= -38; i<= 71;i+=10)
    {
         for (int j = 118; j<=140; j+=20)
            tree(i,j,0);
    }

    for(int i= -33; i<= 75;i+=10)
    {
         for (int j = 128; j<=150; j+=20)
            tree(i,j,0);
    }
    for(int i = -15; i<= 108; i+=15)
        tree(-15,i,0);

    for(int i = -15; i<= 108; i+=15)
        tree(48,i,0);

    for(int i = -5; i<= 38; i+=15)
        tree(i,-18,0);

    cube_obj(1000,1000,.5, -483.5,-450,-0.1,-1);
    glPopMatrix();
}

//score_board
void win_board(double elapsed_time)
{
    if(Game_level == 1)
    {
        glPushMatrix();
        cube_obj(8,.1,4.2,13,114,4.8,0);
        std::string text,time;
        int elap_time = (int)elapsed_time;
        time = to_string(elap_time);

        text = "Level ";
        drawText(16.8,113.5,8, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = to_string(Game_level);
        drawText(17.3,113.5,8, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = "You have completed in ";
        drawText(15.4,113.5,7, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        drawText(17.3,113.5,7, GLUT_BITMAP_TIMES_ROMAN_24, time.data(), time.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = " seconds";
        drawText(17.7,113.5,7, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = "Best Score: ";
        drawText(16.4,113.5,6.5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = to_string(Best_score);
        drawText(17.4,113.5,6.5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = " seconds";
        drawText(17.8,113.5,6.5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = "Press Enter to continue...";
        drawText(17.8,113.5,5.5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        glPopMatrix();
    }
    else if(Game_level == 2)
    {
        glPushMatrix();
        cube_obj(.1,8,5,-21,46.5,4,0);
        std::string text,time;
        int elap_time = (int)elapsed_time;
        time = to_string(elap_time);

        text = "Level ";
        drawText(-20.5,49.5,8, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = to_string(Game_level);
        drawText(-20.5,50.4,8, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = "You have completed in ";
        drawText(-20.5,48,7, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        drawText(-20.5,51,7, GLUT_BITMAP_TIMES_ROMAN_24, time.data(), time.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = " seconds";
        drawText(-20.5,51.5,7, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = "Best Score: ";
        drawText(-20.5,49,6, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = to_string(Best_score);
        drawText(-20.5,51,6, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = " seconds";
        drawText(-20.5,51.5,6, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        text = "Press Enter to go back home.";
        drawText(-20.5,49.5,5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),250.0/255.0,253.0/255.0,15.0/255.0);

        glPopMatrix();
    }

}

//starting page
void display1()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();
    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    cube_obj(.5,30,30,0,0,0,9);
    cube_obj(.5,10,5,-1,10,10,0);
    cube_obj(.5,10,5,-1.5,10,15,10);

    std::string text;
    text = "Press enter to Start!!!";
    drawText(-1.5,17.5,13, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//key options
void display2()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;
    text = "KEY OPTIONS";
    drawText(-1.5,17.5,24, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "E";
    drawText(-1.5,21.5,22, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Fixed/Unfixed your position";
    drawText(-1.5,15.5,22, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Up Arrow";
    drawText(-1.5,21.5,21, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Forward in North or South Direction";
    drawText(-1.5,15.5,21, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Down Arrow";
    drawText(-1.5,21.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Backward in North or South Direction";
    drawText(-1.5,15.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Left Arrow";
    drawText(-1.5,21.5,19, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = " ->  Forward in East or West Direction";
    drawText(-1.5,15.5,19, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Right Arrow";
    drawText(-1.5,21.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Backward in East or West Direction";
    drawText(-1.5,15.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "E then Up Arrow";
    drawText(-1.5,21.5,17, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Looking Upwards";
    drawText(-1.5,15.5,17, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "E then Down Arrow";
    drawText(-1.5,21.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Looking Downwards";
    drawText(-1.5,15.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "E then Left Arrow";
    drawText(-1.5,21.5,15, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Turning Left";
    drawText(-1.5,15.5,15, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "E then Right Arrow";
    drawText(-1.5,21.5,14, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Turning Right";
    drawText(-1.5,15.5,14, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "H";
    drawText(-1.5,21.5,13, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Show Map";
    drawText(-1.5,15.5,13, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "K";
    drawText(-1.5,21.5,12, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Show collected key";
    drawText(-1.5,15.5,12, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "1";
    drawText(-1.5,21.5,11, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Sunlight off/on";
    drawText(-1.5,15.5,11, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "2";
    drawText(-1.5,21.5,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    text = "->  Moonlight off/on";
    drawText(-1.5,15.5,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Press Enter to continue...";
    drawText(-1.5,10,7, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//Level 1 hints and ins
void display3()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "Level 1";
    drawText(-1.5,17.5,24, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Hints";
    drawText(-1.5,17.5,22, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "The wind is rising...";
    drawText(-1.5,20.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "The path is sparkling...";
    drawText(-1.5,20.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "The walls will guide you to your destination...";
    drawText(-1.5,20.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Instructions";
    drawText(-1.5,17.5,13, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Currently you are facing North.";
    drawText(-1.5,20.5,11, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "You can check the map atmost 2 times.";
    drawText(-1.5,20.5,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "You have 240 seconds to solve the maze.";
    drawText(-1.5,20.5,9, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Do not go to the blackhole unless you want to lost in between time and space.";
    drawText(-1.5,20.5,8, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Press Enter to continue...";
    drawText(-1.5,10,5, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//level 1
void display4()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    glPushMatrix();

    if(hint == 0)
    {
        gluLookAt(Exval, Eyval, Ezval, Cxval, Cyval, Czval, 0, 0, 1);
    }
    else
    {
        gluLookAt(Exval2, Eyval2, Ezval2, Cxval2, Cyval2, Czval2, 1, 0, 0);
        glPushMatrix();
        doll(Exval,Eyval,Ezval);
        glPopMatrix();
    }

    light2();

    if(game_won == 1)
    {
        if(take_time == 1)
        {
            Play_time = (t_time / 1000)- passed_time;
            passed_time = t_time/1000;
            Exval = 17;
            Eyval = 109;
            Cxval = 17;
            Cyval = 119;
            FILE *fptr;
            fptr = fopen("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\scores\\level_1.txt","w");
            if((int)Play_time< (int)Best_score)
            {
                Best_score = (int)Play_time;
            }
            fprintf(fptr,"%d", Best_score);
            fclose(fptr);
            take_time = 0;
        }

        glPushMatrix();
        win_board(int(Play_time));
        glPopMatrix();

        glPushMatrix();
        doll(17,115 + doll_y,3);
        glPopMatrix();

    }
    else{
        duration = (t_time / 1000)- passed_time;
        if(duration>= 240)
        {
            disp = 5;
        }
    }


    maze1();
    if(Light1 == true)
    {
        sky(7);
    }
    else{
        sky(17);
    }
    surrounding();


    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//time up
void display5()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "Level 1";
    drawText(-1.5,17.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Game Over!!!";
    drawText(-1.5,17.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Your time is up.";
    drawText(-1.5,18,17, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Press Enter to continue...";
    drawText(-1.5,10,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//inside the blackhole
void display6()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    cube_obj(.5,5000,3000,500,-2485,-1685,16);
    if(count_down>= 0)
    {
        std::string text;
        int remain_time = (int)count_down;
        text = to_string(remain_time);

        drawText(-1.5,17.5,15, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
        Sleep(1);
        count_down -= 0.1;
    }
    else{
        std::string text;

        text = "You are dead!!!";
        drawText(-1.5,17.5,15, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

        text = "Press Enter to continue...";
        drawText(-1.5,17.5,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);
    }

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}


//background story
void display7()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "Background Story";
    drawText(-1.5,17.5,24, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Once you went to take a tour of galaxy. While returning, ";
    drawText(-1.5,21.5,21, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "you forgot the path to your planet and landed on a different";
    drawText(-1.5,21.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "planet named Earth. The people of earth are scared of you. ";
    drawText(-1.5,21.5,19, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "So, they have caged you in a Labyrinth. While spending ";
    drawText(-1.5,21.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "time on Earth, you have gathered knowledge about earth and ";
    drawText(-1.5,21.5,17, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "learnt many languages from satellite using your power.";
    drawText(-1.5,21.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Press Enter to continue...";
    drawText(-1.5,10,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPushMatrix();
    doll2(-10+doll_y,15,12.5);
    glPopMatrix();

    glPopMatrix();


    glFlush();
    glutSwapBuffers();
}


//Level choose
void display8()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "Levels";
    drawText(-1.5,16.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Level 1 -> Press 1 to Play";
    drawText(-1.5,18.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Level 2 -> Press 2 to Play";
    drawText(-1.5,18.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//level 2 hints and ins
void display9()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "Level 2";
    drawText(-1.5,17.5,26, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "People have captured you again.";
    drawText(-1.5,19.5,25, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Hints";
    drawText(-1.5,17.5,23, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Blue is the ocean...";
    drawText(-1.5,19.5,22, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Blue is the sky...";
    drawText(-1.5,19.5,21, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Blue is the key...";
    drawText(-1.5,19.5,20, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "That will not let you die...";
    drawText(-1.5,19.5,19, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Numbers are just not numbers...";
    drawText(-1.5,20.5,17, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "The right number will give you freedom...";
    drawText(-1.5,20.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Instructions";
    drawText(-1.5,17.5,13, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Currently you are facing West.";
    drawText(-1.5,20.5,11, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "You can check the map at most 3 times inside the maze.";
    drawText(-1.5,20.5,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "You have 300 seconds to complete this level.";
    drawText(-1.5,20.5,9, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "The wrong key nullifies the right key.";
    drawText(-1.5,20.5,8, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Do not go to the blackhole without the right key.";
    drawText(-1.5,20.5,7, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Press Enter to continue...";
    drawText(-1.5,10,4, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//level 2 maze
void display10()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    glPushMatrix();
    if(hint == 0)
    {
        gluLookAt(Exval, Eyval, Ezval, Cxval, Cyval, Czval, 0, 0, 1);
    }
    else
    {
        gluLookAt(Exval2, Eyval2, Ezval2, Cxval2, Cyval2, Czval2, 1, 0, 0);
        glPushMatrix();
        doll(Exval,Eyval,Ezval);
        glPopMatrix();
    }


    light2();

    maze2();

    if(Light1 == true)
    {
        sky(7);
    }
    else{
        sky(17);
    }

    surrounding();

    if(game_won == 1)
    {
        if(take_time == 1)
        {
            //time interval
            Play_time = (t_time / 1000)- passed_time;
            passed_time = t_time/1000;

            //lookat fixed
            Exval = -13;
            Eyval = 50;
            Cxval = -26;
            Cyval = 50;

            //score write
            FILE *fptr;
            fptr = fopen("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\scores\\level_2.txt","w");
            if((int)Play_time< (int)Best_score)
            {
                Best_score = (int)Play_time;
            }
            fprintf(fptr,"%d", Best_score);
            fclose(fptr);

            //flag value change
            take_time = 0;
        }

        glPushMatrix();
        win_board(int(Play_time));
        glPopMatrix();
    }
    else{
        duration = (t_time / 1000)- passed_time;
        if(duration>= 300)
        {
            disp = 5;
        }
    }
    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}


// level 2 room
void display11()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    glPushMatrix();

    gluLookAt(Exval, Eyval, Ezval, Cxval, Cyval, Czval, 0, 0, 1);

    light2();

    maze3();


    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//Collected Key
void display12()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "Gathered Right key";
    drawText(-1.5,17.5,17, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = to_string(blue_key);
    drawText(-1.5,15.5,15, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

//The end
void display13()
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 90, w/h, 0.1, 2000);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    glViewport(0,0,w,h);

    light1();

    glPushMatrix();

    gluLookAt(-15, 15, 15, 5, 15, 15, 0, 0, 1);

    std::string text;

    text = "THE END";
    drawText(-1.5,17.5,18, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "You have reached your home safely";
    drawText(-1.5,20.5,16, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    text = "Press ESC to exit";
    drawText(-1.5,10,10, GLUT_BITMAP_TIMES_ROMAN_24, text.data(), text.size(),210.0/255.0,191.0/255.0,55.0/255.0);

    glPopMatrix();

    glFlush();
    glutSwapBuffers();
}


void display()
{
    t_time = glutGet(GLUT_ELAPSED_TIME);
    if(disp == 1)
    {
        display1();
    }
    else if(disp == 2)
    {
        display2();
    }
    else if(disp == 3)
    {
        display3();
    }
    else if(disp == 4)
    {
        display4();
    }
    else if(disp == 5)
    {
        display5();
    }
    else if(disp == 6)
    {
        display6();
    }
    else if(disp == 7)
    {
        display7();
    }
    else if(disp == 8)
    {
        display8();
    }
    else if(disp == 9)
    {
        display9();
    }
    else if(disp == 10)
    {
        display10();
    }
    else if(disp == 11)
    {
        display11();
    }
    else if(disp == 12)
    {
        display12();
    }
    else if(disp == 13)
    {
        display13();
    }
}

int isObstacle(double Ex, double Ey)
{
    if(Game_level == 1)
    {
        if((Ex <= -7.4 && Ex>= -9.0) || (Ex>=40.3)|| (Ex >=-8 && Ex <=23.1 && Ey >=7.9 && Ey<=8.6) ||
           (Ex >=7.9 && Ex <=8.6 && Ey >=-8 && Ey<=-0.6) || (Ex >=7.9 && Ex <=8.6 && Ey >=2.4 && Ey<=8.1) ||
           (Ex >=29.1 && Ex <=40.5 && Ey >=7.9 && Ey<=8.6) || (Ex >=-0.1 && Ex <=14.1 && Ey >=15.9 && Ey<=16.6) ||
           (Ex >=25.9 && Ex <=33.1 && Ey >=15.9 && Ey<=16.6) || (Ey<=-7.4) ||(Ex >=32.9 && Ex <=33.6 && Ey >=15.9 && Ey<=20.6) ||
           (Ex >=19.9 && Ex <=33.6 && Ey >=19.9 && Ey<=20.6) || (Ex >=27.9 && Ex <=41 && Ey >=27.9 && Ey<=28.6) ||
           (Ex >=19.9 && Ex <=20.6 && Ey >=19.9 && Ey<=36.6) || (Ex >=19.9 && Ex <=33.6 && Ey >=32.9 && Ey<=30.6) ||
           (Ex >=32.9 && Ex <=33.6 && Ey >=35.9 && Ey<=48.6) || (Ex >=26.9 && Ex <=33.6 && Ey >=47.9 && Ey<=48.6) ||
           (Ex >=26.9 && Ex <=27.6 && Ey >=39.9 && Ey<=48.6) || (Ex >=18.9 && Ex <=27.6 && Ey >=39.9 && Ey<=40.6) ||
           (Ex >=18.9 && Ex <=23.6 && Ey >=45.9 && Ey<=46.6) || (Ex >=22.9 && Ex <=23.6 && Ey >=45.9 && Ey<=56.6) ||
           (Ex >=22.9 && Ex <=33.6 && Ey >=53.9 && Ey<=56.6) || (Ex >=32.9 && Ex <=33.6 && Ey >=55.9 && Ey<=60.6) ||
           (Ex >=20.9 && Ex <=33.6 && Ey >=59.9 && Ey<=60.6) || (Ex >=20.9 && Ex <=33.6 && Ey >=67.9 && Ey<=68.6) ||
           (Ex >=32.9 && Ex <=33.6 && Ey >=67.9 && Ey<=78.6) || (Ex >=20.9 && Ex <=41.6 && Ey >=77.9 && Ey<=78.6) ||
           (Ex >=20.9 && Ex <=21.6 && Ey >=59.9 && Ey<=54.6) || (Ex >=12.9 && Ex <=21.6 && Ey >=53.9 && Ey<=54.6) ||
           (Ex >=12.9 && Ex <=13.6 && Ey >=53.9 && Ey<=88.6) || (Ex >=20.9 && Ex <=21.6 && Ey >=67.9 && Ey<=78.6) ||
           (Ex >=-0.1 && Ex <=26.1 && Ey >=87.9 && Ey<=88.6) || (Ex >=29.9 && Ex <=33.6 && Ey >=87.9 && Ey<=88.6) ||
           (Ex >=-0.1 && Ex <=30.6 && Ey >=95.9 && Ey<=97.6) || (Ex >=19.9 && Ex <=20.6 && Ey >=87.9 && Ey<=96.6) ||
           (Ex >=32.9 && Ex <=33.6 && Ey >=87.9 && Ey<=100.6) || (Ex >=-0.1 && Ex <=33.6 && Ey >=99.9 && Ey<=100.6) ||
           (Ex >=-7.4 && Ex <=13.1 && Ey >=107.9 && Ey<=108.6) || (Ex >=20.9 && Ex <=41 && Ey >=107.9 && Ey<=108.6) ||
           (Ex >=-0.1 && Ex <=1 && Ey >=15.9 && Ey<=30.6) || (Ex >=-0.1 && Ex <=12.6 && Ey >=29.9 && Ey<=30.6) ||
           (Ex >=-0.1 && Ex <=6.6 && Ey >=37.9 && Ey<=38.6) || (Ex >=-0.1 && Ex <=6.6 && Ey >=51.9 && Ey<=52.6) ||
           (Ex >=-7.4 && Ex <=6.6 && Ey >=59.9 && Ey<=60.6) || (Ex >=5.9 && Ex <=6.6 && Ey >=37.9 && Ey<=60.6) ||
           (Ex >=11.9 && Ex <=12.6 && Ey >=29.9 && Ey<=68.6) || (Ex >=-0.1 && Ex <=12.6 && Ey >=67.9 && Ey<=68.6) ||
           (Ex >=-0.1 && Ex <=1 && Ey >=67.9 && Ey<=88.6) || (Ex >=-7.4 && Ex <=0 && Ey >=97.9 && Ey<=98.6) ||
           (Ex >=-0.1 && Ex <=1 && Ey >=95.9 && Ey<=100.6)
           )
        {
            return 1;
        }
        else if((Ex >=12.9 && Ex <=21.6 && Ey >=73.5 && Ey<=75.6))
        {
            if(slide_h17 == 2)
                return 0;
            else
                return 1;
        }
        else return 0;
    }
    else if(Game_level == 2)
    {
        if(disp == 10)
        {
            if((Ex <= -7.4 && Ex>= -9.0) || (Ex>=40.3)|| (Ex >=7.9 && Ex <=23.1 && Ey >=7.9 && Ey<=8.6) ||
               (Ex >=29.1 && Ex <=40.5 && Ey >=7.9 && Ey<=8.6) || (Ex >=-0.1 && Ex <=14.1 && Ey >=15.9 && Ey<=16.6) ||
               (Ex >=25.9 && Ex <=33.1 && Ey >=15.9 && Ey<=16.6) || (Ey<=-7.4) ||(Ex >=32.9 && Ex <=33.6 && Ey >=15.9 && Ey<=20.6) ||
               (Ex >=19.9 && Ex <=33.6 && Ey >=19.9 && Ey<=20.6) || (Ex >=27.9 && Ex <=41 && Ey >=27.9 && Ey<=28.6) ||
               (Ex >=19.9 && Ex <=20.6 && Ey >=19.9 && Ey<=36.6) || (Ex >=19.9 && Ex <=33.6 && Ey >=32.9 && Ey<=30.6) ||
               (Ex >=32.9 && Ex <=33.6 && Ey >=35.9 && Ey<=48.6) || (Ex >=26.9 && Ex <=33.6 && Ey >=47.9 && Ey<=48.6) ||
               (Ex >=26.9 && Ex <=27.6 && Ey >=39.9 && Ey<=48.6) || (Ex >=18.9 && Ex <=27.6 && Ey >=39.9 && Ey<=40.6) ||
               (Ex >=18.9 && Ex <=23.6 && Ey >=45.9 && Ey<=46.6) || (Ex >=22.9 && Ex <=23.6 && Ey >=45.9 && Ey<=56.6) ||
               (Ex >=22.9 && Ex <=33.6 && Ey >=53.9 && Ey<=56.6) || (Ex >=32.9 && Ex <=33.6 && Ey >=55.9 && Ey<=60.6) ||
               (Ex >=20.9 && Ex <=33.6 && Ey >=59.9 && Ey<=60.6) || (Ex >=20.9 && Ex <=33.6 && Ey >=67.9 && Ey<=68.6) ||
               (Ex >=32.9 && Ex <=33.6 && Ey >=67.9 && Ey<=78.6) || (Ex >=20.9 && Ex <=41.6 && Ey >=77.9 && Ey<=78.6) ||
               (Ex >=20.9 && Ex <=21.6 && Ey >=59.9 && Ey<=54.6) || (Ex >=12.9 && Ex <=13.6 && Ey >=53.9 && Ey<=88.6) ||
               (Ex >=20.9 && Ex <=21.6 && Ey >=67.9 && Ey<=78.6) || (Ex >=-0.1 && Ex <=26.1 && Ey >=87.9 && Ey<=88.6) ||
               (Ex >=29.9 && Ex <=33.6 && Ey >=87.9 && Ey<=88.6) || (Ex >=-0.1 && Ex <=30.6 && Ey >=95.9 && Ey<=97.6) ||
               (Ex >=32.9 && Ex <=33.6 && Ey >=87.9 && Ey<=100.6) || (Ex >=-0.1 && Ex <=33.6 && Ey >=99.9 && Ey<=100.6) ||
               (Ex >=-7.4 && Ex <=13.1 && Ey >=107.9 && Ey<=108.6) || (Ex >=20.9 && Ex <=41 && Ey >=107.9 && Ey<=108.6) ||
               (Ex >=-0.1 && Ex <=1 && Ey >=15.9 && Ey<=30.6) || (Ex >=-0.1 && Ex <=12.6 && Ey >=29.9 && Ey<=30.6) ||
               (Ex >=-0.1 && Ex <=6.6 && Ey >=37.9 && Ey<=38.6) || (Ex >=-0.1 && Ex <=6.6 && Ey >=51.9 && Ey<=52.6) ||
               (Ex >=-7.4 && Ex <=6.6 && Ey >=59.9 && Ey<=60.6) || (Ex >=5.9 && Ex <=6.6 && Ey >=37.9 && Ey<=60.6) ||
               (Ex >=11.9 && Ex <=12.6 && Ey >=29.9 && Ey<=68.6) || (Ex >=-0.1 && Ex <=12.6 && Ey >=67.9 && Ey<=68.6) ||
               (Ex >=-0.1 && Ex <=1 && Ey >=67.9 && Ey<=88.6) || (Ex >=-0.1 && Ex <=1 && Ey >=95.9 && Ey<=100.6)
               )
            {
                return 1;
            }
            else return 0;
        }
        else if(disp == 11)
        {
            if(Ex >=31 || Ex <=3  ||  Ey >=38 || Ey<=15)
            {
                return 1;
            }
            else return 0;
        }
    }

}


void myKeyboardFunc( unsigned char key, int x, int y )
{
    switch ( key )
    {
    // enter key
    case 13:
        if(disp == 1)
        {
            //starting to key
            disp = 2;
        }
        else if(disp == 2)
        {
            // key to backstory
            doll_y= 0.0;
            disp = 7;
        }
        else if(disp == 7)
        {
            //to level choose
            disp = 8;
        }
        else if(disp == 3)
        {
            //hints to level 1
            disp = 4;
            passed_time = t_time /1000;
            Exval=-4;
            Eyval=1;
            Ezval = 5;
            Cxval=4;
            Cyval=1;
            Czval = 5;
            Exval2= 20;
            Eyval2=50;
            Ezval2 = 65;
            Cxval2=20;
            Cyval2=50;
            Czval2 = 0;
            eyeFix = false;
            handrot = 130.0;
            game_won = 0;
            wheel_rot = 0;
            doll_y = 0.0;
            hint = 0;
            hint_checked = 2;
            take_time = 1;

            FILE *fptr;
            fptr = fopen("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\scores\\level_1.txt","r");
            fscanf(fptr,"%d", &Best_score);
            fclose(fptr);
        }
        else if(disp==4)
        {
            //to level 2 hints
            if(game_won==1)
            {
                Game_level = 2;
                disp = 9;
            }
        }
        else if(disp == 5)
        {
            // back to hints level
            if(Game_level == 1)
                disp = 3;
            else if(Game_level == 2)
                disp = 9;
        }
        else if(disp == 6 && count_down<0)
        {
            //restart level
            if(Game_level == 1)
            {
                disp = 4;
                Exval=-4;
                Eyval=1;
                Ezval = 5;
                Cxval=4;
                Cyval=1;
                Czval = 5;
                count_down = 11;
                eyeFix = false;
                handrot =130;
            }
            else if(Game_level == 2)
            {
                disp = 10;
                Exval=26.5;
                Eyval=0;
                Ezval = 5;
                Cxval=26.5;
                Cyval=20;
                Czval = 5;
                key_array[1] = 0;
                key_array[2] = 0;
                count_down = 11;
                eyeFix = false;
            }
        }
        else if(disp == 9)
        {
            // to level 2
            disp = 10;
            passed_time = t_time /1000;
            Exval=26.5;
            Eyval=0;
            Ezval = 5;
            Cxval=26.5;
            Cyval=20;
            Czval = 5;
            Exval2= 20;
            Eyval2=50;
            Ezval2 = 65;
            Cxval2=20;
            Cyval2=50;
            Czval2 = 0;
            hint = 0;
            hint_checked = 3;
            take_time = 1;
            eyeFix = false;
            game_won = 0;
            black_rot = 1;
            wheel_rot = 1;
            blue_key = 0;
            ending_animation=0;

            FILE *fptr;
            fptr = fopen("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\scores\\level_2.txt","r");
            fscanf(fptr,"%d", &Best_score);
            fclose(fptr);
        }
        else if(disp == 10)
        {
            //animation start
            if(game_won == 1)
            {
                Exval = -22;
                Eyval = 50;
                Ezval = 5;
                Cxval = -35;
                Cyval = 50;
                Czval = 5;
                ending_animation = 1;
            }
        }
        break;

    //level choose, sun light off, moonlight on and sunlight on and moonlight off
    case 49:
        if(disp == 8)
        {
            disp = 3;
            Game_level = 1;
        }
        else if(Light1 == true && (disp == 4 || disp == 10 || disp ==11))
        {
            Light1 = false;
            Light2 = true;
            if(disp == 4)
                Light3 = true;
        }
        else if(Light1 == false && (disp == 4 || disp == 10 || disp == 11))
        {
            Light1 = true;
            Light2 = false;
            if(disp == 4)
                Light3 = false;
        }
        break;

    case 50:
        if(disp == 8)
        {
            disp = 9;
            Game_level = 2;
        }
        else if(Light2 == true && (disp == 4 || disp == 10))
        {
            Light2 = false;
        }
        else if(Light2 == false && (disp == 4 || disp == 10))
        {
            Light2 = true;
        }
        break;

    // spotlight on/off
    case 51:
        if(Light3 == true && disp == 4)
        {
            Light3 = false;
        }
        else if(Light3 == false && disp == 4)
        {
            Light3 = true;
        }
        break;


    //ambient off on
    case 'A':
    case 'a':
        if(amb == true)
            amb = false;
        else
            amb = true;
        break;

    //diffuse off on
    case 'D':
    case 'd':
        if(dif == true)
            dif = false;
        else
            dif = true;
        break;

    //specular off on
    case 'S':
    case 's':
        if(spec == true)
            spec = false;
        else
            spec = true;
        break;

    // eye fix or not
    case 'E':
    case 'e':
        if(eyeFix == false)
        {
            eyeFix = true;
        }
        else
        {
            eyeFix = false;
        }
        break;

    // hint
    case 'h':
    case 'H':
        if (hint == 0 && hint_checked > 0)
        {
            hint = 1;
            hint_checked -= 1;
        }
        else
        {
            hint = 0;
        }
        break;

    //button press
    case 'x':
    case 'X':
        if(near_button == 1 && disp == 4)
        {
            button_press = 1;
        }
        else
        {
            button_press = 0;
        }
        break;

    // show key
    case 'k':
    case 'K':
        if(disp == 10)
        {
            disp = 12;
        }
        else if(disp == 12)
        {
            disp = 10;
        }

    // wired version
    case 'w':
    case 'W':
        wired=!wired;
        break;

    // angular rotate doll body
    case 't':
        anglex = ( anglex + 3 ) % 360;
        break;
    case 'T':
        anglex = ( anglex - 3 ) % 360;
        break;

    case 'y':
        angley = ( angley + 3 ) % 360;
        break;
    case 'Y':
        angley = ( angley - 3 ) % 360;
        break;

    case 'z':
        anglez = ( anglez + 3 ) % 360;
        break;
    case 'Z':
        anglez = ( anglez - 3 ) % 360;
        break;

    //zoom in doll
    case '+':
        doll_y -= 0.5;
        break;

    // zoom out doll
    case '-':
        doll_y +=.5;
        break;

    case 27:	// Escape key
        exit(1);
        break;
    }
}

void mySpecialKeyFunc(int key, int x, int y)
{
    switch(key)
    {
    // rotate left without changing eye position / move forward in y-axis
    case GLUT_KEY_LEFT:
        if(eyeFix == true && game_won == 0)
        {
            //camera turning
            double newCXval, newCYval;
            rangle += 0.1;
            newCXval = Cxval*cos(rangle) - Cyval*sin(rangle);
            newCYval = Cxval*sin(rangle) + Cyval*cos(rangle);
            Cxval = newCXval ;
            Cyval = newCYval;
            rangle = 0;
        }
        else if(eyeFix == false && game_won == 0)
        {
            double newEYval = Eyval + Cyval*0.01;
            double newEXval = Exval + Cxval*0.01;
            double newCYval = Cyval + Cyval*0.01;
            double newCXval = Cxval + Cxval*0.01;

            if(isObstacle(newEXval, newEYval) == 0)
            {
                Eyval = newEYval;
                Cyval = newCYval;
            }
        }
        break;

    // turn right without changing eye position / move backward in y axis
    case GLUT_KEY_RIGHT:
        if(eyeFix == true && game_won == 0)
        {
            double newCXval, newCYval;
            rangle -= 0.1;
            newCXval = Cxval*cos(rangle) - Cyval*sin(rangle);
            newCYval = Cxval*sin(rangle) + Cyval*cos(rangle);
            Cxval = newCXval;
            Cyval = newCYval;
            rangle = 0;
        }
        else if(eyeFix == false && game_won == 0)
        {
            double newEYval = Eyval - Cyval*0.01;
            double newEXval = Exval - Cxval*0.01;
            double newCYval = Cyval - Cyval*0.01;
            double newCXval = Cxval - Cxval*0.01;

            if(isObstacle(newEXval, newEYval) == 0)
            {
                Eyval = newEYval;
                Cyval = newCYval;
            }
        }
        break;

    // looking downwards / move backward in x axis
    case GLUT_KEY_DOWN:
        if(eyeFix == true && game_won == 0)
        {
            Czval-=0.2;
        }
        else if(eyeFix == false && game_won == 0 )
        {
            double newEYval = Eyval - Cyval*0.01;
            double newEXval = Exval - Cxval*0.01;
            double newCYval = Cyval - Cyval*0.01;
            double newCXval = Cxval - Cxval*0.01;

            if(isObstacle(newEXval, newEYval) == 0)
            {
                Exval = newEXval;
                Cxval = newCXval;
            }
        }
        break;

    // looking upward / move forward in x-axis
    case GLUT_KEY_UP:
        if(eyeFix == true && game_won == 0 )
        {
            Czval+=0.2;
        }
        else if(eyeFix == false && game_won == 0)
        {
            double newEYval = Eyval + Cyval*0.01;
            double newEXval = Exval + Cxval*0.01;
            double newCYval = Cyval + Cyval*0.01;
            double newCXval = Cxval + Cxval*0.01;

            if(isObstacle(newEXval, newEYval) == 0)
            {
                Exval = newEXval;
                Cxval = newCXval;
            }
        }
        break;
    }
}


void animate()
{
    //doll walking to ufo and ufo taking off
    if(ending_animation == 1 && game_won == 1 && disp == 10)
    {
        if(doll_x<=32.5)
        {
            doll_x+= 0.1;
        }
        else
        {
            if(ufo_z<=60)
            {
                ufo_z+= 0.5;
            }
            else
            {
                disp = 13;
                ending_animation = 0;
            }
        }

        if(leg == 0)
        {
            legrot -= 5;
            if(legrot<=-20)
                leg = 1;
        }
        else if(leg == 1)
        {
            legrot += 5;
            if(legrot >=20)
                leg = 0;
        }
    }


    //doll hand
    if(game_won == 1 && disp == 4)
    {
        if(hand == 0)
        {
            handrot -= 5.0;
            if(handrot<=80.0)
                hand = 1;
        }
        if(hand== 1)
        {
            handrot += 5.0;
            if(handrot>=130)
                hand = 0;
        }
    }
    //wall button
    if((Exval>=13 && Exval<=21) && (Eyval >=67.5 && Eyval<=73.5) && Ezval == 5 && disp == 4)
    {
        near_button = 1;
    }
    else
    {
        near_button = 0;
    }

    //slide the wall
    if(button_press == 1 && disp == 4 && (Exval>=13 && Exval<=21) && (Eyval >=67.5 && Eyval<=73.5) && Ezval == 5 )
    {
        if(button_x<= 0.1)
        {
            button_x+=.01;
            slide_h17 = 0;
        }
        else{
            slide_h17 = 1;
        }
    }

    if(slide_h17 == 1)
    {
        if(slide_x<=21.5)
        {
            slide_x += .1;
        }
        else
        {
            slide_h17 = 2;
        }
    }

    //win
    if((Exval>=13 && Exval<=21) && (Eyval >=109 )&& disp == 4)
    {
        game_won = 1;
    }
    else if((Exval <=-13) && disp == 10)
    {
        game_won = 1;
    }
    else
    {
        game_won = 0;
    }


    // dynamic wheel
    if((Eyval>=8 && Eyval<= 16) && (Exval>=32 && Exval<=41) && disp == 4)
    {
        wheel_rot = 1;
    }
    else if((Eyval>=8 && Eyval<= 16) && (Exval<32) && disp == 4)
    {
        wheel_rot = 0;
    }

    if(wheel_rot == 1 && (disp == 4 || disp ==10))
    {
        if(rotval<360)
        {
            rotval +=30;
        }
        else
        {
            rotval = rotval - 360.0;
        }
    }

    //blackhole
    if(black_rot == 1 && (disp ==4 || disp == 10 || disp == 11))
    {
        if(rotval2<360)
        {
            rotval2 += 30;
        }
        else
        {
            rotval2 = rotval2 - 360.0;
        }
    }

    //inside blackhole of starting maze
    if((Exval>=12.6 && Exval<=19.8) && (Eyval>=29 && Eyval <= 35))
    {
        if(disp == 10 && blue_key == 0)
            disp = 6;
        else if(disp == 10 && blue_key == 1)
        {
            disp =11;
        }
        else if(disp == 4)
            disp = 6;
    }

    // inside the blackholes of the room
    if((Exval>=3.5 && Exval<=9.5) && (Eyval<=19) && disp == 11)
    {
        disp = 10;
        Exval = -12;
        Eyval = 50;
        Ezval = 5;
        Cxval = -25;
        Cyval = 50;
        Czval = 5;
        hint_checked = 0;
    }
    else if((Exval>=10.5 && Exval<=16.5) && (Eyval<=17.5) && disp == 11)
    {
        disp = 6;
    }
    else if((Exval>=24.5 && Exval<=30.5) && (Eyval<=17.5) && disp == 11)
    {
        disp = 6;
    }
    if((Exval>=17.5 && Exval<=23.5) && (Eyval<=19) && disp == 11)
    {
        disp = 10;
        Exval = -11;
        Eyval = 50;
        Ezval = 5;
        Cxval = -25;
        Cyval = 50;
        Czval = 5;
        hint_checked = 0;
    }

    glutPostRedisplay();

}


int main (int argc, char **argv)
{

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowPosition(50,0);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("LABYRINTH");


    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\Cement_wall.bmp",1);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\ceiling-texture2.bmp",2);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\soil.bmp",3);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\green-grass-textures.bmp",4);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\grass_flower.bmp",5);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\green-grass_flower.bmp",6);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\sky2.bmp",7);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\rose.bmp",8);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\LABYRINTH.bmp",9);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\title2.bmp",10);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\tree-texture.bmp",11);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\leaf.bmp",12);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\color_wheel.bmp",13);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\Cement_wall3.bmp",14);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\black-hole.bmp",15);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\space.bmp",16);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\night_sky.bmp",17);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\eye.bmp",18);
    LoadTexture("F:\\class_lecture\\4-2\\CSE 4208 Graphics lab\\Labyrinth\\textures\\pink_wall.bmp",19);

    glClearColor(0,0,0,1);


    glShadeModel( GL_SMOOTH );
    glEnable( GL_DEPTH_TEST );
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);

    glutSpecialFunc (mySpecialKeyFunc);
    glutKeyboardFunc(myKeyboardFunc);
    glutDisplayFunc(display);
    glutIdleFunc(animate);

    glutMainLoop();

    return 0;
}

