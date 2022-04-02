#include<GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include<math.h>
#include<bits/stdc++.h>
#include "../include/Variable.h"

using namespace std;


void light1()
{
    GLfloat no_light[]  = {0, 0, 0, 1.0};

    GLfloat light_ambient[]  = {.60, .60, .60, 1.0};
    GLfloat light_diffuse[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position[] = { -30,15,15, 1.0 };

    glLightfv( GL_LIGHT4, GL_POSITION, light_position);
    glLightfv( GL_LIGHT4, GL_AMBIENT, light_ambient);
    glLightfv( GL_LIGHT4, GL_DIFFUSE, light_diffuse);
    glLightfv( GL_LIGHT4, GL_SPECULAR, light_specular);
    glEnable( GL_LIGHT4);
}

void light2()
{
    GLfloat no_light[] = { 0.0, 0.0, 0.0, 1.0 };

    GLfloat global_ambient[] = { 0.0, 0.0, 0.0, 0.0 };

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,global_ambient);

    //Sunlight
    GLfloat light_ambient[]  = {.80, .80, .80, 1.0};
    GLfloat light_diffuse[]  = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position[] = { 20,50,80, 1.0 };

    glLightfv( GL_LIGHT0, GL_POSITION, light_position);

    //moonlight
    GLfloat light_ambient1[]  = { .40, .40, .40, 1.0 };
    GLfloat light_diffuse1[]  = { 80.0/255.0, 104.0/255.0, 134.0/255.0, 1.0 };
    GLfloat light_specular1[] = { .50, .50, .50, 1.0 };
    GLfloat light_position1[] = { 20.0, 50.0, 80.0, 1.0 };

    glLightfv( GL_LIGHT1, GL_POSITION, light_position1);

    GLfloat light_ambient2[]  = {0.5, 0.25, 0, 1.0};
    GLfloat light_diffuse2[]  = { 1.0, 0.5, 0, 1.0 };
    GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position2[] = { 7,-7,10, 1.0 };

    glLightfv( GL_LIGHT2, GL_POSITION, light_position2);

    // turning ambient light off or on
    if(amb == false)
    {
        glLightfv( GL_LIGHT0, GL_AMBIENT, no_light);
        glLightfv( GL_LIGHT1, GL_AMBIENT, no_light);
        glLightfv( GL_LIGHT2, GL_AMBIENT, no_light);
    }
    else
    {
        glLightfv( GL_LIGHT0, GL_AMBIENT, light_ambient);
        glLightfv( GL_LIGHT1, GL_AMBIENT, light_ambient1);
        glLightfv( GL_LIGHT2, GL_AMBIENT, light_ambient2);
    }

    // turning diffusion light off or on
    if(dif == false)
    {
        glLightfv( GL_LIGHT0, GL_DIFFUSE, no_light);
        glLightfv( GL_LIGHT1, GL_DIFFUSE, no_light);
        glLightfv( GL_LIGHT2, GL_DIFFUSE, no_light);
    }
    else
    {
        glLightfv( GL_LIGHT0, GL_DIFFUSE, light_diffuse);
        glLightfv( GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
        glLightfv( GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
    }

    // turning specular light off or on
    if(spec == true)
    {
        glLightfv( GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv( GL_LIGHT1, GL_SPECULAR, light_specular1);
        glLightfv( GL_LIGHT2, GL_SPECULAR, light_specular2);
    }
    else
    {
        glLightfv( GL_LIGHT0, GL_SPECULAR, no_light);
        glLightfv( GL_LIGHT1, GL_SPECULAR, no_light);
        glLightfv( GL_LIGHT2, GL_SPECULAR, no_light);
    }

    GLfloat spot_direction[] = { 0.0, 0.0, -1.0 };
    glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, spot_direction);
    glLightf( GL_LIGHT2, GL_SPOT_CUTOFF, 65);
    glLightf( GL_LIGHT2, GL_SPOT_EXPONENT, 2.0);



    if(Light1 == true && (disp == 4 || disp ==10 || disp == 11))
    {
        glDisable(GL_LIGHT4);
        glEnable( GL_LIGHT0);
    }
    else if(Light1 == false)
    {
        glDisable(GL_LIGHT0);
    }

    if(Light2 == true)
    {
        glEnable( GL_LIGHT1);
    }
    else if(Light2 == false)
    {
        glDisable(GL_LIGHT1);
    }

    if(Light3 == true)
    {
        glEnable( GL_LIGHT2);
    }
    else if(Light3 == false)
    {
        glDisable(GL_LIGHT2);
    }
}

