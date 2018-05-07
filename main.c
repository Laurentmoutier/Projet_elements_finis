//projet 
/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat    
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>



int main(void)
{  
	//ajout du devoir 4 
	int    n = 15;
	double radius = 0.1;
	double mass = 0.1;
	double radiusIn = 0.5;
	double radiusOut = 2.0;
	double dt = 1e-1;
	double tEnd = 8.0;
	double tol = 1e-6;
	double t = 0;
	double iterMax = 100;
	femGrains* theGrains = femGrainsCreateSimple(n, radius, mass, radiusIn, radiusOut);
    //fin de l'ajout du devoir 4

	printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    F-B-I : Full solver - Band solver - Iterative solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    char meshFileName[] = "../data/triangles_166_xopt.txt";
    
    // Pour Windows, remplacer l'argument :
    // ("../data/rect_quad_1601.txt") 
    // par :
    // ("..\\data\\rect_quad_1601.txt") 
    //
    // Sorry for the inconvenience :-)
    // On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !

    
    
    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    

    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("MECA1120 : homework 5 ");
    glfwMakeContextCurrent(window);

	//ajout du devoir 4
	int theRunningMode = 1.0;
	float theVelocityFactor = 0.25;
	//fin de l'ajout du devoir 4

    do 
    {
        int testConvergence,w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        glfwGetFramebufferSize(window,&w,&h);
      
        if (option == 1) {
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotField(theProblem->mesh,theProblem->soluce);   }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
        glColor3f(0.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
      
        if (solverType != newSolverType || renumType != newRenumType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
       
        glfwSwapBuffers(window);
        glfwPollEvents();

		//ajout du devoir 4:
		int i, w, h;
		double currentTime = glfwGetTime();

		glfwGetFramebufferSize(window, &w, &h);
		glfemReshapeWindows(radiusOut, w, h);
		for (i = 0; i < theGrains->n; i++) {
			glColor3f(1, 0, 0);
			glfemDrawDisk(theGrains->x[i], theGrains->y[i], theGrains->r[i]);
		}
		glColor3f(0, 0, 0); glfemDrawCircle(0, 0, radiusOut);
		glColor3f(0, 0, 0); glfemDrawCircle(0, 0, radiusIn);
		char theMessage[256];
		sprintf(theMessage, "Time = %g sec", t);
		glColor3f(1, 0, 0); glfemDrawMessage(20, 460, theMessage);
		glfwSwapBuffers(window);
		glfwPollEvents();

		if (t < tEnd && theRunningMode == 1) {
			printf("Time = %4g : ", t);
			//
			// A decommenter pour pouvoir progresser pas par pas
			//          printf("press CR to compute the next time step >>");
			//          char c= getchar();
			//
			femGrainsUpdate(theGrains, dt, tol, iterMax);
			t += dt;
		}

		while (glfwGetTime() - currentTime < theVelocityFactor) {
			if (glfwGetKey(window, 'R') == GLFW_PRESS)
				theRunningMode = 1;
			if (glfwGetKey(window, 'S') == GLFW_PRESS)
				theRunningMode = 0;
		}
		//fin de l'ajout du devoir 4

    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

