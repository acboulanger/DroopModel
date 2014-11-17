/**
 * \file droop0D.c
 * \brief Programme de test.
 * \author AC Boulanger
 * \version 0.1
 * \date 20/05/2010
 *
 * Test programm of Droop model in 0D - use of data from
 * paper "A mechanistic investigation of the algae growth "Droop" model" 
 * V. lemesle - L. Mailleret
 */



#include <stdio.h>
#include <stdlib.h>

/**
 * \struct Model
 * \brief Contains the variables of the Droop model
 *
 */

struct Model
{
  double x;/*!< Algae mass*/
  double q;/*!< Intra cellular quota*/
  double s;/*!< Nutrient mass*/
};

/*void init(struct Model *model,double x0, double q0, double s0);
double mu(double q, double qm);
double rho(double s);
double f1(double xn, double qn, double sn, double d, double sin);
double f2(double xn, double qn, double sn, double qm);
double f3(double xn, double qn, double sn, double qm, double d);
void simulation(struct Model *model, double dt, double d, double sin, double qm);
void exportGnuplot(const char *filename, double *sn, double *qn, double *xn, int taille);*/


/**
  \fn void exportGnuplot(const char *filename, double *sn, double *qn, double *xn, int taille)
  \brief Exportation of the data for Gnuplot
  \param filename : name of the txt file for data storage
  \param sn : array for s storage
  \param qn : array for q storage
  \param xn : array for x storage
  \param taille : size of previous arrays
  */
void exportGnuplot(const char *filename, double *sn, double *qn, double *xn, int taille)
{
  FILE* fichier = NULL;
  fichier = fopen(filename, "w+");
  int i;
  if (fichier != NULL)
  {
      //fprintf(fichier, "#s \t q \t x \n");
      for(i=0;i<taille;i++)
      {
	fprintf(fichier, "%f \t %f \t %f \n", sn[i], qn[i], xn[i]);
      }
      fclose(fichier);
  }
  else
  {
      printf("Could not open the file.");
  }
}



/**
  \fn void init(struct Model *model,double x0, double q0, double s0)
  \brief Initialisation of the model
  \param model : pointer on the model structure
  \param x0 : initial condition on x
  \param q0 : initial condition on q
  \param s0 : initial condition on s
 */
void init(struct Model *model,double x0, double q0, double s0)
{
  (*model).x=x0;
  (*model).q=q0;
  (*model).s=s0;
}



/**
  \fn double mu(double q, double qm)
  \brief Computation of the growth rate
  \param q : current intra cellular quota
  \param qm : subsistence quota
  */
double mu(double q, double qm)
{
  return 1.7*(1-qm/q);
}


/**
  \fn double rho(double s)
  \brief Computation of the nutrient uptake rate
  \param s : current nutrient mass
  */
double rho(double s, double q, double rhobarre, double ks, double ql)
{
  return (rhobarre*s/(s+ks)*(1-q/ql));
}


/**
  \fn double f1(double xn, double qn, double sn, double d, double sin)
  \brief Update function of s
  \param xn : current algae mass
  \param qn : current intra cellular quota
  \param sn : current nutrient mass
  \param d : diffusion coefficient
  \param sin : nutrient inflow
*/
double f1(double xn, double qn, double sn, double rhobarre, double ks, double ql)
{
  double res = 0;
  res = -rho(sn, qn, rhobarre, ks, ql)*xn;
  return res;
}


/**
  \fn double f2(double xn, double qn, double sn, double qm)
  \brief Update function of q
  \param xn : current algae mass
  \param qn : current intra cellular quota
  \param sn : current nutrient mass
  \param qm : subsistence quota
*/
double f2(double xn, double qn, double sn, double qm, double rhobarre, double ks, double ql)
{
  double res = 0;
  res = rho(sn, qn, rhobarre, ks, ql) - mu(qn,qm)*qn;
  return res;
}



/**
  \fn double f3(double xn, double qn, double sn, double qm, double d)
  \brief Update function of x
  \param xn : current algae mass
  \param qn : current intra cellular quota
  \param sn : current nutrient mass
  \param qm : subsistence quota
  \param d : diffusion coefficient
*/
double f3(double xn, double qn, double sn, double qm, double d)
{
  double res = 0;
  res = mu(qn,qm)*xn - d*xn;
  return res;
}


/**
  \fn void simulation(struct Model *model, double dt, double d, double sin, double qm)
  \brief One step simulation with Heun method
  \param model : pointer on the model structure
  \param dt : time step
  \param d : diffusion coefficient
  \param sin : nutrient inflow
  \param qm : subsistence quota
 */
void simulation(struct Model *model, double dt, double d, double sin, double qm, double rhobarre, double ks, double ql)
{
  double stilde;
  double qtilde;
  double xtilde;
  
  // Computation of the first step
  stilde = (*model).s + dt*f1((*model).x, (*model).q, (*model).s, rhobarre, ks, ql);
  qtilde = (*model).q + dt*f2((*model).x, (*model).q, (*model).s, qm, rhobarre, ks, ql);
  xtilde = (*model).x + dt*f3((*model).x, (*model).q, (*model).s, qm,d);
  
  //Computation of the second step
  (*model).s = (*model).s + dt/2*(f1((*model).x, (*model).q, (*model).s, rhobarre, ks, ql) + f1(xtilde, qtilde, stilde, rhobarre, ks, ql));
  (*model).q = (*model).q + dt/2*(f2((*model).x, (*model).q, (*model).s, qm, rhobarre, ks, ql) + f2(xtilde, qtilde, stilde, qm, rhobarre, ks, ql));
  (*model).x = (*model).x + dt/2*(f3((*model).x, (*model).q, (*model).s, qm,d) + f3(xtilde, qtilde, stilde, qm,d));
}


/**
 * \fn int main (void)
 * \brief Program entry.
 *
 * \return EXIT_SUCCESS - Arrêt normal du programme.
 */
int main()
{
  printf("Definition of the parameters...");
  /** Parameters of the simulation */
  double dt = 0.1;				//Time step
  int nbStep = 10000; 				//Number of iterations
  double *x_tab = NULL;				//Storage of x for each time step
  double *s_tab = NULL;				//Storage of s for each time step
  double *q_tab = NULL;				//Storage of q for each time step
  x_tab = malloc(nbStep * sizeof(double));
  s_tab = malloc(nbStep * sizeof(double));
  q_tab = malloc(nbStep * sizeof(double));
  
  char *filename = "drooplux.txt";   //File where data are stored s q x
  
  double x0 = 0.01; 				//Initial condition x
  double q0 = 0.5;				//Initial condition q
  double s0 = 0.00001;				//Initial condition s
  
  double sin = 10;				//Inflow of nutrients
  double d = 0.01;			//Dilution coefficient
  double qm = 0.05;				//Subsistence quota
  double rhobarre = 0.073;
  double ks = 0.0012;
  double ql = 0.25;
  printf("OK \n");
  
  
  printf("Initialisation...");
  /** Initialisation */
  struct Model model;
  init(&model,x0,q0,s0);
  printf("OK \n");

  printf("Time loop...");
  /** Time loop */
  int i;
  for(i = 0;i < nbStep;i++) // for each time step
  {
    simulation(&model,dt,d,sin,qm,rhobarre,ks,ql);
    x_tab[i] = model.x;
    q_tab[i] = model.q;
    s_tab[i] = model.s;
  }//end for
  printf("OK \n");
  
  /** Export Gnuplot */
  printf("Export Gnuplot...");
  exportGnuplot(filename ,s_tab, q_tab, x_tab, nbStep); 
  printf("OK \n");
  
  
  return EXIT_SUCCESS;
}