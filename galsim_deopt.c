#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "graphics.h"

struct particle
{
  double x;
  double y;
  double mass;
  double vx;
  double vy;
  //double brightness;
};
typedef struct particle particle_t;

struct quad
{
  double xstart; //the x coordinate of the top left corner
  double ystart; //the y coordinate of the top left corner
  double size; //the width / height of the quadrant
  double x; //the x position of the centre of mass of this quadrant
  double y; //the y position of the centre of mass of this quadrant
  double m; //the total mass of this quadrant
  particle_t* particle; //If there is only one particle, this is a pointer to it
  struct quad* topleft; //a pointer to the topleft sub-quadrant
  struct quad* topright; //a pointer to the topright sub-quadrant
  struct quad* bottomleft; //a pointer to the bottomleft sub-quadrant
  struct quad* bottomright; //a pointer to the bottomright sub-quadrant
};
typedef struct quad galquad_t;

const int windowWidth = 800;
float circleRadius=0.002, circleColor=0;
float maxBrightness=0;
float maxSize=0;
float minSize=1e10;

//x and y are the coordinates of the top left corner of the quadrant that we are adding particle particle to
//If *tree == NULL then make a new quadrant at this point and put particle into it.
//If *tree != NULL the x and y coords are the same as *tree->xstart and *tree->ystart
void add_to_quadtree(galquad_t** tree, particle_t* particle, double x, double y, double size)
{
  double half = size/2;
  //printf("add_to_quadtree(tree %p, particle %p, %f, %f, %f)\n", *tree, particle, x, y, size);
  if((*tree)==NULL)//new quad to be created here
  {
    //printf("Tree = NULL, creating new quad\n");
    (*tree) = (galquad_t*)malloc(sizeof(galquad_t));
    (*tree)->particle = particle;
    (*tree)->xstart = x;
    (*tree)->ystart = y;
    (*tree)->size = size;
    (*tree)->x = particle->x; //unless another particle is added then the x and y positions are the same as those for the particle
    (*tree)->y = particle->y;
    (*tree)->m = particle->mass; //and the mass is also the same (this quadrant acts as a single particle if that is all that is in it
    (*tree)->topleft = NULL;
    (*tree)->bottomleft = NULL;
    (*tree)->topright = NULL;
    (*tree)->bottomright = NULL;
    return;
  }
  if((*tree)->particle != NULL) //this quad contains one particle - move it to a subquad
  {
    //printf("Tree contains a particle, moving it to a subtree\n");
    if((*tree)->particle->x < x + half && (*tree)->particle->y < y + half)//top left
      add_to_quadtree(&((*tree)->topleft), (*tree)->particle, x, y, half);
    else if((*tree)->particle->x < x + half)//bottom left
      add_to_quadtree(&((*tree)->bottomleft), (*tree)->particle, x, y+half, half);
    else if((*tree)->particle->y < y + half)//top right
      add_to_quadtree(&((*tree)->topright), (*tree)->particle, x+half, y, half);
    else
      add_to_quadtree(&((*tree)->bottomright), (*tree)->particle, x+half, y+half, half);
    (*tree)->particle = NULL;
    (*tree)->x = 0; //for a quad with subquads, these will be added later once the tree is complete.
    (*tree)->y = 0;
    (*tree)->m = 0;
  }
  //now we have dealt with any existing particle in this quad we need to add the new particle
  //printf("Now adding the new particle to one of the subtrees\n");
  if(particle->x < x + half && particle->y < y + half) //top left
    add_to_quadtree(&((*tree)->topleft), particle, x, y, half);
  else if(particle->x < x + half) //bottom left
    add_to_quadtree(&((*tree)->bottomleft), particle, x, y+half, half);
  else if(particle->y < y + half) //top right
    add_to_quadtree(&((*tree)->topright), particle, x+half, y, half);
  else //bottom right
    add_to_quadtree(&((*tree)->bottomright), particle, x+half, y+half, half);
}

void clear_quadtree(galquad_t* tree)
{
  if(tree->topleft != NULL)
  {
    clear_quadtree(tree->topleft);
    free(tree->topleft);
    tree->topleft = NULL;
  }
  if(tree->bottomleft != NULL)
  {
    clear_quadtree(tree->bottomleft);
    free(tree->bottomleft);
    tree->bottomleft = NULL;
  }
  if(tree->topright != NULL)
  {
    clear_quadtree(tree->topright);
    free(tree->topright);
    tree->topright = NULL;
  }
  if(tree->bottomright != NULL)
  {
    clear_quadtree(tree->bottomright);
    free(tree->bottomright);
    tree->bottomright = NULL;
  }
}

void initialize_graphics(char* prog, particle_t* particles, double* br, int N)
{
  int i;
  InitializeGraphics(prog, windowWidth, windowWidth);
  for(i=0; i<N; ++i) {
    if(br[i] > maxBrightness) maxBrightness = br[i];
    if(particles[i].mass > maxSize) maxSize = particles[i].mass;
    if(particles[i].mass < minSize) minSize = particles[i].mass;
  }
  //set the top of the colour scale to more than the highest value we will actually use so
  //that we don't get any black or very dark grey spheres
  SetCAxes(0,maxBrightness*1.5);
}
void update_graphics(particle_t* particles, double* br, int N)
{
  ClearScreen();
  int i;
  for(i=0; i<N; ++i)
    DrawCircle(particles[i].x, particles[i].y, 1, 1, circleRadius*(1+(particles[i].mass-minSize)/(maxSize-minSize)), br[i]);
  Refresh();
  usleep(20000);
}
void close_graphics()
{
  FlushDisplay();
  CloseDisplay();
}

void print_particles(particle_t* particles, double* br, int N)
{
  int i;
  for(i=0; i<N; ++i)
  {
    printf("Particle %d: position (%2.4f, %2.4f); mass %2.4f; velocity (%2.4f, %2.4f); brightness %2.4f\n", i, particles[i].x, particles[i].y, particles[i].mass, particles[i].vx, particles[i].vy, br[i]);
  }
}

//Calculate theta for a given particle and tree
double find_theta(particle_t *p, galquad_t* t)
{
  //double half = t->size / 2;
  double dx, dy;
  dx = t->xstart + t->size / 2 - p->x; //x distance between centre of quad and particle
  dy = t->ystart + t->size / 2 - p->y; //y distance between centre of quad and particle
  return t->size / sqrt(pow(dx, 2)+pow(dy, 2));
}


//Calculate the forces acting on the given particle resulting from the given (sub)tree of galaxy quadrants
//and add them to whatever is passed in fx and fy
void force(particle_t* p, galquad_t* t, double* fx, double* fy, double theta_max, double G)
{
  double epsilon = 1e-3;
  double factor;
  double dx, dy;
  //double addx=0, addy=0;

  //first deal with the case where the tree we are passed only contains this particle.
  if(t->particle == p)
  {
    //printf("Tree contains only the current particle. No force applies\n");
    return;
  }
  

  //next case, tree only contains one particle (but not this one). Tree is treated as a point mass.
  if(t->particle != NULL)
  {
    dx = (p->x - t->x);
    dy = (p->y - t->y);
    factor = sqrt(pow(dx,2) + pow(dy,2)) + epsilon;
    //factor = t->m/(pow(factor, 3));
    *fx += (p->x - t->x)*t->m/(pow(factor, 3))*G;
    *fy += (p->y - t->y)*t->m/(pow(factor, 3))*G;
    return;
    //printf("Adding force of a single particle to sum: dx=%f, dy=%f, factor=%f, m=%f, fx=%f, fy=%f\n", dx, dy, factor, t->m, addx, addy);
  }

  
  //last case, tree contains one to four sub-trees. For each subtree we need to check whether
  //theta > theta_max (in which case treat as a point mass), or not (in which case recurse).
  //printf("Tree has subtrees. Take them in turn\n");
  if(t->particle != p && t->topleft != NULL)
  {
    //printf("Adding force of topleft tree to current particle\n");
    if(find_theta(p, t->topleft) > theta_max)
    {
      force(p, t->topleft, fx, fy, theta_max, G);
    }
    else //treat tree as a point mass
    {
      dx = (p->x - t->topleft->x);
      dy = (p->y - t->topleft->y);
      factor = sqrt(pow(dx,2) + pow(dy,2)) + epsilon;
      //factor = t->topleft->m/(pow(factor, 3));
      *fx += (p->x - t->x)*t->m/(pow(factor, 3))*G;
      *fy += (p->y - t->y)*t->m/(pow(factor, 3))*G;    
    }
  }

  if(t->topright != NULL)
  {
    //printf("Adding force of topright tree to current particle\n");
    if(find_theta(p, t->topright) > theta_max)
    {
      force(p, t->topright, fx, fy, theta_max, G);
    }
    else //treat tree as a point mass
    {
      dx = (p->x - t->topright->x);
      dy = (p->y - t->topright->y);
      factor = sqrt(pow(dx,2) + pow(dy,2)) + epsilon;
      // factor = t->topright->m/(pow(factor, 3));
      *fx += (p->x - t->x)*t->m/(pow(factor, 3))*G;
      *fy += (p->y - t->y)*t->m/(pow(factor, 3))*G;
    
    }
  }
  if(t->bottomleft != NULL)
  {
    //printf("Adding force of bottomleft tree to current particle\n");
    if(find_theta(p, t->bottomleft) > theta_max)
    {
      force(p, t->bottomleft, fx, fy, theta_max, G);
    }
    else //treat tree as a point mass
    {
      dx = (p->x - t->bottomleft->x);
      dy = (p->y - t->bottomleft->y);
      factor = sqrt(pow(dx,2) + pow(dy,2)) + epsilon;
      //factor = t->bottomleft->m/(pow(factor, 3));
      *fx += (p->x - t->x)*t->m/(pow(factor, 3))*G;
      *fy += (p->y - t->y)*t->m/(pow(factor, 3))*G;   
    }
  }
  if(t->bottomright != NULL)
  {
    //printf("Adding force of bottomright tree to current particle\n");
    if(find_theta(p, t->bottomright) > theta_max)
    {
      force(p, t->bottomright, fx, fy, theta_max, G);
    }
    else //treat tree as a point mass
    {
      dx = (p->x - t->bottomright->x);
      dy = (p->y - t->bottomright->y);
      factor = sqrt(pow(dx,2) + pow(dy,2)) + epsilon;
      // factor = t->bottomright->m/(pow(factor, 3));
      *fx += (p->x - t->x)*t->m/(pow(factor, 3))*G;
      *fy += (p->y - t->y)*t->m/(pow(factor, 3))*G;    
    }
  }
  //*fx += addx;
  //*fy += addy;
}

//Perform one step of the simulation for <N> particles <particles>, depicted by <tree> quad tree, with timestep <dt>
//and gravitational constant -<G> (n.b. negative). Where a particle is at a relative distance of theta_max or less
//from a particular quadrant, treat particles in that quadrant as a point mass.
void sim_step(particle_t* particles, galquad_t* tree, int N, const double dt, const double G, const double theta_max)
{
  int i;
  double* fx = (double*)calloc(N,sizeof(double));
  double* fy = (double*)calloc(N,sizeof(double));
  double* ax = (double*)malloc(N*sizeof(double));
  double* ay = (double*)malloc(N*sizeof(double));

  //find forces acting on each particle (ignoring gravity)
  //for(i=0; i<N; ++i)
  //{
    //printf("Calculate forces for particle %d:\n", i);
  //}

  //once we have all the forces, update the particles to their new positions and velocities.
  for(i=0; i<N; ++i)//try unrolling this loop and/or use vector ops
  {
    force(&particles[i], tree, &fx[i], &fy[i], theta_max, G);
    //printf("applying force %f, %f to particle %d\n", fx[i], fy[i], i);
    ax[i] = fx[i]*dt;
    ay[i] = fy[i]*dt;
    particles[i].vx += ax[i];
    particles[i].vy += ay[i];
    particles[i].x += dt*particles[i].vx;
    particles[i].y += dt*particles[i].vy;
  }
  free(fx);
  free(fy);
  free(ax);
  free(ay);
}

//recursive function to populate all the quadtree structures with the total mass and centre of mass
//for that quadrant
void add_masses(galquad_t* tree)
{
  if(tree->particle != NULL)
  { //already dealt with by the tree constructor
    //tree->m = tree->particle->mass;
    //tree->x = tree->particle->x;
    //tree->y = tree->particle->y;
    //printf("Add masses: (sub)tree %p has mass %f centred at (%f, %f) [single particle]\n", tree, tree->m, tree->x, tree->y);
    return;
  }
  tree->m = 0;
  tree->x = 0;
  tree->y = 0;
  if(tree->topleft != NULL)
  {
    add_masses(tree->topleft);
    tree->m = tree->topleft->m;
    tree->x = tree->topleft->m * tree->topleft->x;
    tree->y = tree->topleft->m * tree->topleft->y;
  }
  if(tree->topright != NULL)
  {
    add_masses(tree->topright);
    tree->m += tree->topright->m;
    tree->x += tree->topright->m * tree->topright->x;
    tree->y += tree->topright->m * tree->topright->y;
  }
  if(tree->bottomleft != NULL)
  {
    add_masses(tree->bottomleft);
    tree->m += tree->bottomleft->m;
    tree->x += tree->bottomleft->m * tree->bottomleft->x;
    tree->y += tree->bottomleft->m * tree->bottomleft->y;
  }
  if(tree->bottomright != NULL)
  {
    add_masses(tree->bottomright);
    tree->m += tree->bottomright->m;
    tree->x += tree->bottomright->m * tree->bottomright->x;
    tree->y += tree->bottomright->m * tree->bottomright->y;
  }
  //we have added the masses and weighted distances for up to 4 sub-quads, now rescale
  tree->x = tree->x / (tree->m);
  tree->y = tree->y / (tree->m);
  //printf("Add masses: (sub)tree %p has mass %f centred at (%f, %f)\n", tree, tree->m, tree->x, tree->y);
}

//Perform galaxy simulation for <N> <particles> with brightnesses <br> (only used for <graphics>), over <steps> 
//iterations of timestep <td> with Barnes-Hut threshold <theta_max>
void sim_gal(particle_t* particles, double* br, int N, int steps, double dt, double theta_max, int graphics)
{
  //printf("Particles before calculation:\n");
  //print_particles(particles, br, N);
  int step, i;
  double G = ((double)(-100))/N;
  galquad_t* tree = NULL;
  for(step=0; step<steps; ++step)
  {
    //printf("Beginning step %d:\n", step);
    //in each step we create the quad tree from scratch
    for(i=0; i<N; ++i)
    {
      //printf("Add particle %d (%f, %f, %p) to tree:\n", i, particles[i].x, particles[i].y, &particles[i]);
      add_to_quadtree(&tree, &particles[i], 0, 0, 1);
    }
    //printf("  Created tree.\n");

    //once the tree is complete add mass and centre of mass to each quad
    add_masses(tree);
    //printf("  Added masses to tree.\n");

    //then perform the simulation
    sim_step(particles, tree, N, dt, G, theta_max);
    //printf("  Completed sim step %d.\n", step);

    //update the graphics
    if(graphics && !(step%5))
    {
      update_graphics(particles, br, N);
      //printf("  Updated graphics\n");
    }

    //clear out the quad tree ready for the next step
    clear_quadtree(tree);
    free(tree);
    tree = NULL;
    //printf("  Completed tree clean.\n");
  }
  //printf("Particles after calculation:\n");
  //print_particles(particles, br, N);
}

//write the particles and brightnesses to a file.
int saveoutput(particle_t* particles, double* br, int N)
{
  char* outfile = "result.gal";
  FILE* outfileptr = fopen(outfile, "wb");
  int i;
  const int size = sizeof(particle_t);
  if(!outfileptr) {printf("Failed to open file %s\n", outfile); return -1;};
  for(i=0; i<N; ++i)
  {
    fwrite(&particles[i], size, 1, outfileptr);
    fwrite(&br[i], sizeof(double), 1, outfileptr);
  }
  fclose(outfileptr);
  return 0;
}

int main(int argc, char** argv)
{
  if(argc != 7)
  {
    printf("Usage: %s <N> <filename> <nsteps> <delta_t> <theta_max> graphics\n", argv[0]);
    return -1;
  }

  // For time measurements
  struct timespec t0, t1;
  clock_gettime(CLOCK_REALTIME, &t0);

  int N=atoi(argv[1]);
  particle_t *particles = (particle_t*)malloc(N*sizeof(particle_t));
  double *br = (double*)malloc(N*sizeof(double));
  //Read the (first N*6 lines of the) input file
  FILE* infile = fopen(argv[2], "rb");
  if(!infile) {printf("Failed to open file %s\n", argv[2]); return -1;}

  int i;
  for(i=0; i<N; ++i)
  {
    if (fread(&(particles[i].x), sizeof(double), 1, infile) != 1)
      printf("Could not read x-position\n");
    if (fread(&(particles[i].y), sizeof(double), 1, infile) != 1)
      printf("Could not read y-position\n");
    if (fread(&(particles[i].mass), sizeof(double), 1, infile) != 1)
      printf("Could not read mass\n");
    if (fread(&(particles[i].vx), sizeof(double), 1, infile) != 1)
      printf("Could not read x-velocity\n");
    if (fread(&(particles[i].vy), sizeof(double), 1, infile) != 1) 
      printf("Could not read y-velocity\n");
    //read brightness into a separate array to keep particle memory usage smaller
    if(fread(&(br[i]), sizeof(double), 1, infile) != 1)
      printf("Could not read brightness\n");
  }
  fclose(infile);
  //print_particles(particles, br, N);
  int graphics;
  if((graphics=atoi(argv[6])))
    initialize_graphics(argv[0], particles, br, N);
  sim_gal(particles, br, N, atoi(argv[3]), atof(argv[4]), atof(argv[5]), graphics);
  if(graphics)
    close_graphics();
  int ret = saveoutput(particles, br, N);
  free(particles);
  free(br);

  clock_gettime(CLOCK_REALTIME, &t1);
  long elapsed_time_nsec = (t1.tv_sec-t0.tv_sec)*1e9 + t1.tv_nsec-t0.tv_nsec;
  double elapsed_time_sec = (t1.tv_sec-t0.tv_sec) + (t1.tv_nsec-t0.tv_nsec)/1e9;
  printf("%ld nano sec, %lf sec\n", elapsed_time_nsec, elapsed_time_sec);
  return ret;
}
