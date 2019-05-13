#include <stdio.h>
#include "i3module.h"

/**
 * Follow these steps:
 *    1. call initJre() to create Java runtime environment
 *       call setStderr(filename) to redirect stderr (optional)
 *    2. call initMMC(options, flag) to initialize the mmc module:
 *       flag=1: to run tfa/Amanda
 *       flag=2: to run gen/AtmFlux as generator
 *       flag=3: to run gen/AtmFlux as propagator
 *    3. call propagate(...) for flag=[1/3] or createNext() for flag=2
 *    4. call getNum() to get the number of particles in the returned event
 *    5. call startRead(i) before reading particle #i
 *    6. call endRead() when done reading event particles
 *    7. call endProp() when you are done with the processed event
 *    8. at the end of your program, or when Jre is no longer necessary
 *       call stopJre() to stop the Jre
 */

char *name;
double x, y, z, theta, phi, e, t;

main(){
  int events=10; // number of events to generate
  int n, i;

  char *options[]={
    "-seed=1 -user -sdec -raw -rw=0",
    "-seed=1 -user -Emu=1.e7 -Enm=1.e6 -Ene=1.e6 -sdec -raw",
    "-seed=1 -user -sdec -raw -prop"
  };

  double r, z1=0, z2=0, h1=0, h2=0, nx=0, ny=0, nz=0, e1=0, e2=0, ec=0;

  initJre();

  // setStderr("mmc.log");
  printf("\nRunning %s\n", "tfa/Amanda (1)");
  initMMC(options[0], 1);

  for(n=0; n<events; n++){
    printf("Event %i\n", n);
    setStart(1);
    name="mu-";
    x=1;    // [m]
    y=2;    // [m]
    z=1000; // [m]
    theta=0;
    phi=0;
    e=2000; // [GeV]
    t=0;    // [s]
    propagate(name, x*1.e2, y*1.e2, z*1.e2, 180-theta, phi<180?phi+180:phi-180, e*1.e3, t);

    resultsOut();

    r=getDouble("r", 1);
    e=getDouble("e", 1);
    z1=getDouble("z1", 2);
    z2=getDouble("z2", 2);
    h1=getDouble("h1", 2);
    h2=getDouble("h2", 2);
    nx=getDouble("nx", 2);
    ny=getDouble("ny", 2);
    nz=getDouble("nz", 2);
    e1=getDouble("e1", 2);
    e2=getDouble("e2", 2);
    ec=getDouble("ec", 2);

    r*=1.e-2; e*=1.e-3;
    printf("%g %g %g %g %g %g %g %g %g %g %g %g\n", r, e, z1, z2, h1, h2, nx, ny, nz, e1, e2, ec);
    setStart(0);
    printf("EVENT WEIGHT = %g at %g\n", getDouble("fw", 2), getDouble("hw", 2));
    endProp();
  }
  deleteMMC();

  printf("\nRunning %s\n", "gen/AtmFlux (2)");
  initMMC(options[1], 2);
  for(n=0; n<events; n++){
    createNext();
    resultsOut();
    endProp();
  }
  deleteMMC();

  printf("\nRunning %s\n", "gen/AtmFlux (3)");
  initMMC(options[2], 3);
  for(n=0; n<events; n++){
    printf("Event %i\n", n);

    name="nu_tau";
    x=1;    // [m]
    y=2;    // [m]
    z=-5e6; // [m]
    theta=180;
    phi=0;
    e=1.e8; // [GeV]
    t=0;    // [s]
    propagate(name, x*1.e2, y*1.e2, z*1.e2, 180-theta, phi<180?phi+180:phi-180, e*1.e3, t);

    resultsOut();
    endProp();
  }
  deleteMMC();

  stopJre();
}



resultsOut(){
  int i;

  // get the results

  for(i=0; i<getNum(); i++){
    startRead(i);

    name=getName(1);
    x=getDouble("x", 0);
    y=getDouble("y", 0);
    z=getDouble("z", 0);
    theta=getDouble("theta", 0);
    phi=getDouble("phi", 0);
    t=getDouble("t", 0);
    e=getDouble("e", 0);

    x*=1.e-2; y*=1.e-2; z*=1.e-2; theta=180-theta; phi=phi<180?phi+180:phi-180; e*=1.e-3;
    printf("%s %g %g %g %g %g %g %g\n", name, x, y, z, theta, phi, t, e);

    endRead();
  }
}
