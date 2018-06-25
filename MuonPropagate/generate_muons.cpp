#include "PROPOSAL/Propagator.h"
#include <iostream>
#include <math.h>

vector<double> linspace(double Emin,double Emax,int div){
    vector<double> linpoints;
    double step_lin = (Emax - Emin)/double(div);

    double EE = Emin;
    while (EE <= Emax+0.001){
        linpoints.push_back(EE);
        EE = EE + step_lin;
    }

    return linpoints;
}

vector<double> logspace(double Emin,double Emax,int div){
    vector<double> logpoints;
    double Emin_log,Emax_log;
    if (Emin < 1.0e-5 ) {
        Emin_log = 0.0;
    } else {
        Emin_log = log(Emin);
    }
    Emax_log = log(Emax);

    double step_log = (Emax_log - Emin_log)/double(div);

    double EE = Emin_log;
    while (EE <= Emax_log+0.001){
        logpoints.push_back(exp(EE));
        EE = EE + step_log;
    }

    return logpoints;
}

using namespace std;

int main(){
  const double GeV = 1.0e3;
  const double km = 1.0e5;
  const double meter = 1.0e2;

  EnergyCutSettings* ecut = new EnergyCutSettings(10,1e-3);
  Medium* ICE = new Medium("ice",1.);
  Propagator* ice_prop = new Propagator(ICE,ecut,"mu","resources/tables");//,false,true,true,true,1,12,1,1,1,1,false,2);

  unsigned int number_of_particles = 1e4;
  unsigned int OnePercent = number_of_particles/100;

  std::vector<double> initial_energy_array = logspace(1.0e2*GeV,1.0e7*GeV,100); // in MeV
  std::vector<double> distance_array = linspace(1.5*km,150.0*km,100); // in cm

    for( double initial_energy : initial_energy_array ){
      for( double distance : distance_array ){
        for( unsigned int i = 0; i < number_of_particles; i++ ){
          double ice_distance = distance;
          // then propagate in ice
          ice_prop->GetParticle()->SetEnergy(initial_energy);
          ice_prop->GetParticle()->SetPhi(0);
          ice_prop->GetParticle()->SetTheta(0);
          ice_prop->GetParticle()->SetX(0);
          ice_prop->GetParticle()->SetY(0);
          ice_prop->GetParticle()->SetZ(0);
          ice_prop->GetParticle()->SetT(0);
          ice_prop->GetParticle()->SetPropagatedDistance(0);
          ice_prop->GetParticle()->SetParticleId(0);

          double final_energy = ice_prop->GetParticle()->GetEnergy();

          std::cout << initial_energy << " " << distance << " " << final_energy << std::endl;
        }
      }
    }

  return 0;
}
