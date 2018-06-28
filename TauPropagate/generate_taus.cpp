#include "PROPOSAL/Propagator.h"
#include <iostream>
#include <math.h>


using namespace std;

std::vector<double> linspace(double Emin,double Emax,int div){
    std::vector<double> linpoints;
    double step_lin = (Emax - Emin)/double(div);

    double EE = Emin;
    while (EE <= Emax+0.001){
        linpoints.push_back(EE);
        EE = EE + step_lin;
    }

    return linpoints;
}

std::vector<double> logspace(double Emin,double Emax,int div){
    std::vector<double> logpoints;
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

int main(int argc, char* argv[]){
  //const double PeV = 1.0e9;
  const double GeV = 1.0e3;
  const double km = 1.0e5;
  const double meter = 1.0e2;

  PROPOSAL::Propagator* ice_prop = new PROPOSAL::Propagator(PROPOSAL::TauMinusDef::Get(),"resources/config_ice.json"); //,false,true,true,true,1,12,1,1,1,1,false,2);//,false,true,true,true,1,12,1,1,1,1,false,2);
  //ice_prop->set_seed(seed);
  std::cout << "Made the Propagator" << std::endl;
  unsigned int number_of_particles = 1e4;
  unsigned int OnePercent = number_of_particles/100;
  double initial_energy;
  initial_energy = atof(argv[1]);
  //std::vector<double> initial_energy_array = logspace(1e10*GeV,1e13*GeV,50); // in MeV
  std::vector<double> distance_array = linspace(1.5*km,150.0*km,50); // in cm
  PROPOSAL::Vector3D position(0., 0., 0.);
  PROPOSAL::Vector3D direction(1., 1., 1.);

  std::cout << initial_energy << std::endl;
//    for( double initial_energy : initial_energy_array ){
    for( double distance : distance_array ){
      for( unsigned int i = 0; i < number_of_particles; i++ ){
        double ice_distance = distance;
          // then propagate in ice

        ice_prop->GetParticle().SetEnergy(initial_energy);
          //ice_prop->GetParticle().SetPhi(0);
        //ice_prop->GetParticle()->SetTheta(0);
        //ice_prop->GetParticle()->SetX(0);
        //ice_prop->GetParticle()->SetY(0);
        //ice_prop->GetParticle()->SetZ(0);
        //ice_prop->GetParticle()->SetT(0);
        ice_prop->GetParticle().SetPosition(position);
        ice_prop->GetParticle().SetDirection(direction);

        ice_prop->GetParticle().SetPropagatedDistance(0);
        ice_prop->GetParticle().SetParticleId(0);
        //std::cout << "particle set" << std::endl;

        ice_prop->Propagate(ice_distance);


        //std::cout<< "Propaaaa" << std::endl;
        double final_energy = ice_prop->GetParticle().GetEnergy();

        double final_distance = ice_prop->GetParticle().GetPropagatedDistance();

        std::cout << final_distance << " " << final_energy << std::endl;
        }
      }
  return 0;
}
