#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Output.h"
#include "process_muons.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>

#include <boost/histogram.hpp>
#include <boost/histogram/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/histogram/histogram_ostream_operators.hpp>

vector<double> linspace(double Emin, double Emax, unsigned int n_edges){
    vector<double> linpoints;
    double step_lin = (Emax - Emin)/double(n_edges-1);

    double EE = Emin;
    while (EE <= Emax+0.0001){
        linpoints.push_back(EE);
        EE = EE + step_lin;
    }

    return linpoints;
}

std::vector<double> logspace(double Emin_log, double Emax_log, unsigned int n_edges, double index=10.) {
    std::vector<double> logpoints;
    double step_log = (Emax_log - Emin_log)/double(n_edges-1);
    double EE = Emin_log;
    while (EE <= Emax_log+0.0001){
        logpoints.push_back(std::pow(index,EE));
        EE = EE + step_log;
    }
    return logpoints;
}

double overburden(double cos_theta, double depth=9.882075327949499, double r_E=32387.868420987306, double m=5.0677309374099995e-3) {
    double d = depth;
    double r = r_E;
    double z = r-d;
    return (std::sqrt(std::pow(z,2)*std::pow(cos_theta,2)+d*(2*r-d))-z*cos_theta)/m;
}

int main(int argc, char * argv[])
{
    const double GeV = 1.0e3;
    const double km = 1.0e5;
    const double meter = 1.0e2;

    unsigned int energy_index;
    std::stringstream ss;
    std::string s(argv[1]);
    ss.str(s);
    ss >> energy_index;
    //energy *= GeV;

    int seed;
    ss.clear();
    s = std::string(argv[2]);
    ss.str(s);
    ss >> seed;

    int n_muons;
    ss.clear();
    s = std::string(argv[3]);
    ss.str(s);
    ss >> n_muons;

    EnergyCutSettings* ecut = new EnergyCutSettings(10,1e-4);
    Medium* ICE = new Medium("ice",1.);
    Propagator* ice_prop = new Propagator(ICE,ecut,"mu","resources/tables");
    ice_prop->set_seed(seed);

    ss.str("");
    ss.clear();
    ss << "./" << energy_index;
    std::string dir_name = ss.str();
    mkdir(dir_name.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH);
    ss <<"/" << seed << ".txt";

    std::string file_name = ss.str();

    std::cout << file_name << std::endl;
    
    std::ofstream f(file_name);
    //boost::archive::text_oarchive out_archive(f);
    std::vector<double> e_edges = logspace(0, 11, 88*3+1);
    std::vector<double> e_grid = e_edges;

    std::adjacent_difference(e_grid.begin(), e_grid.end(), e_grid.begin(), [](double a, double b)->double{return (a+b)/2.0;});
    std::reverse(e_grid.begin(), e_grid.end());
    std::reverse(e_edges.begin(), e_edges.end());
    e_grid.pop_back();
    std::cout << "e_grid: " << e_grid.size() << std::endl;

    std::vector<double> L_grid = linspace(-0.05, 1, 1050+1);
    std::transform(L_grid.begin(), L_grid.end(), L_grid.begin(), [](double a)->double{return overburden(a);});
    std::reverse(L_grid.begin(), L_grid.end());
    
    std::cout << "L_grid: " << L_grid.size() << std::endl;
    
    muons::MuonHistogrammer<double> muon_hist(e_edges, e_grid, L_grid, energy_index);
    
    int id = 0;

    std::function<void(Particle*,int,std::pair<double,std::string>,double)> callback = [&] (Particle *particle, int secondary_id, pair<double,string> energy_loss, double distance) {
        double d = std::sqrt(std::pow(particle->GetX(), 2) + std::pow(particle->GetY(), 2) + std::pow(particle->GetZ(), 2))/meter;
        double dE = energy_loss.first/GeV;
        double E = particle->GetEnergy()/GeV;
        if(secondary_id != id) {
            id = secondary_id;
            muon_hist.reset(energy_index);
        }
        muon_hist.process_checkpoint(muons::MuonCheckpoint<double>(d, dE, E));
    };

    Output::getInstance().SetSecondaryCallback(callback);

    for(int i=0; i<n_muons; ++i)
    {
        ice_prop->GetParticle()->SetEnergy(e_grid[energy_index]*GeV);
        ice_prop->GetParticle()->SetPhi(0);
        ice_prop->GetParticle()->SetTheta(0);
        ice_prop->GetParticle()->SetX(0);
        ice_prop->GetParticle()->SetY(0);
        ice_prop->GetParticle()->SetZ(0);
        ice_prop->GetParticle()->SetT(0);
        ice_prop->GetParticle()->SetPropagatedDistance(0);
        ice_prop->GetParticle()->SetParticleId(i);

        ice_prop->Propagate(1e9*meter);

        double final_energy = ice_prop->GetParticle()->GetEnergy();

        double final_distance = ice_prop->GetParticle()->GetPropagatedDistance();
        Output::getInstance().ClearSecondaryVector();
        if(i%100 == 0) {
            std::cout << i << std::endl;
        }
    }

    namespace bh = boost::histogram;
    using namespace bh::literals;

    //out_archive << muon_hist.hist;
    //out_archive << muon_hist.zero_hist;
    f << muon_hist << std::endl;
}

