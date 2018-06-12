#include <vector>
#include <algorithm>
#include <iostream>
#include <boost/histogram.hpp>
#include <boost/histogram/serialization.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/fusion/include/mpl.hpp>

namespace bh = boost::histogram;

namespace muons {

template<typename T>
struct MuonCheckpoint {
    T distance;
    T deltaE;
    T E;
    MuonCheckpoint(T d, T de, T e) : distance(d), deltaE(de), E(e) {};
};

const static double min_energy = .1056583745;

template<typename T>
class MuonHistogrammer{
    std::vector<T> const & energy_grid;
    std::vector<T> const & energy_edges;
    std::vector<T> const & distance_grid;
    std::vector<unsigned int> check_energy_indices;
    std::vector<unsigned int> check_distance_indices;
    std::vector<T> check_distance_offsets;

    public:

    friend std::ostream &operator<<(std::ostream & os, MuonHistogrammer<T> const & obj) {
        std::vector<unsigned int> dims;
        unsigned int d0 = obj.hist.axis(boost::mpl::int_<0>()).bins();
        unsigned int d1 = obj.hist.axis(boost::mpl::int_<1>()).bins();
        unsigned int d2 = obj.hist.axis(boost::mpl::int_<2>()).bins();
        //unsigned int d0 = dims[0];
        //unsigned int d1 = dims[1];
        //unsigned int d2 = dims[2];
        os << "(";
        for(unsigned int i=0; i<d0; ++i) {
            os << "(";
            for(unsigned int j=0; j<d1; ++j) {
                os << "(";
                for(unsigned int k=0; k<d2; ++k) {
                    os << obj.hist.value(i,j,k);
                    os << ",";
                }
                if(d2 > 0)
                    os << "\b";
                os << ")";
                os << ",";
            }
            if(d1 > 0)
                os << "\b";
            os << ")";
            os << ",";
        }
        if(d0 > 0)
            os << "\b";
        os << ")";
        os << std::endl;
        os << "(";
        for(unsigned int i=0; i<d0; ++i) {
            os << "(";
            for(unsigned int j=0; j<d1; ++j) {
                os << obj.zero_hist.value(i, j);
                os << ",";
            }
            if(d1 > 0)
                os << "\b";
            os << ")";
            os << ",";
        }
        if(d0 > 0)
            os << "\b";
        os << ")";
        return os;
    };

    boost::histogram::histogram<std::integral_constant<int, 0>, boost::mpl::vector<boost::histogram::axis::integer<int>, boost::histogram::axis::integer<int>, boost::histogram::axis::variable<T> >, boost::histogram::adaptive_storage<std::allocator> > hist;
    boost::histogram::histogram<std::integral_constant<int, 0>, boost::mpl::vector<boost::histogram::axis::integer<int>, boost::histogram::axis::integer<int> >, boost::histogram::adaptive_storage<std::allocator> > zero_hist;

    unsigned int counter = 0;
    unsigned int zero_counter = 0;

    private:
    unsigned int max_E_index;
    MuonCheckpoint<T> previous_checkpoint;
    MuonCheckpoint<T> next_checkpoint;
    T previous_energy;
    
    public:

    MuonHistogrammer(std::vector<T> const & e_edges, std::vector<T> const & e_grid, std::vector<T> const & d_grid, unsigned int E_index) : 
        energy_edges(e_edges),
        energy_grid(e_grid),
        distance_grid(d_grid),
        previous_checkpoint(0, 0, energy_grid[E_index]),
        next_checkpoint(0, 0, energy_grid[E_index]) {
        hist = bh::make_static_histogram(bh::axis::integer<>(0, e_grid.size()-1), bh::axis::integer<>(0, d_grid.size()-1), bh::axis::variable<>(e_edges.begin(), e_edges.end()));
        zero_hist = bh::make_static_histogram(bh::axis::integer<>(0, e_grid.size()-1), bh::axis::integer<>(0, d_grid.size()-1));
        reset(E_index);
    }

    void reset(unsigned int E_index) {
        max_E_index = E_index;
        previous_checkpoint = MuonCheckpoint<T>(0, 0, energy_grid[E_index]);
        next_checkpoint = previous_checkpoint;
        previous_energy = previous_checkpoint.E - previous_checkpoint.deltaE;
        add_check(E_index, 0, 0);
    }

    void add_check(unsigned int E_index, unsigned int L_index, T d_offset) {
        check_energy_indices.push_back(E_index);
        check_distance_indices.push_back(L_index);
        check_distance_offsets.push_back(d_offset);
        max_E_index = std::max(max_E_index, E_index);
    }

    void perform_check(unsigned int check_index) {
        bool done = false;
        unsigned int E_index = check_energy_indices[check_index];
        unsigned int L_index = check_distance_indices[check_index];
        T d_offset = check_distance_offsets[check_index];

        MuonCheckpoint<T> & pc = previous_checkpoint;
        MuonCheckpoint<T> & nc = next_checkpoint;

        T L;
        while((!done) && (L_index < distance_grid.size())) {
            L = distance_grid[L_index];
            if(nc.distance > L+d_offset) {
                T Ef = previous_energy - (previous_energy - nc.E)/(nc.distance - pc.distance)*(L + d_offset - pc.distance);
                add_histogram_entry(E_index, L_index, Ef);
                L_index += 1;
            }
            else {
                done = true;
            }
        }
        check_distance_indices[check_index] = L_index;
    }

    void check_zero(unsigned int check_index) {
        // Look at the L index we are currently on
        // Add zero (or min_energy) to the histogram for the current L and any larger L
        unsigned int E_index = check_energy_indices[check_index];
        unsigned int L_index = check_distance_indices[check_index];
        T d_offset = check_distance_offsets[check_index];

        MuonCheckpoint<T> & nc = next_checkpoint;

        T L;
        while(L_index < distance_grid.size()) {
            L = distance_grid[L_index];
            if(L+d_offset >= nc.distance) {
                add_zero_histogram_entry(E_index, L_index);
                L_index += 1;
            }
        }
        check_distance_indices[check_index] = L_index;
    }

    void update_checks(unsigned int E_index) {
        bool done = false;
        T next_Ei;
        
        MuonCheckpoint<T> & pc = previous_checkpoint;
        MuonCheckpoint<T> & nc = next_checkpoint;

        while((!done) && (E_index < (energy_grid.size() - 1))) {
            next_Ei = energy_grid[E_index+1];
            if(nc.E < next_Ei && previous_energy > next_Ei) {
                T Ef = pc.distance + (previous_energy - next_Ei)*(nc.distance - pc.distance)/(previous_energy - nc.E);
                add_check(E_index+1, 0, Ef);
            }
            if(nc.E-nc.deltaE < next_Ei) {
                E_index += 1;
            }
            else {
                done = true;
            }
        }
    }

    void process_checkpoint(MuonCheckpoint<T> const & checkpoint) {
        previous_checkpoint = next_checkpoint;
        next_checkpoint = checkpoint;
        previous_energy = previous_checkpoint.E - previous_checkpoint.deltaE;
        update_checks(max_E_index);

        unsigned int i=0;
        while(i < check_energy_indices.size()) {
            perform_check(i);
            i += 1;
        }
        if((next_checkpoint.E - min_energy) / min_energy <= 0.001) {
            i = 0;
            while(i < check_energy_indices.size()) {
                check_zero(i);
                i += 1;
            }
        }
    }

    void add_histogram_entry(unsigned int E_index, unsigned int L_index, T Ef) {
        hist.fill(E_index, L_index, Ef);
        ++counter;
    }
    
    void add_zero_histogram_entry(unsigned int E_index, unsigned int L_index) {
        zero_hist.fill(E_index, L_index);
        ++zero_counter;
    }

    void process_checkpoints(std::vector<MuonCheckpoint<T> > const & checkpoints) {
        for(auto checkpoint : checkpoints) {
            if(previous_energy < energy_grid.back()) {
                break;
            }
            process_checkpoint(checkpoint);
        }
    }
};

} // namepsace muons

