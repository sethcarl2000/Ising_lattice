#ifndef IsingLattice_CXX
#define IsingLattice_CXX

#include "IsingLattice.h"
#include <stdio.h> 
#include <iostream> 
#include <cmath> 
#include <array> 
#include <sstream> 
#include <stdexcept> 
//_________________________________________________________________________________________
//_________________________________________________________________________________________
double& SpinLattice_t::operator()(int ix, int iy)
{
    if (iy >= edge_size) { iy += -edge_size; } if (iy < 0) { iy += +edge_size; }
    if (ix >= edge_size) { ix += -edge_size; } if (ix < 0) { ix += +edge_size; }

    return spins[ix*edge_size + iy]; 
}
//_________________________________________________________________________________________
IsingLattice::IsingLattice(const int _lattice_edge_size, double _Beta)
    : lattice_edge_size{_lattice_edge_size}, lattice{_lattice_edge_size, {}}, Beta{_Beta} 
{
    lattice.spins.reserve(lattice_edge_size*lattice_edge_size); 

    //initialize random-number generators
    std::random_device rd; 
    twister = std::mt19937(rd()); 
    real_dist = std::uniform_real_distribution<double>(0., 1.); 
    index_dist = std::uniform_int_distribution<int>(0, lattice_edge_size-1); 

    //start with a cold start
    ColdStart(); 
}
//_________________________________________________________________________________________
IsingLattice::IsingLattice(const SpinLattice_t& _lattice, double _Beta)
    : lattice_edge_size{_lattice.edge_size}, lattice{}, Beta{_Beta} 
{
    //copy the array of spins we were given
    lattice.spins.clear(); 
    std::copy( _lattice.spins.begin(), _lattice.spins.end(), lattice.spins.begin() ); 

    //initialize random-number generators
    std::random_device rd; 
    twister = std::mt19937(rd()); 
    real_dist = std::uniform_real_distribution<double>(0., 1.); 
    index_dist = std::uniform_int_distribution<int>(0, lattice_edge_size-1); 
}
//_________________________________________________________________________________________
void IsingLattice::HotStart()
{
    //initialize spins, then randomly pick all spins
    ColdStart(); 
    for (auto& S : lattice.spins) S = (Rand() > 0.5) ? +1. : -1.;
}
//_________________________________________________________________________________________
void IsingLattice::ColdStart(bool spin_up)
{
    //assign all spins to either up or down
    lattice.spins = std::vector<double>(lattice_edge_size*lattice_edge_size, spin_up ? +1. : -1.); 
}
//_________________________________________________________________________________________
double IsingLattice::FlipUpdate(const long long int n_steps) 
{
    //we're going to sweep over the whole lattice, using the metropolis alg. to deterime if we should flip spins. 

    long long int n_accepted_steps=0; 

    long long int i_step=0; 
    do {
        //pick a random spin
        int ix{RandIndex()}, iy{RandIndex()}; 

        //get the 'central' spin
        double& S  = lattice(ix, iy); 
        
        //this is the 'molecular field'; it's just this spin's neighborhood of spins. 
        double dE = 2.*lattice(ix, iy)*(
            lattice(ix +1, iy   ) +  
            lattice(ix   , iy +1) + 
            lattice(ix -1, iy   ) + 
            lattice(ix   , iy -1)
        ); 
   
        if (dE < 0. || Rand() < std::exp(-Beta*dE) ) { S *= -1.; ++n_accepted_steps; }
        
    } while (++i_step < n_steps); 

    return ((double)n_accepted_steps)/((double)n_steps); 
}
//_________________________________________________________________________________________
double IsingLattice::ClusterUpdate(const long long int n_clusters)
{
    long long int average_clus_size=0; 

    //probability of adding a new block to the cluster
    const double p = 1. - std::exp( -2.*Beta );

    int cluster_size =0; 
    LatticeSite_t cluster[lattice_edge_size*lattice_edge_size]; 

    /// @return 'true' if the given site is inside the cluster, 'false' otherwise
    auto inside_cluster = [&cluster, &cluster_size](const LatticeSite_t& site) {

        for (int i=0; i<cluster_size; i++) { if (site == cluster[i]) return true; }
        return false; 
    }; 
    //___________________________________________________________________________

    long long int i_cluster=0; 
    do {
        //pick random site & add it to the cluster
        LatticeSite_t site{ RandIndex(), RandIndex() }; 

        cluster[0] = site; ++cluster_size;

        const double spin = lattice(site); 

        int last_spin_index{0}, next_spin_index{cluster_size};  

        do {
            //we're adding new spins to the cluster in layers, like an onion. 
            for (int i=last_spin_index; i<next_spin_index; i++) { const auto& site_i = cluster[i]; 

                //check all 4 neighboring spins, to see if they're in cluster or are the same spin
                for (auto site_j : {
                    LatticeSite_t{site_i.x + 1, site_i.y    }, 
                    LatticeSite_t{site_i.x    , site_i.y + 1}, 
                    LatticeSite_t{site_i.x - 1, site_i.y    }, 
                    LatticeSite_t{site_i.x    , site_i.y - 1}
                }) {
                    if      (site_j.x < 0)                  { site_j.x += lattice_edge_size; } 
                    else if (site_j.x >= lattice_edge_size) { site_j.x -= lattice_edge_size; }
                    
                    if      (site_j.y < 0)                  { site_j.y += lattice_edge_size; } 
                    else if (site_j.y >= lattice_edge_size) { site_j.y -= lattice_edge_size; }

                    if (spin*lattice(site_j) < 0.) continue; //skip this site if it's a different sign then the cluster
                    
                    if (inside_cluster(site_j)) continue;    //skip this site if it's already inside the cluster
            
                    //this spin is the correct sign, and it's not already in the cluster 
                    // so this spin is eligible to be added. 
                    if (Rand() < p) cluster[cluster_size++] = site_j; 
                    
                }; 
            } 
            //printf("cluster size: %i\n", cluster_size); 
            
            //now that we've added all the spins we can find, let's start a new round. 
            last_spin_index = next_spin_index; 
            next_spin_index = cluster_size; 

        } while (next_spin_index - last_spin_index > 0); //break the loop if we didn't add any spins

        //flip all the spins in the cluster
        for (int i=0; i<cluster_size; i++) { lattice(cluster[i]) *= -1.; }
        average_clus_size += cluster_size; 
        cluster_size =0;   

    } while (++i_cluster < n_clusters);  

    return ((double)average_clus_size)/((double)n_clusters); 
}
//_________________________________________________________________________________________
double IsingLattice::Energy()
{
    double E = 0.; 
    //scan each row as it it were a 1d ising model.

    //first, scan all columns 
    for (int ix=0; ix<lattice_edge_size; ix++) {
        
        double S_last = lattice(ix, -1); 
        for (int iy=0; iy<lattice_edge_size; iy++) {
            double S = lattice(ix, iy); 
            E += S*S_last; 
            S_last = S; 
        }
    }
    //next, scan all rows 
    for (int iy=0; iy<lattice_edge_size; iy++) {
        
        double S_last = lattice(-1, iy); 
        for (int ix=0; ix<lattice_edge_size; ix++) {
            double S = lattice(ix, iy); 
            E += S*S_last; 
            S_last = S; 
        }
    }

    return -1.*E; 
}
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
void IsingLattice::Print()
{
    for (int ix=0; ix<lattice_edge_size; ix++) {
        for (int iy=0; iy<lattice_edge_size; iy++) {
            std::cout << (lattice(ix,iy) > 0. ? "+" : "-");
        }        
        std::cout << "\n";
    } 
}
//_________________________________________________________________________________________

#endif