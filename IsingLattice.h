
#ifndef IsingLattice_H
#define IsingLattice_H 

#include <random>
#include <vector>

struct LatticeSite_t {
    int x, y; 
    inline bool operator==(const LatticeSite_t&rhs) const {
        return (rhs.x==x) && (rhs.y==y); 
    }
};

/// List of all spins
/// Spins are indexed from right-to-left, top-to-bottom. 
/// for example, in a lattice of 100x100 spins, the spin at site {49,20} 
/// (wtith {0,0} being the top-leftmost spin) would have array index = 20*100 + 49; 
struct SpinLattice_t {
    const int edge_size; 
    std::vector<double> spins; ///list of all spins.
    inline double& operator()(int ix, int iy); 
    inline double& operator()(const LatticeSite_t& s) { return (*this)(s.x, s.y); }
};

/// @brief Class which implements an Ising ring
class IsingLattice {
private: 

    SpinLattice_t lattice; 

    const int lattice_edge_size; //length of the lattice's edge 
    
    double Beta; /// The thermodynamic Beta of the resivoir (specifically, this is the ratio Beta/J)  

    std::mt19937 twister; 
    std::uniform_real_distribution<double> real_dist; 
    std::uniform_int_distribution<int> index_dist; 

    inline double Rand() { return real_dist(twister); }

    inline int RandIndex() { return index_dist(twister); }

public: 

    /// @brief Default constructor for a 'fresh' ring.
    /// @param lattice_edge_size edge size for the square lattice. 
    /// @param Beta 1/T, in units of 1/J. (energy scale)
    IsingLattice(const int lattice_edge_size, double Beta); 

    IsingLattice(const SpinLattice_t& _lattice, double Beta); 

    ~IsingLattice() { lattice.spins.clear(); }; 

    void SetBeta(double _beta) { Beta=_beta; }

    /// @brief Randomize the lattice; each spin is random
    void HotStart(); 

    /// @brief Set the lattice to spin all in one direction. 
    /// @param spin_up true = Start the lattice all spin-up (lower energy when H > 0)
    void ColdStart(bool spin_up=true); 

    /// @brief Consider a fixed number of individual spin-flips at the given temperature / H. 
    /// @param n_flips number of individual flips to consider 
    /// @return fraction of events which are accepted
    double FlipUpdate(const long long int n_flips);

    /// @brief Consider a number of cluster updates using the wolf algorithm 
    /// @param n_clusters number of cluster flips to perform
    /// @return average cluster size
    double ClusterUpdate(const long long int n_clusters); 

    /// @return Pointer to array of spins
    SpinLattice_t SetSpins() const { return lattice; } 

    /// @return Returns energy (in units of J)
    double Energy(); 

    /// @return lattice magnetization
    double Magnetization(); 

    /// @brief compute energy & magnetization 
    /// @param E variable in which to store the energy (units: J)
    /// @param M variable in which to store the magnetization
    void GetEnergyMagnetization(double& E, double& M); 

    /// @brief Initialize the spin-array with a particular value
    /// @param _lattice lattice of spins
    void SetSpins(const SpinLattice_t& _lattice) { 
        lattice.spins.clear(); 
        std::copy(_lattice.spins.begin(), _lattice.spins.end(), lattice.spins.begin() );
    }
    
    /// @brief print out all spins
    void Print(); 
}; 

#endif 