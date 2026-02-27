#include "IsingLattice.h"
#include "IsingLattice.cxx"
#include <stdio.h>
#include <TStopwatch.h> 
#include <TH2D.h>
#include <TCanvas.h>
#include <TStopwatch.h> 
#include <TStyle.h>
#include <TAxis.h> 
#include <TPad.h>
#include <TGraph.h> 
#include <TF1.h> 
#include <TFile.h> 
#include <thread> 
#include <vector> 
#include <mutex> 

using namespace std; 

namespace {
    /// @brief Add data from each bin of 'source' histogram into 'target' histogram
    void Add_histogram(TH2D* source, TH2D *target){
        for (int ix=1; ix<=source->GetXaxis()->GetNbins(); ix++) {
            for (int iy=1; iy<=source->GetYaxis()->GetNbins(); iy++) {
                target->Fill( 
                    source->GetXaxis()->GetBinCenter(ix),
                    source->GetYaxis()->GetBinCenter(iy),
                    source->GetBinContent(ix, iy)
                );
            }
        }
        return; 
    }
    
    //_________________________________________________________________________________
}

int test_ring_cluster(const char* path_outfile="histos_Cluster_test.root")
{
    //the output file where we will save our data
    auto file = new TFile(path_outfile, "RECREATE"); 

    /// the nunber of spins for each testing ring 
    const int lattice_edge_size = 32; 
    /// number of spins
    const int n_spins = lattice_edge_size*lattice_edge_size; 
    /// the number of cluster updates between each temperature step
    const long long int annealing_clusters = 2000; 
    /// cluster updates between measurements
    const long long int clusters_between_measurements = 4; 
    /// number of times to repeat a measurement at each step
    const int n_measurements_per_step = 3e3; 

    IsingLattice lattice(lattice_edge_size, 1.e6);

    lattice.HotStart(); 

    //temperature 
    const int n_bins_x = 80; 

    int energy_resolution = n_spins/2;

    const double T_range[] = {1.5, 3.5}; 
    TH2D* hist_E = new TH2D("h_E", "E vs T;T/J;E/J*N", 
        n_bins_x, T_range[0],T_range[1], 
        energy_resolution, -2.*((double)n_spins) - 0.5,  0. + 0.5
    ); 
    
    TH2D* hist_M = new TH2D("h_M", "M vs T;T/J; M/N", 
        n_bins_x, T_range[0], T_range[1],
        n_spins, -((double)n_spins) - 0.5, ((double)n_spins) + 0.5 
    );
    
    std::mutex mutex_histogram, mutex_stdout; 

    vector<thread> threads; 

    const int n_threads = std::thread::hardware_concurrency(); 

    std::vector<string> thread_progress(n_threads, "0.00"); 
    auto Print_thread_update = [&thread_progress]() {
        printf("\r t: "); 
        for (int i=0; i<thread_progress.size(); i++) {
            printf(" %s", thread_progress[i].data()); 
        }
        cout << flush; 
    };
    Print_thread_update(); 

    TStopwatch timer; 

    const int n_measurements_per_step_thread = (n_measurements_per_step / n_threads);  
    for (int t=0; t<n_threads; t++) {

        //create a thread with its own Ising ring & histogram
        threads.emplace_back([
            n_threads, 
            &thread_progress, 
            &Print_thread_update,
            &mutex_histogram, &mutex_stdout, t, hist_E, hist_M,  
            lattice_edge_size, n_spins,
            n_measurements_per_step_thread, 
            annealing_clusters, clusters_between_measurements,
            n_measurements_per_step
        ]{
            mutex_histogram.lock();
            auto hist_E_thread = (TH2D*)hist_E->Clone(Form("h_E_t%i",t)); 
            auto hist_M_thread = (TH2D*)hist_M->Clone(Form("h_M_t%i",t)); 
            mutex_histogram.unlock(); 
            
            //start the ring at 'zero temp' 
            IsingLattice lattice(lattice_edge_size, 10000.); 
            lattice.HotStart();  //flip whether the spins start 'up' or 'down'

            //create & anneal our ising ring
            auto x_axis = hist_E_thread->GetXaxis(); 
            
            //go to each X-bin, and use the corresponding temperature to create an ising ring at that temp. 
            for (int ix=1; ix<=x_axis->GetNbins(); ix++) {

                double T = x_axis->GetBinCenter(ix); 

                mutex_stdout.lock();
                thread_progress[t] = Form("%.2f", ((double)ix-1)/((double)x_axis->GetNbins())); 
                Print_thread_update(); 
                mutex_stdout.unlock(); 
            
                //now, anneal this for a while at this temp. 
                lattice.SetBeta(1./T); 
                lattice.ClusterUpdate(annealing_clusters); 
            
                //now, measure the energy a bunch of times. 
                for (long long int i=0; i<n_measurements_per_step_thread + (t < (n_measurements_per_step % n_threads) ? 1 : 0); i++) {

                    double E, M; 
                    
                    lattice.GetEnergyMagnetization(E, M); 

                    hist_E_thread->Fill( T, E ); 
                    hist_M_thread->Fill( T, M ); 

                    lattice.ClusterUpdate(clusters_between_measurements); 
                }
            }
            
            //once it's done, merge this thread's histogram with the 'main' one. 
            mutex_histogram.lock(); 
            thread_progress[t] = "done";
            Print_thread_update(); 
            Add_histogram(hist_E_thread, hist_E); 
            Add_histogram(hist_M_thread, hist_M); 
            mutex_histogram.unlock(); 

            delete hist_E_thread;
        }); 
    }
    for (auto& thread : threads) thread.join(); 

    const double elapsed_time = timer.RealTime(); 

    printf(
        "\nTime elapsed for %i spins: \n"
        "       overall time              %.3f s \n"
        "       time per spin             %.3f ms / spin\n"
        "       time per T-step * spin    %.3f ms / spin * temp\n", 
        n_spins, 
        elapsed_time, 
        1000.*elapsed_time /((double)n_spins),
        1000.*elapsed_time /((double)n_bins_x*n_spins)
    ); 

    //Energy 
    new TCanvas; 
    gStyle->SetPalette(kBird); 
    gStyle->SetOptStat(0);

    hist_E->DrawCopy("colz"); 
    

    //Magnetization 
    new TCanvas; 
    hist_M->DrawCopy("col"); 

    //write the histogram / graphs to the TFile 
    hist_E->Write("hist_E_cluster"); 
    
    file->Close(); 

    return 0; 
}
