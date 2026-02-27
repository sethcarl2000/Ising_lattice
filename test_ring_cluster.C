#include "IsingLattice.h"
#include "IsingLattice.cxx"
#include <stdio.h>
#include <TStopwatch.h> 
#include <TH2D.h>
#include <TCanvas.h>
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

int test_ring_cluster(const char* path_outfile="histos_Metropolis.root")
{
    //the output file where we will save our data
    auto file = new TFile(path_outfile, "RECREATE"); 

    /// the nunber of spins for each testing ring 
    const int lattice_edge_size = 16; 
    /// number of spins
    const int n_spins = lattice_edge_size*lattice_edge_size; 
    /// the number of cluster updates between each temperature step
    const long long int annealing_clusters = 5000; 
    /// cluster updates between measurements
    const long long int flips_between_measurements = 5; 
    /// number of times to repeat a measurement at each step
    const int n_measurements_per_step = 30e3; 

    IsingLattice lattice(lattice_edge_size, 1.e6);

    lattice.HotStart(); 

    //temperature 
    const int n_bins_x = 80; 

    int energy_resolution = n_spins/4;

    const double T_range[] = {0., 6.}; 
    TH2D* hist_E = new TH2D("h_E", "E vs T;T/J;E/J*N", n_bins_x, T_range[0],T_range[1], energy_resolution, -2. -(0.5/((double)energy_resolution)), 0. +(0.5/((double)energy_resolution))); 
    
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

    vector<double> 
        E2(n_bins_x, 0.), 
        E(n_bins_x, 0.),
        flip_acceptance_rate(n_bins_x, 0.), 
        T_array(n_bins_x, 0.);

    const int n_measurements_per_step_thread = (n_measurements_per_step / n_threads);  
    for (int t=0; t<n_threads; t++) {

        //create a thread with its own Ising ring & histogram
        threads.emplace_back([
            n_threads, 
            &thread_progress, 
            &Print_thread_update,
            &mutex_histogram, &mutex_stdout, t, hist_E, 
            lattice_edge_size, n_spins,
            n_measurements_per_step_thread, 
            annealing_flips, flips_between_measurements,
            n_measurements_per_step, &E2, &E, &flip_acceptance_rate, &T_array
        ]{
            mutex_histogram.lock();
            auto hist_E_thread = (TH2D*)hist_E->Clone(Form("h_E_t%i",t)); 
            mutex_histogram.unlock(); 
            
            //start the ring at 'zero temp' 
            IsingLattice lattice(lattice_edge_size, 10000.); 
            lattice.ColdStart(t % 2 == 0);  //flip whether the spins start 'up' or 'down'

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
                lattice.FlipUpdate(annealing_flips); 
            
                double sum_E2_t{0.}, sum_E_t{0.}, sum_avg_accept_rate_t{0.};  
            
                //now, measure the energy a bunch of times. 
                for (long long int i=0; i<n_measurements_per_step_thread + (t < (n_measurements_per_step % n_threads) ? 1 : 0); i++) {

                    double E = lattice.Energy()/((double)n_spins); 

                    sum_E2_t += E*E; 
                    sum_E_t  += E; 

                    hist_E_thread->Fill( T, E ); 

                    //now, re-sweep and measure again
                    sum_avg_accept_rate_t += lattice.FlipUpdate(flips_between_measurements); 
                }

                mutex_histogram.lock(); 
                T_array[ix-1] = T; 
                E2[ix-1]                   += sum_E2_t              / ((double)n_measurements_per_step); 
                E[ix-1]                    += sum_E_t               / ((double)n_measurements_per_step); 
                flip_acceptance_rate[ix-1] += sum_avg_accept_rate_t / ((double)n_measurements_per_step); 
                mutex_histogram.unlock(); 
            }
            
            //once it's done, merge this thread's histogram with the 'main' one. 
            mutex_histogram.lock(); 
            thread_progress[t] = "done";
            Print_thread_update(); 
            for (int ix=1; ix<=hist_E->GetXaxis()->GetNbins(); ix++) {
                for (int iy=1; iy<=hist_E->GetYaxis()->GetNbins(); iy++) {
                    hist_E->Fill( 
                        hist_E->GetXaxis()->GetBinCenter(ix),
                        hist_E->GetYaxis()->GetBinCenter(iy),
                        hist_E_thread->GetBinContent(ix, iy) / ((double)n_measurements_per_step)
                    );
                }
            }
            mutex_histogram.unlock(); 

            delete hist_E_thread;
        }); 
    }
    for (auto& thread : threads) thread.join(); 


    //Energy 
    new TCanvas; 
    gStyle->SetPalette(kBird); 
    gStyle->SetOptStat(0);

    hist_E->DrawCopy("colz"); 

    auto g_E = new TGraph(n_bins_x, T_array.data(), E.data()); 
    g_E->SetMarkerStyle(kOpenCircle);
    g_E->SetMarkerSize(1.0);  
    g_E->Draw("SAME P"); 


    //Heat Capacity
    double pts_C[n_bins_x], pts_X[n_bins_x]; 
    
    auto vector_sum = [](const vector<double>& v) { double sum=0.; { for (double x : v) sum+=x; } return sum; }; 

    for (int i=0; i<n_bins_x; i++) {
        pts_C[i] = n_spins*(E2[i] - E[i]*E[i])/(T_array[i]*T_array[i]); 
    }

    new TCanvas; 
    auto g_C = new TGraph(n_bins_x, T_array.data(), pts_C); 
    g_C->SetTitle("Heat Capacity vs T;T/J;C/J");  
    g_C->SetMarkerStyle(kOpenCircle);
    g_C->Draw(); 
    

    //Average Cluster Size 
    new TCanvas; 
    auto g_FAR = new TGraph(n_bins_x, T_array.data(), flip_acceptance_rate.data()); 
    g_FAR->SetTitle("Average 'flip' acceptance rate vs. T;T/J;avg. flip acceptance rate"); 
    g_FAR->SetMarkerStyle(kOpenCircle); 
    g_FAR->Draw(); 

    //write the histogram / graphs to the TFile 
    hist_E->Write("hist_E_cluster"); 
    g_E->Write("graph_E_cluster");
    g_C->Write("graph_C_cluster");
    g_FAR->Write("graph_flipAcceptRate"); 
    file->Close(); 

    return 0; 
}
