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
#include <thread> 
#include <vector> 
#include <mutex> 

    
using namespace std; 

namespace {
    
}; 

int test_ring()
{
    /// the nunber of spins for each testing ring 
    const int lattice_edge_size = 16; 
    /// the number of steps from a 'hot start' to get the rings to work
    const long long int annealing_steps = lattice_edge_size*6e2; 
    /// steps betewen a measurement of each quantity
    const long long int measurement_steps = lattice_edge_size*lattice_edge_size*20; 
    /// number of times to repeat a measurement at each step
    const int n_measurements_per_step = 16e4; 

    IsingLattice lattice(lattice_edge_size, 0.01);

    lattice.HotStart(); 

    cout << "\n\nstart --- " << endl; 
    lattice.Print();

    for (int i=0; i<10; i++) {

        lattice.ClusterUpdate(1e3); 

        cout << "\n\niteration " << i << endl; 
        lattice.Print(); 
    }

    //temperature 
    const int n_bins_x = 80; 

    int energy_resolution = lattice_edge_size*lattice_edge_size/4;

    const double T_range[] = {0., 6.}; 
    TH2D* hist_E = new TH2D("h_E", "E vs T;T/J;E/J*N", n_bins_x, T_range[0],T_range[1], energy_resolution, -1. -(0.5/((double)energy_resolution)), 0. +(0.5/((double)energy_resolution))); 
    
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

    vector<double> 
        E2(n_bins_x, 0.), 
        E(n_bins_x, 0.); 

    const int n_measurements_per_step_thread = (n_measurements_per_step / n_threads);  
    for (int t=0; t<n_threads; t++) {

        //create a thread with its own Ising ring & histogram
        threads.emplace_back([
            n_threads, 
            &thread_progress, 
            &Print_thread_update,
            &mutex_histogram, &mutex_stdout, t, hist_E, 
            n_measurements_per_step_thread, 
            annealing_steps, 
            measurement_steps, &E2, &E, &M2, &M, &T_array
        ]{
            mutex_histogram.lock();
            auto hist_E_thread = (TH2D*)hist_E->Clone(Form("h_E_t%i",t)); 
            auto hist_M_thread = (TH2D*)hist_M->Clone(Form("h_M_t%i",t));
            mutex_histogram.unlock(); 
            
            //start the ring at 'zero temp' 
            IsingLattice ring(n_spins, 10000., 0.); 
            ring.ColdStart(t % 2 == 0);  //flip whether the spins start 'up' or 'down'

            //create & anneal our ising ring
            auto x_axis = hist_E_thread->GetXaxis(); 
            
            //go to each X-bin, and use the corresponding temperature to create an ising ring at that temp. 
            for (int ix=1; ix<=x_axis->GetNbins(); ix++) {

                double T = x_axis->GetBinCenter(ix); 

                if ((ix-1) % (x_axis->GetNbins()/10) == 0 && ix>1) {
                    mutex_stdout.lock();
                    thread_progress[t] = Form("%.2f", ((double)ix-1)/((double)x_axis->GetNbins())); 
                    Print_thread_update(); 
                    mutex_stdout.unlock(); 
                }

                //now, anneal this for a while at this temp. 
                ring.SetBeta(1./T); 
                ring.Sweep(annealing_steps); 
            
                double sum_E2_t{0.}, sum_E_t{0.}, sum_M2_t{0.}, sum_M_t{0.}; 
            
                //now, measure the energy a bunch of times. 
                for (long long int i=0; i<n_measurements_per_step_thread + (t < (n_measurements_per_step % n_threads) ? 1 : 0); i++) {

                    double E = ring.Energy()/((double)n_spins); 
                    double M = ring.Magnetization()/((double)n_spins); 

                    sum_E2_t += E*E; 
                    sum_E_t  += E; 
                    sum_M2_t += M*M; 
                    sum_M_t  += M; 

                    hist_E_thread->Fill( T, E ); 
                    hist_M_thread->Fill( T, M ); 

                    //now, re-sweep and measure again
                    ring.Sweep(measurement_steps); 
                }

                mutex_histogram.lock(); 
                T_array[ix-1] = T; 
                E2[ix-1] += sum_E2_t    / ((double)n_measurements_per_step); 
                E[ix-1]  += sum_E_t     / ((double)n_measurements_per_step); 
                M2[ix-1] += sum_M2_t    / ((double)n_measurements_per_step); 
                M[ix-1]  += sum_M_t     / ((double)n_measurements_per_step); 
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
            for (int ix=1; ix<=hist_M->GetXaxis()->GetNbins(); ix++) {
                for (int iy=1; iy<=hist_M->GetYaxis()->GetNbins(); iy++) {
                    hist_M->Fill(
                        hist_M->GetXaxis()->GetBinCenter(ix),
                        hist_M->GetYaxis()->GetBinCenter(iy),
                        hist_M_thread->GetBinContent(ix, iy) / ((double)n_measurements_per_step)
                    );
                }   
            }
            mutex_histogram.unlock(); 

            delete hist_E_thread;  
            delete hist_M_thread; 
        }); 
    }

    for (auto& thread : threads) thread.join(); 



    double x_test = 0.001;
    cout << tf1_X->Eval(x_test) << endl; 

    new TCanvas; 
    gStyle->SetPalette(kBird); 
    gStyle->SetOptStat(0);

    hist_E->Draw("colz"); 

    auto g_E = new TGraph(n_bins_x, T_array.data(), E.data()); 
    g_E->SetMarkerStyle(kOpenCircle);
    g_E->SetMarkerSize(1.0);  
    g_E->Draw("SAME P"); 

    tf1_E->SetLineStyle(kSolid); 
    tf1_E->Draw("SAME"); 

    new TCanvas; 
    hist_M->Draw("colz"); 

    auto g_M = new TGraph(n_bins_x, T_array.data(), M.data()); 
    g_M->SetMarkerStyle(kOpenSquare); 
    g_M->SetMarkerSize(0.3);  
    g_M->Draw("SAME P"); 

    double pts_C[n_bins_x], pts_X[n_bins_x]; 
    
    auto vector_sum = [](const vector<double>& v) { double sum=0.; { for (double x : v) sum+=x; } return sum; }; 

    for (int i=0; i<n_bins_x; i++) {
        pts_C[i] = n_spins*(E2[i] - E[i]*E[i])/(T_array[i]*T_array[i]); 
        pts_X[i] = n_spins*(M2[i] - M[i]*M[i])/T_array[i]; 
    }

    new TCanvas; 
    tf1_C->SetLineStyle(kSolid); 
    tf1_C->Draw(); 

    auto g_C = new TGraph(n_bins_x, T_array.data(), pts_C); 
    g_C->SetTitle("Heat Capacity vs T;T/J;C/J");  
    g_C->SetMarkerStyle(kOpenCircle);
    g_C->Draw("SAME P"); 

    new TCanvas; 
    tf1_X->SetLineStyle(kDashed); 
    tf1_X->Draw(); 

    auto g_X = new TGraph(n_bins_x, T_array.data(), pts_X); 
    g_X->SetTitle("#chi vs T;T/J;#chi"); 
    g_X->SetMarkerStyle(kOpenCircle);
    g_X->Draw("SAME P"); 

    return 0; 
}
