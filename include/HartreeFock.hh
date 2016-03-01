#ifndef __HARTREE_H__
#define __HARTREE_H__


#include <iostream>
#include <armadillo>
#include <fstream>
#include <algorithm>
#include "State.hh"

class HartreeFock {
  public:
    HartreeFock();
    ~HartreeFock();

    void Run(std::string const &in_mat_file_name, std::string const &in_sp_file_name,
             std::string const &out_file_name, int num_particles);
    unsigned int num_states;
//    static const int NUM_PARTICLES = 8;
    static const int MAX_ITERATIONS = 50;
    static const int MAX_STATES = 100;

  private:
    //Probably not necessary
    //double CalcHarmonicOscillatorPotential();
    //double CalcWoodsSaxon();

    bool IsConverged() const;

    //Write energies and eigenvectors to file
    void SaveToFile(std::string const & file_name, arma::mat &, std::vector<State>&) const;

    //Fills the 80x80 matrix_elements array
    void ReadMatrixElements(std::string const &in_mat_file_name);
    unsigned int ReadSingleParticleStates(std::string const &in_sp_file_name);
    void FillDensityMatrix(std::vector<int>, int num_particles);

    std::vector<State> single_particle_states;
    arma::mat hamiltonian; 
    arma::mat density_matrix;
    arma::vec energies;
    arma::vec prev_energies;//Need these to know when converged
    arma::mat eigenvectors;

    //Note given matrix elements are already anti-symmetrized
    double matrix_elements[MAX_STATES][MAX_STATES][MAX_STATES][MAX_STATES];

};
//ring size 7.5

#endif
