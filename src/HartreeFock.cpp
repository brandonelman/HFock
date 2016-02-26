#include "HartreeFock.hh" 

using namespace arma;

HartreeFock::HartreeFock(){
  for (int i =0; i < NUM_STATES; i++){
    for (int j = 0; j < NUM_STATES; j++){
      for (int k = 0; k < NUM_STATES; k++){
        for (int l = 0; l < NUM_STATES; l++){
          matrix_elements[i][j][k][l] = 0.0;
        }
      }
    }
  }

  hamiltonian    = zeros(NUM_STATES,NUM_STATES);
  density_matrix = zeros(NUM_STATES,NUM_STATES);
  eigenvectors   = zeros(NUM_STATES,NUM_STATES);
  energies       = zeros(NUM_STATES,1);
  prev_energies  = zeros(NUM_STATES,1);
  for (int i =0; i < 16; i++){
      density_matrix(i,i) = 1;
  }
}

HartreeFock::~HartreeFock(){
}

void HartreeFock::Run(std::string const &in_mat_file_name, std::string const &in_sp_file_name,
                      std::string const &out_file_name){
  const double HAMILTONIAN_THRESHOLD = 10e-6;
  int iteration = 0;
  mat harmonic_oscillator_energies = zeros(NUM_STATES,NUM_STATES);
  std::cout << "Reading single particle states from file: " << in_sp_file_name << std::endl;
  ReadSingleParticleStates(in_sp_file_name);
  if (single_particle_states.size() != NUM_STATES){
    std::cout << "ERROR: Failed to get correct number of single particle states!" << std::endl;
    std::cout << "Expect: " << NUM_STATES << "\tFound: " << single_particle_states.size() << std::endl;
    return;
  }
  for (unsigned int i = 0; i < single_particle_states.size(); i++){
    harmonic_oscillator_energies(i,i) = single_particle_states.at(i).GetEnergy();
    prev_energies(i) = harmonic_oscillator_energies(i,i); 
  }
  std::cout << "Reading two body matrix elements from file: " << in_mat_file_name << std::endl;
  ReadMatrixElements(in_mat_file_name);


  double single_particle_potential = 0;
  std::cout << "Starting iterations..."<< std::endl;
  while(iteration < MAX_ITERATIONS){
    std::cout << "Completed iteration: " << iteration << std::endl; 
    for (int alpha = 0; alpha < NUM_STATES; alpha++){
      for (int beta = 0; beta < NUM_STATES; beta++){
        single_particle_potential = 0;
        for (int gamma = 0; gamma < NUM_STATES; gamma++){
          for (int delta = 0; delta < NUM_STATES; delta++){
            single_particle_potential += density_matrix(gamma,delta)*matrix_elements[alpha][gamma][beta][delta];
          }//loop over delta
        }//loop over gamma
        
        //Hamiltonian is Hermitian
        if (alpha == 64 && beta == 68){
          std::cout << "single particle potential: " << single_particle_potential << std::endl;
          std::cout << "harmonic oscillator energies: " << harmonic_oscillator_energies(alpha,beta) << "\n\n\n";
        }
        
        hamiltonian(alpha,beta) = harmonic_oscillator_energies(alpha,beta) + single_particle_potential;
        hamiltonian(beta,alpha) = harmonic_oscillator_energies(alpha,beta) + single_particle_potential;
        if (fabs(hamiltonian(alpha,beta)) < HAMILTONIAN_THRESHOLD){
          hamiltonian(alpha,beta) = 0;
          hamiltonian(beta,alpha) = 0;
        }

      }//loop over beta
    }//loop over alpha
    eig_sym(energies,eigenvectors,hamiltonian);
    if  (IsConverged()){
      std::cout << "Converged after " << iteration<< "!" << std::endl;
      break;
    }
    FillDensityMatrix();
    prev_energies = energies;
    iteration++;
  }//while iteration < MAX_ITERATIONS and not converged
  SaveToFile(out_file_name, harmonic_oscillator_energies, single_particle_states);
}


bool HartreeFock::IsConverged() const {
  const double THRESHOLD = pow(10,-8);
  vec diff = energies - prev_energies;
  double sum = 0;
  for (int i = 0; i < NUM_STATES; i++){
    sum += fabs(diff(i)); 
  }
  sum /= NUM_STATES;
//std::cout << "sum = " << sum << std::endl;
//std::cout << "diff = " << diff << std::endl;

  if (sum < THRESHOLD){
    return true;
  }
  else {
    return false;
  }
}

void HartreeFock::SaveToFile(std::string const &file_name, mat &harmonic_oscillator_energies, std::vector<State> &states) const{
  std::ofstream proton_out_file;
  std::ofstream neutron_out_file;
  std::string proton_file_name = "proton_"+file_name;
  std::string neutron_file_name = "neutron_"+file_name;

  proton_out_file.open(proton_file_name.c_str());
  if (!proton_out_file.is_open()){
    std::cout << "ERROR: Failed to open output file!" << std::endl;
  }
  neutron_out_file.open(neutron_file_name.c_str());
  if (!neutron_out_file.is_open()){
    std::cout << "ERROR: Failed to open output file!" << std::endl;
  }

  proton_out_file << "HO Energies\tHF Energies\tOrbital#\tnlj\tmj\ttz" << std::endl;
  neutron_out_file << "HO Energies\tHF Energies\tOrbital#\tnlj\tmj\ttz" << std::endl;
  std::vector<std::string> shell_names;// = {"s","p","d","f","g"};
  shell_names.push_back("s");
  shell_names.push_back("p");
  shell_names.push_back("d");
  shell_names.push_back("f");
  shell_names.push_back("g");
  shell_names.push_back("h");

  for (int i = 0; i < NUM_STATES; i++){
    if (states.at(i).tz2 < 0){
      proton_out_file << harmonic_oscillator_energies(i,i) <<"\t"<<energies(i) << "\t" <<states.at(i).state_index<<"\t"
                      << states.at(i).n << shell_names.at(states.at(i).l) << states.at(i).j2 << "/2\t" 
                      << states.at(i).mj2 <<"/2\t" << states.at(i).tz2 << "/2" << std::endl; 
    }
    else if (states.at(i).tz2 > 0){
      neutron_out_file << harmonic_oscillator_energies(i,i) <<"\t"<<energies(i) << "\t" <<states.at(i).state_index<<"\t" 
                        << states.at(i).n << shell_names.at(states.at(i).l) << states.at(i).j2 << "/2\t" 
                        << states.at(i).mj2 <<"/2\t" << states.at(i).tz2 << "/2" << std::endl; 
    }
  }

  std::ofstream mat_file("out_mat_file.dat");
  mat_file << hamiltonian;
//for (int i = 0; i < NUM_STATES; i++){
//  for (int j = 0; j < NUM_STATES; j++){
//    mat_file << hamiltonian<< std::endl;
//  }
//}
}

void HartreeFock::ReadMatrixElements(std::string const &file_name){
  std::ifstream input_file;
  input_file.open(file_name.c_str());
  if (!input_file.is_open()){
    std::cout << "ERROR: Failed to open matrix elements file" << std::endl;
    return;
  }

  int a = -1;
  int b = -1;
  int c = -1; 
  int d = -1;
  double matrix_element =-1;
  std::string line;
  std::cout << "Starting to parse matrix elements..." << std::endl;
  while(std::getline(input_file,line)){
    sscanf(line.c_str(), "%d %d %d %d %lf", &a,&b,&c,&d,&matrix_element); 
    //minus one because want to be referenced from 0
    matrix_elements[a-1][b-1][c-1][d-1] = matrix_element;   
  }
  std::cout << "Finished parsing matrix elements" << std::endl;
  return;
}

void HartreeFock::ReadSingleParticleStates(std::string const &file_name){
  std::ifstream input_file;
  input_file.open(file_name.c_str());
  if (!input_file.is_open()){
    std::cout << "Failed to open matrix elements file" << std::endl;
    return;
  }

  int state_index =-1;
  int n = -1;
  int l = -1; 
  int j2 = -1;
  int mj2 = -1;
  int tz2 = -1;
  std::string line;
  while(std::getline(input_file,line)){
    sscanf(line.c_str(), "Orbit number: %d %d %d %d %d %d", &state_index,&n,&l,&j2,&mj2,&tz2); 
    State state(state_index,n,l,j2,mj2,tz2);
    state.Print();
    single_particle_states.push_back(state);
  }
  return;
}

void HartreeFock::FillDensityMatrix(){
  density_matrix = zeros(NUM_STATES,NUM_STATES);
//density_matrix = eigenvectors*eigenvectors.t();
  for (int gamma = 0; gamma < NUM_STATES; gamma++){
    for (int delta = 0; delta < NUM_STATES; delta++){
      for(int state = 0; state < 16; state++){
        density_matrix(gamma,delta) += eigenvectors(delta,state)*eigenvectors(gamma,state);
      }//state
    }//delta
  }//gamma
}


