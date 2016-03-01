#include "HartreeFock.hh" 

using namespace arma;

HartreeFock::HartreeFock(){
  for (int i =0; i < MAX_STATES; i++){
    for (int j = 0; j < MAX_STATES; j++){
      for (int k = 0; k < MAX_STATES; k++){
        for (int l = 0; l < MAX_STATES; l++){
          matrix_elements[i][j][k][l] = 0.0;
        }
      }
    }
  }

  num_states     = 0;
}

HartreeFock::~HartreeFock(){
}

void HartreeFock::Run(std::string const &in_mat_file_name, std::string const &in_sp_file_name,
                      std::string const &out_file_name){
  const double HAMILTONIAN_THRESHOLD = 10e-20;//Matrix elements only have precision of 10^-6
  int iteration = 0;

  std::cout << "Reading single particle states from file: " << in_sp_file_name << std::endl;

  num_states = ReadSingleParticleStates(in_sp_file_name);
  std::cout << "Found " << num_states << " states." << std::endl;
  if (num_states == 0){
    std::cout << "Failed to get correct number of states" << std::endl;
    return;
  }
  hamiltonian    = zeros(num_states,num_states);
  density_matrix = zeros(num_states,num_states);
  eigenvectors   = zeros(num_states,num_states);
  energies       = zeros(num_states,1);
  prev_energies  = zeros(num_states,1);
  for (int i =0; i < NUM_PARTICLES; i++){
      density_matrix(i,i) = 1;
  }

  mat harmonic_oscillator_energies = zeros(num_states,num_states);
  for (unsigned int i = 0; i < single_particle_states.size(); i++){
    harmonic_oscillator_energies(i,i) = single_particle_states.at(i).GetEnergy();
    prev_energies(i) = harmonic_oscillator_energies(i,i); 
  }
  std::cout << "Reading two body matrix elements from file: " << in_mat_file_name << std::endl;
  ReadMatrixElements(in_mat_file_name);


  double single_particle_potential = 0;
  std::cout << "Starting iterations..."<< std::endl;
  while(iteration < MAX_ITERATIONS){
    std::cout << "Completed iteration: " << iteration << "\xd"; 
    for (unsigned int alpha = 0; alpha < num_states; alpha++){
      for (unsigned int beta = 0; beta < num_states; beta++){
        single_particle_potential = 0;
        for (unsigned int gamma = 0; gamma < num_states; gamma++){
          for (unsigned int delta = 0; delta < num_states; delta++){
            single_particle_potential += density_matrix(gamma,delta)*matrix_elements[alpha][gamma][beta][delta];
          }//loop over delta
        }//loop over gamma
        
        //Hamiltonian is Hermitian
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
      std::cout << std::endl;
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
  for (unsigned int i = 0; i < num_states; i++){
    sum += fabs(diff(i)); 
  }
  sum /= num_states;
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

  for (unsigned int i = 0; i < num_states; i++){
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

unsigned int HartreeFock::ReadSingleParticleStates(std::string const &file_name){
  std::ifstream input_file;
  input_file.open(file_name.c_str());
  if (!input_file.is_open()){
    std::cout << "Failed to open matrix elements file" << std::endl;
    return 0;
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
  return single_particle_states.size();
}

void HartreeFock::FillDensityMatrix(){
//TODO: Check if density matrix is initialized correctly
  density_matrix = zeros(num_states,num_states);
  for (unsigned int gamma = 0; gamma < num_states; gamma++){
    for (unsigned int delta = 0; delta < num_states; delta++){
      for(unsigned int state = 0; state < NUM_PARTICLES; state++){
        density_matrix(gamma,delta) += eigenvectors(delta,state)*eigenvectors(gamma,state);
      }//state
    }//delta
  }//gamma
}


