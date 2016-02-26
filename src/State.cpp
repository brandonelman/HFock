#include "State.hh"
#include <iostream>
#include <vector>
#include <iomanip>
State::State(){
  state_index = 0;
  n = 1;
  l = 0; j2 = 1;
  mj2 = 1;
  tz2 = 1;
}

State::~State(){
}

State::State(int state_index, int n, int l, int j2, int mj2, int tz2){
  this->state_index = state_index; 
  this->n =  n;
  this->l = l;
  this->j2 = j2;
  this->mj2 = mj2;
  this->tz2 = tz2;
}
double State::GetEnergy(){
  const double HBARW = 10.;//MeV
  return (2*n+l+(3./2.))*HBARW;//MeV
}

void State::Print(){
  std::vector<std::string> shell_names;
  shell_names.push_back("s");
  shell_names.push_back("p");
  shell_names.push_back("d");
  shell_names.push_back("f");
  shell_names.push_back("g");
  shell_names.push_back("h");

  std::cout << std::setw(2) << state_index << "\t" << std::setw(1) << n << shell_names.at(l)
            << j2 << "/2\t"<< mj2 << "/2\t" << tz2 << "/t" << std::endl;
}
