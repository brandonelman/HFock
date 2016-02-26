
class State {
  public:
    State();
    ~State();

    State(int state_index, int n, int l, int j2, int mj2, int tz2);
    

    //Calculates harmonic oscillator energy for number of nodes N, 
    //angular momentum quantum number L, and given  frequency W.
    //Returns (2*N+L+3/2)hbar*W
    double GetEnergy();

    

    int state_index;
    int n;//number of nodes in wave function (starting from 0);
    int l;//orbital angular momentum

    //These are all twice the value of the quantity
    //because they're all half-integers
    int j2;//total angular momentum
    int mj2;//total angular momentum z-projection
    int tz2;//isospin z-projection  
    void Print();
};
