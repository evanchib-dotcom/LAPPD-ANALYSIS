#include <iostream>
#include <random>
#include <vector>
#include <fstream>

class electron {
public:
    double position;
    double velocity;
    double  acceleration;
   
 electron () {
     std::random_device rd;                          // Non-deterministic seed
    std::mt19937 gen(rd());                         // Mersenne Twister engine
    std::uniform_int_distribution<> dist(-5, 5); 
    this -> position = 0 ;
    this -> velocity = (dist(gen));
    this -> acceleration = 0.1;
 }

    void updatePosition(double timestep) {
        this -> position = position + (velocity * timestep);
    }

    void updateVelocity (double timestep) {
        this -> velocity = velocity + (acceleration * timestep);
    }
};

int main() { 
    std::vector<electron> electrons;
double timestep = 0.1;

for(int i = 0; i < 50 ; i++) {
  electrons.push_back(electron());      
}

for(int j = 0; j < 1000; j++) {
    for(int k = 0; k < electrons.size(); k++) {
        electrons[k].updatePosition(timestep);
        electrons[k].updateVelocity(timestep); 
         
            if(electrons[k].position <= -5 || electrons[k].position >= 6) {

                std:: cout << k << " th electron" << std :: endl << "Finished at:  " << j * timestep << std::endl ;
                electrons[k] = electron();
            }
    }

}

}