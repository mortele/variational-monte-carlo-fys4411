#pragma once

class Particle {
public:
    Particle();
    void setPosition(double* position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    double* getPosition() { return this->position; }

private:
    int     numberOfDimensions;
    double* position;
};

