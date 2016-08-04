#pragma once

class NeuralNetwork
{
public:
    explicit NeuralNetwork(){};
    virtual void Train();
    template<typename T>
    virtual T Simulate() = 0;
};