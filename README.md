# ScalarEccentric

In order to use the FEW with scalar field
We need to replace the original ode.cc.
Then we can get the orbital evolution or gravitational waveforms

To calculate the Fisher information matrix.
We firstly calculate the orbital evolution for given parameters, and then calculate the 
partial derivative of orbital parameters such as p(t), e(t) and so on with the initial parameters.
Then we use the Mathematica code to calculate the partial derivative of detector's signal with the initial parameters.
Finally, we calculate the Fiser information matrix with ipython.
