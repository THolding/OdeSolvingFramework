#include "runge_kutta_ck_solver.hpp"
#include <iostream>
#include <cmath>
#include "model_definition.hpp"
#include "model_driver.hpp"

//Constructor
RungeKuttaCKSolver::RungeKuttaCKSolver(ModelDriver* _model) : model(_model),
    precision(RungeKuttaCKSolver::defaultPrecision)
{  }

//Performs one 'dumb' step for each equation. Writes to _yOut and _yErr.
void RungeKuttaCKSolver::calc(const std::vector<double>& _yCurrent, const std::vector<double>& _derivs, const double  _timestep, std::vector<double>& _yOut, std::vector<double>& _yErr)
{
    //Changes to variable names
    //y == _yCurrent
    //h == _timestep
    //dydt == _derivs
    //youy == _yOut
    //yErr == yErr

    //Get pointer to model definition.
    const ModelDefinition* modelDefinition = model->get_model_definition();

    //places to store intermediate steps
    std::vector<double> k2(model->get_num_compartments()), k3(model->get_num_compartments()), k4(model->get_num_compartments()), k5(model->get_num_compartments()), k6(model->get_num_compartments()), ytemp(model->get_num_compartments());

    //first step: simply dydt[0:i]???
    //2nd step
    for(unsigned int i=0; i<model->get_num_compartments(); i++)
        ytemp[i] = _yCurrent.at(i)+b21*_timestep*_derivs.at(i); //Calculate temporary point to do next step from.
    modelDefinition->calc_derivatives(ytemp, k2, model->get_params()); //Store next set of derivatives (from temporary point just calculated).

    //3rd step
    for(unsigned int i=0; i<model->get_num_compartments(); i++)
        ytemp[i] = _yCurrent.at(i)+_timestep*(b31*_derivs.at(i)+b32*k2[i]); //Calculate temporary point to do next step from.
    modelDefinition->calc_derivatives(ytemp, k3, model->get_params()); //Store next set of derivatives (from temporary point just calculated).

    //4th step
    for(unsigned int i=0; i<model->get_num_compartments(); i++)
        ytemp[i] = _yCurrent.at(i)+_timestep*(b41*_derivs.at(i)+b42*k2[i]+b43*k3[i]); //Calculate temporary point to do next step from.
    modelDefinition->calc_derivatives(ytemp, k4, model->get_params()); //Store next set of derivatives (from temporary point just calculated).

    //5th step
    for(unsigned int i=0; i<model->get_num_compartments(); i++)
        ytemp[i] = _yCurrent.at(i)+_timestep*(b51*_derivs.at(i)+b52*k2[i]+b53*k3[i]+b54*k4[i]); //Calculate temporary point to do next step from.
    modelDefinition->calc_derivatives(ytemp, k5, model->get_params()); //Store next set of derivatives (from temporary point just calculated).

    //6th step
    for(unsigned int i=0; i<model->get_num_compartments(); i++)
        ytemp[i] = _yCurrent.at(i)+_timestep*(b61*_derivs.at(i)+b62*k2[i]+b63*k3[i]+b64*k4[i]+b65*k5[i]); //Calculate temporary point to do next step from.
    modelDefinition->calc_derivatives(ytemp, k6, model->get_params()); //Store next set of derivatives (from temporary point just calculated).

    //Now write output.
    //Update output vector (passed by reference).
    for(unsigned int i=0; i<model->get_num_compartments(); i++)
        _yOut.at(i) = _yCurrent.at(i)+_timestep*(c1*_derivs.at(i)+c3*k3[i]+c4*k4[i]+c6*k6[i]);

    //Update error vector (passed by reference).
    for(unsigned int i=0;i<model->get_num_compartments();i++)
        _yErr.at(i) = _timestep*(dc1*_derivs.at(i)+dc3*k3[i]+dc4*k4[i]+dc5*k5[i]+dc6*k6[i]);
}

//Assesses error for the 'dumb' step (given by calc()) and decides whether the step size needs to be altered.
//Modifies _yCurrent, _stepsize and _nextStepsize.
void RungeKuttaCKSolver::step(std::vector<double>& _yCurrent, const std::vector<double>& _derivs, double *_stepsize, double *_nextStepsize/*, const std::vector<double>& _yScale*/)
{
    //std::vector<double> yTemp(model->num_compartments()), yErr(model->num_compartments()); //written to by the calc function. Holds proposed values for y's next step and their corresponding error.

	//double errMax; //Use to hold the maximum error calculated in any of the proposed y values.
	//double tempStepsize; //Current step size to try

	//tempStepsize = *_stepsize; //Initially the current stepsize should be stepsize.

    //TODO: is the below line actually doing anything??? yTemp and yErr values just overwritten later before they're used...
    //calc(_yCurrent, _derivs, *_stepsize, yTemp, yErr); //Perform a trial step, storing proposed y values in yTemp and corresponding errors in yErr.
    std::vector<double> yErr(model->get_num_compartments());
    calc(_yCurrent, _derivs, *_stepsize, _yCurrent, yErr); //Perform fixed step size step.
    *_nextStepsize = *_stepsize;

    /*
	///Start adaptive timestep bit...
	for(;;) //Searching for an adequate timestep (yeilding acceptable error).
	{
	    if (*_stepsize > 0.1) //Maximum step size.
            *_stepsize = 0.1;

        //Try the current timestep size. Store error in yErr.
		calc(_yCurrent, _derivs, *_stepsize, yTemp, yErr);

		//Calculate error.
		errMax= 0.0;
		for(unsigned int i=0; i<model->num_compartments(); i++)
            errMax = std::max(errMax, std::fabs(yErr[i]/(_yScale.at(i)))); //Choose largest (scaled) error, in any of the y variables.
		errMax /= precision; //Scale error by required tolerance.

		//Check error measure
		if(errMax<=1.0) //if the computed maximum possible error is smaller than 1 break the for loop. Stepsize is accepted.
            break;
        else //Stepsize not accepted. Greater accuracy required. Make stepsize smaller.
        {
            tempStepsize = 0.9*(*_stepsize)*pow(errMax,-0.25); //proposed stepsize

            //Check proposed stepsize against other conditions (at least a than factor of 10 change).
            *_stepsize = (*_stepsize>=0.0 ? std::max(tempStepsize,0.1*(*_stepsize)) : std::min(tempStepsize, 0.1*(*_stepsize))); //Accommodate forward and backward movement in time. No more than a factor of ten.
        }
	}


    //Determine size of step to start with next time.
	if(errMax > 1.89E-4) //If the calculated maximum encountered error is greater than a magic number make it smaller
        *_nextStepsize = 0.9*(*_stepsize)*std::pow(errMax,-0.2);
	else
        *_nextStepsize = 5.0*(*_stepsize); //Otherwise, make it bigger!
	///End adaptive timestep bit
	*/
	//for(unsigned int i=0; i<model->num_compartments(); i++) //Now that an adequate timestep has been verified, and new y values found for it, copy the temporary y values across to the output...
    //    _yCurrent.at(i) = yTemp[i];
}

//Drives the stepper, solving between the start and stop times.
void RungeKuttaCKSolver::solve()
{
    const ModelDefinition* modelDefinition = model->get_model_definition();
	double t = model->get_start_time(); //current time
	double dt = model->get_dt(); //current size of time step
	double dtnext; // UNUSED - artefact of adaptive timestep. next time step size (calculated at runtime for adaptive time step algorithm)
	std::vector<double> derivs(model->get_num_compartments()); //Point-gradient of at t or S, I and R etc. (elements of y).
	std::vector<double> yScale(model->get_num_compartments()); //How error is scaled to fit desired accuracy for each y variable.

	std::vector<double> yCurrent = model->get_current_values();

	//************* start integration ***********
	while (true) //Main loop. bool CompartmentalModel::stopConditionMet() now handles stopping. //t < model->getStopTime()) //TODO: will not work with negative timesteps...
	{
	    //calculate derivs
		modelDefinition->calc_derivatives(yCurrent, derivs, model->get_params());

		//calculate error scaling
		/*for (unsigned int i=0; i<model->get_num_compartments(); i++) //For each equation/variable
		{
		    //Calculate scaling factors for each variable
			if(yCurrent.at(i) != 0.0) //if there's something in y[i].
                yScale[i]= fabs(yCurrent.at(i))+fabs(derivs[i]*dt); //set error threshold scaling? Constant fractional errors (except very near y == 0 inwhich case slope dominates).
			else // y[i]==0.0...
                yScale[i]=1;
		}*/

		step(yCurrent, derivs, &dt, &dtnext/*, yScale*/); //Call stepper to compute next output.

		t+=dt; //Increment time.
		//dt=dtnext; //Load next timestep.

		//Hook back into model (and 'commit' changes).
        model->update(t, yCurrent);

        //Check time step isn't below minimum threshold. If it is, print error message and stop!
		/*if (dt<1E-7)
		{
		    std::cout << "dt = " << dt << " -> ";
			std::cout << "step size underflow at t = " << t << "!\n";
			break;
		}*/

		//Check stop condition hasn't been met.
		if (model->stop_condition_met(t, yCurrent))
		{
		    if (t >= model->get_max_time())
                std::cout << "Stop condition met at t= " << t << ". WARNING: May not have converged to equilibrium state...\n";
            break;
		}

	} //End main while loop
}
