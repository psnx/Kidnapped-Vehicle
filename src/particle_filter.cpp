/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	std::default_random_engine gen;
	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);

	for (int i=0; i<num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(particle.weight);
	}
	is_initialized = true;

	cout<< "Initialized" << endl;

	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine generator;
	for (int i=0; i<num_particles; i++)
	{
		double x_, y_, theta_; //the predicted values
		if (yaw_rate=0)
		{
			x_ = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y_ = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			theta_ = particles[i].theta; // no change

		} else {
			 // if not 0 yaw rate
			 x_ = particles[i].x + (velocity/yaw_rate) * (sin(particles[i].theta+yaw_rate*delta_t)-sin(theta_));
			 y_ = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta)-cos(theta_ + yaw_rate*delta_t));
			 theta_ = particles[i].theta + yaw_rate * delta_t;
		}

		std::normal_distribution<double> N_x(x_, std_pos[0]);
		std::normal_distribution<double> N_y(y_, std_pos[1]);
		std::normal_distribution<double> N_theta(theta_, std_pos[2]);
		
		particles[i].x = N_x(generator);
		particles[i].y = N_y(generator);
		particles[i].theta = N_theta(generator);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double bestX, bestY;
	double minDist = 1000;	
	int association;

	for (int i = 0; i < observations.size(); i++) //each observation from the particles vintage point.
	{
		LandmarkObs observed = observations[i];
		for (int j = 0; j < predicted.size(); j++)
		{			
			double d = dist(predicted[j].x, predicted[j].y, observed.x, observed.y);
			if  (d < minDist)	
			{
				association = predicted[j].id;
				minDist = d;
			}
			observations[i].id = association;
		}
	}
}

std::vector<LandmarkObs> ParticleFilter::transformObservations(Particle p, std::vector<LandmarkObs> observations)
{		
	std::vector<LandmarkObs> trans_observations; //the transformed observations
	LandmarkObs obs, trans_obs;	
		
	for (int j = 0; j < observations.size(); j++)
	{	
		obs = observations[j];
		trans_obs.x = p.x + obs.x*cos(p.theta) - obs.y*sin(p.theta);
		trans_obs.y = p.y + obs.x*sin(p.theta) + obs.y*cos(p.theta);			
		trans_observations.push_back(trans_obs);		
	}			
	return trans_observations;
}

std::vector<LandmarkObs> ParticleFilter::predictObservations(Map m)
{
	// Range calculation not implemented due to small size of map
	LandmarkObs landmark;
	std::vector<LandmarkObs> relevantLandmarks;

	for (int i = 0; i < m.landmark_list.size(); i++){
		
		landmark.id = m.landmark_list[i].id_i;
		landmark.x = m.landmark_list[i].x_f;
		landmark.y = m.landmark_list[i].y_f;		
		relevantLandmarks.push_back(landmark);
	}
	return relevantLandmarks;
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) 
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_nor-mal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	vector<LandmarkObs> tr_observations;
	vector<LandmarkObs> predictedObservations;
	
	for (int i = 0; i <num_particles; i++)
	{
		Particle p = particles[i];
		//transform actual observations to map space
		tr_observations = transformObservations(p, observations);
		//predict observations by extracting relevant map data
		predictedObservations = predictObservations(map_landmarks);		
		dataAssociation(predictedObservations, tr_observations);

		for (int j=0; j < tr_observations.size(); j++)
		{
			p.associations.push_back(tr_observations[j].id);
			p.sense_x.push_back(tr_observations[j].x);
			p.sense_y.push_back(tr_observations[j].y);

			double measX = (tr_observations[j].x);
			double measY= (tr_observations[j].y);

			double muX = map_landmarks.landmark_list[tr_observations[j].id].x_f;
			double muY = map_landmarks.landmark_list[tr_observations[j].id].y_f;			

			long double multiplier = 1 / (2*M_PI*std_landmark[0]*M_PI*std_landmark[1])*
			exp(-						
			(
				pow(measX - muX , 2.0)/(2*pow(std_landmark[0], 2.0)) + 
				pow(measY - muY , 2.0)/(2*pow(std_landmark[1], 2.0))
			) );
			p.weight*=multiplier;		

			
		}

		

	}	
}




void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	default_random_engine generator;
	discrete_distribution<int> distr(weights.begin(), weights.end());
	vector<Particle> resampled;

	for (int i = 0; i<num_particles; i++)
	{
		resampled.push_back(particles[distr(generator)]);
	}
	particles = resampled;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
