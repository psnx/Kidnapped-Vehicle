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

	cout << "Initialized" << endl;
	cout << "0 particle: " << particles[0].id << endl;
	
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	double x_ , y_, theta_; // the new predicted values	

	for (int i=0; i<num_particles; i++)
	{
		double x_, y_, theta_; //the predicted values
		if (yaw_rate == 0)
		{
			x_ = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y_ = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			theta_ = particles[i].theta; // no change

		} else {
			 // if not 0 yaw rate
			 x_ = particles[i].x + (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta));
			 y_ = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			 theta_ = particles[i].theta + yaw_rate * delta_t;
		}

		std::normal_distribution<double> N_x(x_, std_pos[0]);
		std::normal_distribution<double> N_y(y_, std_pos[1]);
		std::normal_distribution<double> N_theta(theta_, std_pos[2]);
		
		//state transition
		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);

		//particles[i].associations.clear();		
	}
	cout << "predict done" << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double bestX, bestY;
	double minDist = 100;
	int association=0;

	for (int i = 0; i < observations.size(); i++) //each observation from the particles vintage point.
	{		
		minDist = 100;
		for (int j = 0; j < predicted.size(); j++)
		{
						
			double d = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if  (d < minDist)	
			{
				association = predicted[j].id;
				minDist = d;				
			}				
		}
	observations[i].id=association;
	}
}

std::vector<LandmarkObs> ParticleFilter::transformObservations(Particle p, std::vector<LandmarkObs> observations)
{		
	std::vector<LandmarkObs> trans_observations; //the transformed observations
	LandmarkObs obs, trans_obs;	
		
	for (auto& obs : observations)
	{			
		trans_obs.x = p.x + (obs.x * cos(p.theta) - obs.y * sin(p.theta));
		trans_obs.y = p.y + (obs.x * sin(p.theta) + obs.y * cos(p.theta));			
		trans_observations.push_back(trans_obs);		
	}			
	return trans_observations;
}

std::vector<LandmarkObs> ParticleFilter::predictObservations(Particle p, double range,  Map m)
{
	// TODO: range calculation
	LandmarkObs landmark;
	std::vector<LandmarkObs> relevantLandmarks;	

	for (auto &lm : m.landmark_list)
	{
		if (dist(p.x, p.y, lm.x_f, lm.y_f) <= range) 
		{
			landmark.id = lm.id_i;
			landmark.x = lm.x_f;
			landmark.y = lm.y_f;		
			relevantLandmarks.push_back(landmark);			
		}		
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
	vector<int> associations;
	vector<double> s_x, s_y;
	//weights.clear();
	
	for (auto &p: particles)
	{		
		//transform actual observations to map space
		tr_observations = transformObservations(p, observations);
		//predict observations by extracting relevant map data
		predictedObservations = predictObservations(p, sensor_range, map_landmarks);		
		dataAssociation(predictedObservations, tr_observations);
		
		associations.clear();
		s_x.clear();
		s_y.clear();
		p.weight = 1;
		

		for (auto &t_obs : tr_observations)
		{				
			associations.push_back(t_obs.id);
			s_x.push_back(t_obs.x);
			s_y.push_back(t_obs.y);									

			double measX = (t_obs.x);
			double measY= (t_obs.y);

			// numbering of rows in map starts with 1, hence association -1
			double muX = map_landmarks.landmark_list[t_obs.id-1].x_f;
			double muY = map_landmarks.landmark_list[t_obs.id-1].y_f;			
			long double mul = 1 / (2.*M_PI*std_landmark[0]*M_PI*std_landmark[1])*			
				exp(-1*						
				(
					pow(measX - muX , 2.0)/(2*pow(std_landmark[0], 2.0))+
					pow(measY - muY , 2.0)/(2*pow(std_landmark[1], 2.0))
				) );
				p.weight *= mul;
			
			//weights.push_back(p.weight);	
			/*
			cout << "particle id: " << p.id;			
			cout << "\t mul : " << mul;	
			cout << "\t measX: " << measX;	
			cout << "\t muX: " << muX << endl;
			*/		
		}	
		p = SetAssociations(p, associations, s_x, s_y);
		weights.clear();
		for (auto& p : particles){		
		weights.push_back(p.weight);
		}
		
							

	}	
	/*diagnostic
	cout << "after update weights" << endl;
	for (auto& obs : tr_observations)
	{
		cout << "Transformed observations...";
		cout << " \t id: " << obs.id << "\t";	
		cout << " \t x: " << obs.x; 
		cout << " \t y: " << obs.y << endl;
		
		cout << "Coordinated from map......";
		cout << "\t id: " << map_landmarks.landmark_list[obs.id-1].id_i << "\t";
		cout << "\t x: " << map_landmarks.landmark_list[obs.id-1].x_f;
		cout << "\t y: " << map_landmarks.landmark_list[obs.id-1].y_f << endl;
		cout << "-------------------------------------------" << endl;
	}
	cout <<"n ass: " << tr_observations.size() << endl;
	*/

	/*weight print
	for (auto& p : particles){
		cout << "particle: \t" << p.id << " w: \t" << p.weight << endl;
	}
*/
	
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
