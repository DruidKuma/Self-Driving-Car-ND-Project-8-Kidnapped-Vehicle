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
#include "helper_functions.h"

using namespace std;

// Set the number of particles. Initialize all particles to first position (based on estimates of 
// x, y, theta and their uncertainties from GPS) and all weights to 1. 
// Add random Gaussian noise to each particle.
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  num_particles = 1000;

  // Gaussian noise
  normal_distribution<double> d_x(0, std[0]);
  normal_distribution<double> d_y(0, std[1]);
  normal_distribution<double> d_theta(0, std[2]);

  // init particles
  for (unsigned int i = 0; i < num_particles; i++) {
    Particle p;
    p.id = i;
    p.x = x + d_x(gen);
    p.y = y + d_y(gen);
    p.theta = theta + d_theta(gen);
    p.weight = 1.0;

    particles.push_back(p);
  }

  is_initialized = true;
}

// Add measurements to each particle and add random Gaussian noise.
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

  	// Gaussian noise
  	normal_distribution<double> d_x(0, std_pos[0]);
  	normal_distribution<double> d_y(0, std_pos[1]);
  	normal_distribution<double> d_theta(0, std_pos[2]);

  	for (auto& p : particles) {

  		// if yaw_rate is about zero
    	if (fabs(yaw_rate) < 0.00001) {  
      		p.x += velocity * delta_t * cos(p.theta);
      		p.y += velocity * delta_t * sin(p.theta);
    	} else {
      		p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
      		p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
      		p.theta += yaw_rate * delta_t;
    	}

    	// add noise
    	p.x += d_x(gen);
    	p.y += d_y(gen);
    	p.theta += d_theta(gen);
  	}

}

// Find the predicted measurement that is closest to each observed measurement and assign the 
// observed measurement to this particular landmark.
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	
	for (auto& o : observations) {
		
		//init min distance and id of nearest landmark
		double min_distance = numeric_limits<double>::max();
		int landmark_id = -1;

		//find nearest landmark
		for (auto& p : predicted) {
			double distance = dist(o.x, o.y, p.x, p.y);

			if (distance < min_distance) {
				min_distance = distance;
				landmark_id = p.id;
			}
		}

		o.id = landmark_id;
	}

}

// Update the weights of each particle using a mult-variate Gaussian distribution.
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	for (auto& p : particles) {

		// get particle's corrdinates and theta
		double p_x = p.x;
		double p_y = p.y;
		double p_theta = p.theta;

		// predicted landmarks in the sensor range of the particle
		vector<LandmarkObs> predictions;
		for (auto& lm : map_landmarks.landmark_list) {
			if (fabs(lm.x_f - p_x) <= sensor_range && fabs(lm.y_f - p_y) <= sensor_range) {
				predictions.push_back(LandmarkObs{ lm.id_i, lm.x_f, lm.y_f });
			}
		}

		// transform observations' coordinates (from vehicle to map)
		vector<LandmarkObs> transformed_observations;
		for (auto& o : observations) {
			double t_x = cos(p_theta) * o.x - sin(p_theta) * o.y + p_x;
			double t_y = sin(p_theta) * o.x + cos(p_theta) * o.y + p_y;
			transformed_observations.push_back(LandmarkObs{ o.id, t_x, t_y });
		}

		// associate predicted and observed measurements
		dataAssociation(predictions, transformed_observations);

		// reset particle's weight
		p.weight = 1.0;

		for (auto& obs : transformed_observations) {

			// get coordinates of associated prediction				
			double pred_x, pred_y;
			for (auto& pred : predictions) {
				if(pred.id == obs.id) {
					pred_x = pred.x;
					pred_y = pred.y;
					break;
				}
			}

			double s_x = std_landmark[0];
			double s_y = std_landmark[1];

			// calculate new weight with multivariate Gaussian and multiply with total weight
			p.weight *= (1/(2*M_PI*s_x*s_y)) * exp( -( pow(pred_x-obs.x,2)/(2*pow(s_x, 2)) + (pow(pred_y-obs.y,2)/(2*pow(s_y, 2))) ) );
		}
	}
}

// Resample particles with replacement with probability proportional to their weight.
void ParticleFilter::resample() {

  	// collect all particle weights
  	vector<double> weights;
  	for (auto& p : particles) {
  		weights.push_back(p.weight);
  	}

  	// generate random starting index for resampling wheel
  	uniform_int_distribution<int> uniintdist(0, num_particles-1);
  	auto index = uniintdist(gen);

  	// get max weight
  	double max_weight = *max_element(weights.begin(), weights.end());

  	// uniform random distribution [0.0, max_weight)
  	uniform_real_distribution<double> unirealdist(0.0, max_weight);

  	// choose new particles
  	double beta = 0.0;
  	vector<Particle> new_particles;
  	for (int i = 0; i < num_particles; i++) {
    	double randomd = unirealdist(gen);
    	beta += randomd * 2.0;

    	while (beta > weights[index]) {
      		beta -= weights[index];
      		index = (index + 1) % num_particles;
    	}
    
    	new_particles.push_back(particles[index]);
  	}

  	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y) {
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  	vector<int> v = best.associations;
  	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
  	vector<double> v = best.sense_x;
  	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
  	vector<double> v = best.sense_y;
  	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}