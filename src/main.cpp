#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "MPC.h"

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

int main() {

	uWS::Hub h;

	// MPC is initialized here!
	MPC mpc;

	h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
										 uWS::OpCode opCode) {
		// "42" at the start of the message means there's a websocket message event.
		// The 4 signifies a websocket message
		// The 2 signifies a websocket event
		string sdata = string(data).substr(0, length);
		std::cout << sdata << std::endl;
		if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
			string s = hasData(sdata);
			if (s != "") {
				auto j = json::parse(s);
				string event = j[0].get<string>();
				if (event == "telemetry") {
					// j[1] is the data JSON object
					//read next 6 waypoints from the lane detection
					vector<double> ptsx = j[1]["ptsx"];
					vector<double> ptsy = j[1]["ptsy"];

					//read current 6 states
					double px = j[1]["x"];
					double py = j[1]["y"];
					double psi = j[1]["psi"];
					double v = j[1]["speed"];
					double delta = j[1]["steering_angle"];
					double a = j[1]["throttle"];

					// std::cout << "px: " << px << std::endl;
					// std::cout << "py: " << py << std::endl;
					// std::cout << "psi: " << psi << std::endl;
					// std::cout << "v: " << v << std::endl;
					// std::cout << "delta: " << delta << std::endl;
					// std::cout << "a: " << a << std::endl;

					//convert waypoints to vehicle coordinates from global coordinates
					//compute a 3rd order polynomial to model trajectory
					const size_t num_of_waypoints = ptsx.size();

					Eigen::VectorXd waypoints_x_carcoord(num_of_waypoints);
					Eigen::VectorXd waypoints_y_carcoord(num_of_waypoints);

					for (size_t i = 0; i < num_of_waypoints; i++) {
						double diff_x = ptsx[i] - px;
						double diff_y = ptsy[i] - py;
						waypoints_x_carcoord[i] = diff_x * cos(-psi) - diff_y * sin(-psi);
						waypoints_y_carcoord[i] = diff_y * cos(-psi) + diff_x * sin(-psi);
					}

					const int poly_order = 3;
					VectorXd new_coeffs = polyfit(waypoints_x_carcoord, waypoints_y_carcoord, poly_order);
					// for (int i = 0; i < new_coeffs.size(); i++) {
					//   std::cout << "k: " << new_coeffs[i] << ".  ";
					// }
					// std::cout << "\n";

					/**
					 * TODO: Calculate steering angle and throttle using MPC.
					 * Both are in between [-1, 1].
					 */

					//compute first period vehicle states using these readings and model
					Eigen::VectorXd state(mpc.num_of_states);

					// current CTE is fitted polynomial (road curve) evaluated at px = 0.0
					double cte = new_coeffs[0];
					// current heading error is the tangent to the road curve at px = 0.0
					double epsi = -atan(new_coeffs[1]);

					double next_px = 0.0 + v * mpc.dt;
					double next_py = 0.0;
					double next_psi = 0.0 - v * delta / mpc.Lf * mpc.dt;
					double next_v = v + a * mpc.dt;
					double next_cte = cte + v * sin(epsi) * mpc.dt;
					double next_epsi = epsi + v * (-delta) / mpc.Lf * mpc.dt;

					state << next_px, next_py, next_psi, next_v, next_cte, next_epsi;

					std::vector<double> solution = mpc.Solve(state, new_coeffs);
					std::cout << "solution size: " << solution.size() << std::endl; 

					double steer_value;
					double throttle_value;

					steer_value = solution[0] / deg2rad(25);
					throttle_value = solution[1];

					json msgJson;
					// NOTE: Remember to divide by deg2rad(25) before you send the 
					//   steering value back. Otherwise the values will be in between 
					//   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
					msgJson["steering_angle"] = steer_value;
					msgJson["throttle"] = throttle_value;

					// Display the MPC predicted trajectory 
					vector<double> mpc_x_vals;
					vector<double> mpc_y_vals;
					for ( size_t i = 2; i < solution.size(); i++ ) {
						if ( i % 2 == 0 ) {
							mpc_x_vals.push_back( solution[i] );
						} else {
							mpc_y_vals.push_back( solution[i] );
						}
					}

					/**
					 * TODO: add (x,y) points to list here, points are in reference to 
					 *   the vehicle's coordinate system the points in the simulator are 
					 *   connected by a Green line
					 */

					msgJson["mpc_x"] = mpc_x_vals;
					msgJson["mpc_y"] = mpc_y_vals;

					// Display the waypoints/reference line
					vector<double> next_x_vals;
					vector<double> next_y_vals;

					/**
					 * TODO: add (x,y) points to list here, points are in reference to 
					 *   the vehicle's coordinate system the points in the simulator are 
					 *   connected by a Yellow line
					 */
					const double poly_step = 5.0;
					const int poly_num = 10;

					for (size_t i = 0; i < poly_num; i++) {
						double x_coord = poly_step * i;
						double y_coord = polyeval(new_coeffs, x_coord);
						next_x_vals.push_back(x_coord);
						next_y_vals.push_back(y_coord);
					}

					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;


					auto msg = "42[\"steer\"," + msgJson.dump() + "]";
					std::cout << msg << std::endl;
					// Latency
					// The purpose is to mimic real driving conditions where
					//   the car does actuate the commands instantly.
					//
					// Feel free to play around with this value but should be to drive
					//   around the track with 100ms latency.
					//
					// NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
					std::this_thread::sleep_for(std::chrono::milliseconds(100));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
				}  // end "telemetry" if
			} else {
				// Manual driving
				std::string msg = "42[\"manual\",{}]";
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			}
		}  // end websocket if
	}); // end h.onMessage

	h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
		std::cout << "Connected!!!" << std::endl;
	});

	h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
												 char *message, size_t length) {
		ws.close();
		std::cout << "Disconnected" << std::endl;
	});

	int port = 4567;
	if (h.listen(port)) {
		std::cout << "Listening to port " << port << std::endl;
	} else {
		std::cerr << "Failed to listen to port" << std::endl;
		return -1;
	}
	
	h.run();
}