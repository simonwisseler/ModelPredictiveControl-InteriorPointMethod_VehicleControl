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
                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];
                    double px = j[1]["x"];
                    double py = j[1]["y"];
                    double psi = j[1]["psi"]; // angle in rad
                    // DEBUG
                    // std::cout << "psi: " << psi << std::endl;
                    double v = j[1]["speed"];
                    
                    // map waypoints (ptsx, ptsy) into car's coordinate system
                    // (car is at origin of coordinate system (x = 0, y = 0, psi = 0))
                    Eigen::VectorXd waypoints_x(6);
                    Eigen::VectorXd waypoints_y(6);
                    
                    for (unsigned int i = 0; i < ptsx.size(); i++) {
                        double dx = ptsx[i] - px;
                        double dy = ptsy[i] - py;
                        // counter-clockwise rotation "matrix" but using minus psi makes rotation clockwise
                        waypoints_x(i) = dx * cos(-psi) - dy * sin(-psi);
                        waypoints_y(i) = dx * sin(-psi) + dy * cos(-psi);
                    }
                    
                    auto coeffs = polyfit(waypoints_x, waypoints_y, 3);
                    // x = 0 as car at origin of its coordinate system
                    double cte = polyeval(coeffs, 0);
                    // at x = 0, first derivative of polynomial simplifies to coeffs[1] as all other summands multiplied by 0
                    // atan returns angle expressed in radians; minus atan as angle of car (0 in its coordinate system) is smaller than it should be
                    double psi_error = -atan(coeffs[1]);
                    
                    Eigen::VectorXd state(6);
                    state << 0, 0, 0, v, cte, psi_error;
                    
                    auto vars = mpc.Solve(state, coeffs);
                    
                    double steer_value = vars[0];
                    double throttle_value = vars[1];
                    
                    json msgJson;
                    // NOTE: Remember to divide by deg2rad(25) before you send the
                    //   steering value back. Otherwise the values will be in between
                    //   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
                    //
                    msgJson["steering_angle"] = steer_value / deg2rad(25);
                    msgJson["throttle"] = throttle_value;
                    
                    // MPC predicted trajectory
                    // points are in reference to vehicle's coordinate system and connected by a green line in the simulator
                    vector<double> mpc_x_vals;
                    vector<double> mpc_y_vals;
                    
                    for (unsigned int i = 2; i < vars.size(); i++) {
                        if (i % 2 == 0) {
                            mpc_x_vals.push_back(vars[i]);
                        } else {
                            mpc_y_vals.push_back(vars[i]);
                        }
                    }
                    
                    msgJson["mpc_x"] = mpc_x_vals;
                    msgJson["mpc_y"] = mpc_y_vals;
                    
                    
                    // waypoints/reference line
                    //points are in reference to vehicle's coordinate system and connected by a yellow line in the simulator
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;
                    
                    for (unsigned int i = 0; i < 6; i++){
                        next_x_vals.push_back(waypoints_x(i));
                        next_y_vals.push_back(waypoints_y(i));
                    }
                    
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;
                    
                    
                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                    std::cout << msg << std::endl;
                    // Latency
                    // The purpose is to mimic real driving conditions where
                    //   the car doesn't actuate the commands instantly.
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

