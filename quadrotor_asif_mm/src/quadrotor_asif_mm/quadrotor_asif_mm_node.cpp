#include "ros/ros.h"
#include <mavros/frame_tf.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/SetMode.h>
#include "asif_controller.hpp"

ros::Publisher quad1_control_pub_;
ros::Publisher quad2_control_pub_;

void actuatorTimerCallback(const ros::TimerEvent &);


mavros_msgs::ActuatorControl quad1_actuator_msg;
mavros_msgs::ActuatorControl quad2_actuator_msg;

void actuatorTimerCallback(const ros::TimerEvent &)
{
    quad1_control_pub_.publish(quad1_actuator_msg);
    quad2_control_pub_.publish(quad2_actuator_msg);
}


int main(int argc, char** argv)
{
    ros::init(argc,argv,"quadrotor_asif_mm_node");
    ros::NodeHandle nh;

    quad1_control_pub_ = nh.advertise<mavros_msgs::ActuatorControl>("/uav0/mavros/actuator_control/asif", 1);
    quad2_control_pub_ = nh.advertise<mavros_msgs::ActuatorControl>("/uav1/mavros/actuator_control/asif", 1);

    ros::Timer actuatorTimer = nh.createTimer(ros::Duration(1.0/250), &actuatorTimerCallback);

    asif::control_s controls[2] {{-0.8, -0.8, -0.8, -0.8},{-0.8, -0.8, -0.8, -0.8}};
    asif::quad_states_s posdes[2];

    ros::ServiceClient arm_quad1_client = nh.serviceClient<mavros_msgs::CommandBool>("/uav0/mavros/cmd/arming");
    ros::ServiceClient set_mode_quad1_client = nh.serviceClient<mavros_msgs::SetMode>("/uav0/mavros/set_mode");

    mavros_msgs::CommandBool arm_quad1;
    mavros_msgs::SetMode set_mode_quad1;

    arm_quad1.request.value = 1;
    set_mode_quad1.request.base_mode = 0;
    set_mode_quad1.request.custom_mode = "OFFBOARD";

    ros::ServiceClient arm_quad2_client = nh.serviceClient<mavros_msgs::CommandBool>("/uav1/mavros/cmd/arming");
    ros::ServiceClient set_mode_quad2_client = nh.serviceClient<mavros_msgs::SetMode>("/uav1/mavros/set_mode");

    mavros_msgs::CommandBool arm_quad2;
    mavros_msgs::SetMode set_mode_quad2;

    arm_quad2.request.value = 1;
    set_mode_quad2.request.base_mode = 0;
    set_mode_quad2.request.custom_mode = "OFFBOARD";

    asif::ASIF_Controller asifController(nh);
    ros::Rate rate(250);

    while (ros::ok() && !asifController.sensors_initialized()) {
        ros::spinOnce();
        rate.sleep();
    }

    if (arm_quad1_client.call(arm_quad1)) {
        ROS_INFO("Arm quadrotor 1 result: %i", arm_quad1.response.result);
    }
    if (arm_quad2_client.call(arm_quad2)) {
        ROS_INFO("Arm quadrotor 2 result: %i", arm_quad2.response.result);
    }

    asifController.update_actuator_quads(controls);
    asifController.getQuadsMsg(quad1_actuator_msg, quad2_actuator_msg);
    ROS_INFO("Initializing actuators...");

    if (set_mode_quad1_client.call(set_mode_quad1)) {
        ROS_INFO("Quadrotor 1 set Mode OFFBOARD: %i", set_mode_quad1.response.mode_sent);
    }
    if (set_mode_quad2_client.call(set_mode_quad2)) {
        ROS_INFO("Quadrotor 2 Set Mode OFFBOARD: %i", set_mode_quad2.response.mode_sent);
    }
    rate.sleep();
    double begin = ros::Time::now().toSec();
    ros::Rate cotrol_rate1(250);
    while (ros::ok()) {
        if ((ros::Time::now().toSec() - begin) < 60.0) {
            posdes[0].x = 0.0;
            posdes[0].y = 0.0;
            posdes[0].z = -2.0;

            posdes[1].x = 0.0;
            posdes[1].y = 1.0;
            posdes[1].z = -1.0;

            asifController.geometric_controller(posdes,controls);
            asifController.update_actuator_quads(controls);
            asifController.getQuadsMsg(quad1_actuator_msg, quad2_actuator_msg);
            cotrol_rate1.sleep();
        } else if ((ros::Time::now().toSec() - begin) < 300.0) {
            posdes[0].x = 0.0;
            posdes[0].y = 0.0;
            posdes[0].z = -2.0;

            posdes[1].x = 0.0;
            posdes[1].y = 0.0;
            posdes[1].z = -2.0;
            asifController.geometric_controller(posdes,controls);

            if (!asifController.asif_filter(controls)) {
              asifController.backup_control(controls);
            }
            asifController.update_actuator_quads(controls);
            asifController.getQuadsMsg(quad1_actuator_msg, quad2_actuator_msg);
            cotrol_rate1.sleep();
        } else {
            asifController.backup_control(controls);
            asifController.update_actuator_quads(controls);
            asifController.getQuadsMsg(quad1_actuator_msg, quad2_actuator_msg);
        }

        ros::spinOnce();
    }

    return 0;
}