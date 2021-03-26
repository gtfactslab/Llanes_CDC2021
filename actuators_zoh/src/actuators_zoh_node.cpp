

#include <ros/ros.h>
#include "mavros_msgs/ActuatorControl.h"

ros::Publisher quad1_control_pub_;
ros::Publisher quad2_control_pub_;

ros::Subscriber quad1_control_sub_;
ros::Subscriber quad2_control_sub_;
void actuatorTimerCallback(const ros::TimerEvent &);


mavros_msgs::ActuatorControl quad1_actuator_msg;
mavros_msgs::ActuatorControl quad2_actuator_msg;

void actuatorTimerCallback(const ros::TimerEvent &)
{
    quad1_control_pub_.publish(quad1_actuator_msg);
    quad2_control_pub_.publish(quad2_actuator_msg);
}

void actuator1NewCallback(const mavros_msgs::ActuatorControl::ConstPtr& msg)
{
    quad1_actuator_msg.controls = msg->controls;
}

void actuator2NewCallback(const mavros_msgs::ActuatorControl::ConstPtr& msg)
{
    quad2_actuator_msg.controls = msg->controls;
}

int main(int argc, char** argv)
{
    ros::init(argc,argv,"actuators_zoh_node");
    ros::NodeHandle nh;

    quad1_control_pub_ = nh.advertise<mavros_msgs::ActuatorControl>("/uav0/mavros/actuator_control", 1);
    quad2_control_pub_ = nh.advertise<mavros_msgs::ActuatorControl>("/uav1/mavros/actuator_control", 1);

    quad1_control_sub_ = nh.subscribe<mavros_msgs::ActuatorControl>("/uav0/mavros/actuator_control/asif",1,&actuator1NewCallback);
    quad2_control_sub_ = nh.subscribe<mavros_msgs::ActuatorControl>("/uav1/mavros/actuator_control/asif",1,&actuator2NewCallback);

    ros::Timer actuatorTimer = nh.createTimer(ros::Duration(1.0/250), &actuatorTimerCallback);

    ros::spin();

    return 0;
}
