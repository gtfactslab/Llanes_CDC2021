#pragma once

#include "asif_util.hpp"
#include "workspace.h"
#include <ros/ros.h>
#include "nav_msgs/Odometry.h"
#include <Eigen/Dense>
#include <mavros/frame_tf.h>
#include "mavros_msgs/ActuatorControl.h"
#include "gazebo_msgs/ModelStates.h"
#include "mav_asif_msgs/ASIF.h"
#include "mav_asif_msgs/control.h"

namespace asif {
    static const std::string kDefaultQuad1Name = "uav0";
    static const std::string kDefaultQuad2Name = "uav1";
    static const std::string kDefaultQuadActuatorPubTopic = "/mavros/actuator_control";

    static constexpr float kDefaultControlRate = 100.0F;
    static constexpr float kDefaultFwdTrajTimestep = 0.01;
    static constexpr float kDefaultMass = 1.545F;
    static constexpr float kDefaultGravity = 9.8066;
    static constexpr float kDefaultMomentInertiaX = 0.029125;

    static constexpr float kDefaultGeomSpringSaturation = 0.03F;
    static constexpr float kDefaultGeomSpringConstant = 3.0F;
    static constexpr float kDefaultGeomVelocityDamp = 0.7F;
    static constexpr float kDefaultGeomAngleErrorGain = 2.0F;
    static constexpr float kDefaultGeomRollRateDamp = 0.5F;

    static constexpr float kDefaultAsifSpringSaturation = 0.03F;
    static constexpr float kDefaultAsifSpringConstant = 3.0F;
    static constexpr float kDefaultAsifVelocityDamp = 0.7F;
    static constexpr float kDefaultAsifAngleErrorGain = 3.0F;
    static constexpr float kDefaultAsifRollRateDamp = 1.0F;
    static constexpr float kDefaultAsifAlpha = 100.0F;

    static constexpr float kDefaultSafeDistanceY = 1.0F;
    static constexpr float kDefaultSafeDistanceZ = 1.0F;
    static constexpr float kDefaultQuadLateralMomentArm = 0.22;
    static constexpr float kDefaultQuadLongitudinalMomentArm = 0.13;
    static constexpr float kDefaultRotorMaxThrust = 7.0664;
    static constexpr float kDefaultRotorThrustConstant = 5.84e-6;
    static constexpr float kDefaultRotorMomentConstant = 0.06F;
    static constexpr float kDefaultRotorMaxRotationalVelocity = 1100.0F;
    class ASIF_Controller {
    public:
        explicit ASIF_Controller(ros::NodeHandle& node_handle);

        virtual ~ASIF_Controller();

        void backup_control(asif::control_s controls[2]);
        void geometric_controller(const asif::quad_states_s pos_des[2], asif::control_s controls[2]);
        int asif_filter(asif::control_s controls[2]);
        void update_actuator_quads(asif::control_s controls[2]);
        void getQuadsMsg(mavros_msgs::ActuatorControl &quad1_actuator_msg, mavros_msgs::ActuatorControl &quad2_actuator_msg);
        bool sensors_initialized();
    protected:

    private:
        float gravity_;
        float mass_;
        float moment_inertia_x_;
        float spring_saturation_;
        float spring_constant_;
        float asif_alpha_;
        float k_angle_error_;
        float k_roll_rate_damp_;
        float k_velocity_damp_;
        float safe_distance_y_;
        float safe_distance_z_;
        float quad_lateral_moment_arm_;
        float quad_longitudinal_moment_arm_;
        float quad_max_motor_thrust_;
        float rotor_max_rot_velocity_;
        float rotor_thrust_constant_;
        float rotor_moment_constant_;
        float dt_control_;
        float dt_fwd_traj_;
        int backup_barrier_steps_;
        float geom_k_angle_error_;
        float geom_k_roll_rate_damp_;
        float geom_k_velocity_damp_;
        float geom_spring_saturation_;
        float geom_spring_constant_;

        std::string quad1_name_;
        std::string quad2_name_;

        int quad1_gazebo_model_id_ {};
        int quad2_gazebo_model_id_ {};

        ros::NodeHandle& node_handle_;
        std::string actuator_topic_;
//        ros::Subscriber quad1_state_sub_;
//        ros::Subscriber quad2_state_sub_;
        ros::Subscriber quad_state_sub_;
        ros::Publisher quad1_control_pub_;
        ros::Publisher quad2_control_pub_;
        ros::Publisher quad_asif_pub_;
        ros::Publisher quad_asif_force_moment_pub_;
        ros::Publisher quad_geom_force_moment_pub_;
        mavros_msgs::ActuatorControl quad1_actuator_msg_;
        mavros_msgs::ActuatorControl quad2_actuator_msg_;
        mav_asif_msgs::ASIF quad_asif_msg_;
        mav_asif_msgs::control quad_geom_force_moment_msg_;
        mav_asif_msgs::control quad_asif_force_moment_msg_;
        ros::Timer actuatorTimer;

        bool sensor_data_initialized_ = false;

        OSQPWorkspace workspace_;
        asif::ASIF *asif_p_;
        asif::quad_states_s quads_[2];
//        asif::quad_states_s quads_test_[2];


        int iris0_idx_;
        int iris1_idx_;

        void actuatorTimerCallback(const ros::TimerEvent &);
        float motor_force_to_mixer_input(float force) const;
//        void quad1_state_callback(const nav_msgs::Odometry::ConstPtr& msg);
//        void quad2_state_callback(const nav_msgs::Odometry::ConstPtr& msg);
        void quad_state_callback(const gazebo_msgs::ModelStates::ConstPtr& msg);
        Eigen::Vector3d quaternion_to_rpy_wrap(const Eigen::Quaterniond& q);
        Eigen::Vector3d euler_rates_to_body_rates(const Eigen::Vector3d orientation, const Eigen::Vector3d& euler_rates);
    };
}
